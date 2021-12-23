args <- commandArgs(TRUE)
#input
isocan <- args[1]# "~/long_read/isocan"
genome_pkg <- args[2]
chimera_paf_binded <- args[3]
input_files <- args[4]
sv_df_colname <- args[5]
ref_gtf <- args[6]
sqanti_junc <- args[7]
#output
chimera_output <- args[8] # .RDS
gene_fragment_gtf <- args[9]
#parameter
sv_ambiguity_thres <- as.integer(args[10]) # 100L
sv_internal_exonend_distance_thres <- as.integer(args[11]) #100000L
n_worker <- as.integer(args[12])
canonical_list <- list(c("GT", "AG"), c("GC", "AG"), c("AT", "AC"))

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(rlang)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(furrr)
  pkgload::load_all(isocan, export_all = FALSE)
  pkgload::load_all(genome_pkg, export_all = FALSE)
  pkgname <- fs::path_file(genome_pkg)
  genome <- eval(parse(text = paste0(pkgname, "::", pkgname)))
})
plan(multicore, workers = n_worker)



# input ----
candi_df <- read_tsv(chimera_paf_binded, col_types = "ciiicciicc")
input_files <- read_csv(input_files, col_types = cols(.default = col_character()))

# functions ----

rc <- function(x) {
  x[!is.na(x)] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x[!is.na(x)])))
  x
}

na2false <- function(x) `[<-`(x, is.na(x), FALSE)

get_seq_2_arrow_na <- function(chr, strand, start, end, ref_genome) {
  res_seq <- rep(NA_character_, length(chr))
  gtf_gr <- 
    tibble(seqnames = chr, start = start, end = end, strand = strand, .name_repair = "minimal") %>% 
    filter(start <= end) %>% 
    GRanges()

  gtf_gr <- split(gtf_gr, seq_len(length(gtf_gr)))
  tx_seq <- as.character(extractTranscriptSeqs(ref_genome, gtf_gr))
  res_seq[na2false(start <= end)] <- tx_seq
  res_seq
}

is_canonical <- function(donor_seq, acceptor_seq, canonical = canonical_list) {
  donor_SEQ <- str_to_upper(donor_seq)
  acceptor_SEQ <- str_to_upper(acceptor_seq)
  map(canonical, ~ donor_SEQ == .x[1] & acceptor_SEQ == .x[2]) %>% pmap_lgl(any)
}

combine_seqid_idx <- function(df)  mutate(df, seq_id = str_c(seq_id, idx, sep = "___"), idx = NULL)
separate_seqid_idx <- function(df) separate(df, seq_id, c("seq_id", "idx"), sep = "___") %>% mutate_at(vars("idx"), as.integer)

retrieve_sv_info <- function(candi_df_mod, retain_colnames = NULL, add_sv_number = FALSE) {
  add_sv_number_fun <- function(x, add_sv_number) {
    if (add_sv_number) {
      x %>% 
        group_by(seq_id) %>% 
        mutate(sv_number = row_number()) %>% 
        ungroup()
    } else {
      x
    }
  }
  
  req_cols <- c("seq_id", "idx", "chr", "strand", "start", "end")
  stopifnot(all(retain_colnames %in% colnames(candi_df_mod)))
  retain_colnames <- retain_colnames[!retain_colnames %in% req_cols]
  
  candi_df_mod_tmp <- 
    candi_df_mod %>% 
    mutate(up_pos = if_else(strand == "+", end, start),
           down_pos = if_else(strand == "+", start, end)) %>% 
    group_by(seq_id) %>% 
    arrange(seq_id, idx) %>% 
    mutate(down_pos = lead(down_pos),
           down_chr = lead(chr),
           down_strand = lead(strand)) %>% 
    ungroup()
  
  if (is.null(retain_colnames)) {
    candi_df_mod_tmp %>% 
      select(seq_id, sample, up_chr = "chr", up_strand = "strand", up_pos, down_chr, down_strand, down_pos) %>% 
      drop_na() %>% 
      add_sv_number_fun(add_sv_number)
  } else {
    candi_df_mod_tmp <- mutate(candi_df_mod_tmp, .internal_idx = row_number())
    left_join(
      candi_df_mod_tmp %>% 
        select(.internal_idx, seq_id, sample, up_chr = "chr", up_strand = "strand", up_pos, down_chr, down_strand, down_pos) %>% 
        drop_na() %>% 
        add_sv_number_fun(add_sv_number),
      candi_df_mod_tmp %>% 
        select(.internal_idx, !!!syms(retain_colnames)),
      by = ".internal_idx"
    ) %>% 
      select(-.internal_idx)
  }
}

determine_right_junction_mono <- function(overlap_length_mono, donor_candi_seq_mono, acceptor_candi_seq_mono) {
  if (is.na(overlap_length_mono) || is.na(donor_candi_seq_mono) || is.na(acceptor_candi_seq_mono)) {
    NA_integer_
  } else {
    donor_seq_vec <- str_sub(donor_candi_seq_mono, 0:overlap_length_mono + 1L, 0:overlap_length_mono + 2L)
    acceptor_seq_vec <- str_sub(acceptor_candi_seq_mono, 0:overlap_length_mono + 1L, 0:overlap_length_mono + 2L)
    which(is_canonical(donor_seq_vec, acceptor_seq_vec))
  }
}

determine_right_junction <- function(overlap_length, donor_candi_seq, acceptor_candi_seq, borbose = TRUE) {
  candi_len <- str_length(donor_candi_seq)
  stopifnot(all(candi_len - overlap_length == 2L | is.na(candi_len) | is.na(overlap_length)))
  stopifnot(all(candi_len == str_length(acceptor_candi_seq) | is.na(acceptor_candi_seq)))
  future_pmap(list(overlap_length, donor_candi_seq, acceptor_candi_seq), determine_right_junction_mono, .progress = borbose)
}


modify_candidate_df <- function(candi_df) {
  # preparation
  seq_id_multi <- 
    candi_df %>% 
    count(seq_id) %>% 
    filter(n > 1L) %>% 
    `[[`("seq_id")
  candi_df %>% 
    filter(seq_id %in% seq_id_multi) %>% 
    group_by(seq_id) %>% 
    arrange(query_start) %>% 
    mutate(idx = row_number(),
           acceptor_overlap = lag(query_end) - query_start,
           donor_overlap = lead(acceptor_overlap)) %>% 
    ungroup()
}


modify_sv_junction <- function(candi_df_mod, borbose = TRUE) {
  inform("obtaining sequences...")
  # plus strand
  candi_df_mod_plus <- 
    candi_df_mod %>% 
    filter(strand == "+") %>% 
    mutate(acceptor_seq_candi = get_seq_2_arrow_na(chr, strand, start - 1L, start + acceptor_overlap, genome), #0-based
           donor_seq_candi = get_seq_2_arrow_na(chr, strand, end + 1L - donor_overlap, end + 2L, genome)
    )
  
  # minus strand
  candi_df_mod_minus <- 
    candi_df_mod %>% 
    filter(strand == "-") %>% 
    mutate(acceptor_seq_candi = get_seq_2_arrow_na(chr, strand, end + 1L - acceptor_overlap, end + 2L, genome), #0-based
           donor_seq_candi = get_seq_2_arrow_na(chr, strand, start - 1L, start + donor_overlap, genome)
    )
  
  inform("determining right SV junctions...")
  candi_df_right <- 
    bind_rows(
      candi_df_mod_plus,
      candi_df_mod_minus
    ) %>% 
    group_by(seq_id) %>% 
    arrange(query_start) %>% 
    mutate(donor_right_junction = determine_right_junction(donor_overlap, donor_seq_candi, lead(acceptor_seq_candi), borbose = borbose))
  
  inform("removing SV-undetermined chimeras...")
  candi_df_count <- 
    candi_df_right %>% 
    ungroup() %>% 
    mutate(n_candi = map_int(donor_right_junction, length))
  candi_df_dropout <- 
    candi_df_count %>% 
    filter(n_candi != 1L | donor_overlap < 0L) %>% 
    select(seq_id, n_candi)
  candi_df_selected <- 
    candi_df_count %>% 
    filter(n_candi == 1L & donor_overlap >= 0L) %>% 
    select(seq_id, n_candi)
  
  seq_id_partial_determined <- intersect(candi_df_selected$seq_id, candi_df_dropout$seq_id) %>% unique()
  seq_id_no_candi <- candi_df_dropout$seq_id[candi_df_dropout$n_candi <= 1L] %>% setdiff(seq_id_partial_determined)
  seq_id_multi_candi <- candi_df_dropout$seq_id[candi_df_dropout$n_candi > 1L] %>% setdiff(seq_id_partial_determined)
  
  inform("calculating genomic positions...")
  modify_paf_sv <- function(candi_df_right) {
    undetermined_rows <- !map_lgl(candi_df_right$donor_right_junction, ~ !identical(.x, integer()) && !is.na(.x) && length(.x) == 1L)
    candi_df_right$donor_right_junction[undetermined_rows] <- NA_integer_
    candi_df_right_mod <- 
      candi_df_right %>% 
      unnest(cols = c(donor_right_junction)) %>% 
      mutate(acceptor_right_junction = lag(donor_right_junction)) %>% 
      ungroup() %>% 
      mutate(query_start = if_else(is.na(acceptor_right_junction), query_start, query_start + acceptor_right_junction - 1L),
             query_end = if_else(is.na(donor_right_junction), query_end, query_end + donor_right_junction - donor_overlap - 1L),
             cs = NA_character_,
             sv_type = if_else(is.na(donor_right_junction), NA_character_, "sj_confirmed")
      )
    paf <- 
      bind_rows(
        candi_df_right_mod %>% 
          filter(strand == "+") %>% 
          mutate(start = if_else(is.na(acceptor_right_junction), start, start + acceptor_right_junction - 1L),
                 end = if_else(is.na(donor_right_junction), end, end + donor_right_junction - donor_overlap - 1L)
          ),
        candi_df_right_mod %>% 
          filter(strand == "-") %>% 
          mutate(end = if_else(is.na(acceptor_right_junction), end, end - acceptor_right_junction + 1L),
                 start = if_else(is.na(donor_right_junction), start, start - donor_right_junction + donor_overlap + 1L)
          )
      ) %>% 
      mutate_at(vars("start", "end"), as.integer) %>% 
      select(seq_id, query_length, query_start, query_end, strand, chr, start, end, cs, sample, idx, sv_type) %>% 
      arrange(sample, seq_id, idx)
    
    sv <- retrieve_sv_info(paf, "sv_type", add_sv_number = TRUE)
    
    list(paf %>% select(-sv_type), sv)
  }
  # partial_determined
  partial_determined_res <- 
    candi_df_right %>% 
    filter(seq_id %in% seq_id_partial_determined) %>% 
    modify_paf_sv()
  
  # full_determined
  full_determined_res <- 
    candi_df_right %>% 
    filter(!seq_id %in% candi_df_dropout$seq_id) %>% 
    modify_paf_sv()
  
  inform("done.")
  list(
    paf = full_determined_res[[1]],
    sv = full_determined_res[[2]],
    partial_paf = partial_determined_res[[1]],
    partial_sv = partial_determined_res[[2]] %>% replace_na(list(sv_type = "sj_not_confirmed")),
    no_candidate_id = seq_id_no_candi,
    multi_candidate_id = seq_id_multi_candi
  )
}

judge_sv_match <- function(chimera_sv_candi, mode = c("ambiguity", "internal"), sv_df, sv_ambiguity_thres = 100L, sv_internal_exonend_distance_thres = 100000L) {
  mode <- match.arg(mode)
  by_1 <- c("sample" = "Sample", "up_chr" = "Chr_1", "up_strand" = "Dir_1", "down_chr" = "Chr_2", "down_strand" = "Dir_2")
  by_2 <- c("sample" = "Sample", "up_chr" = "Chr_2", "up_strand" = "Dir_2", "down_chr" = "Chr_1", "down_strand" = "Dir_1")
  
  chimera_sv_joined <- 
    future_map2_dfr(
      list(
        sv_df %>% mutate(Dir_2 = if_else(Dir_2 == "+", "-", "+")) %>% dplyr::rename(Pos_up = "Pos_1", Pos_down = "Pos_2"),
        sv_df %>% mutate(Dir_1 = if_else(Dir_1 == "+", "-", "+"), Inserted_Seq = rc(Inserted_Seq)) %>% dplyr::rename(Pos_up = "Pos_1", Pos_down = "Pos_2"),
        sv_df %>% mutate(Dir_2 = if_else(Dir_2 == "+", "-", "+")) %>% dplyr::rename(Pos_up = "Pos_2", Pos_down = "Pos_1"),
        sv_df %>% mutate(Dir_1 = if_else(Dir_1 == "+", "-", "+"), Inserted_Seq = rc(Inserted_Seq)) %>% dplyr::rename(Pos_up = "Pos_2", Pos_down = "Pos_1")
      ),
      list(by_1, by_1, by_2, by_2),
      ~ inner_join(chimera_sv_candi, .x, by = .y)
    )
  
  if (mode == "ambiguity") {
    chimera_sv_match <- 
      chimera_sv_joined %>% 
      filter(abs(up_pos - Pos_up) <= sv_ambiguity_thres & abs(down_pos - Pos_down) <= sv_ambiguity_thres)
  } else {
    minus2negative <- function(strand) if_else(strand == "+", 1L, -1L)
    chimera_sv_match <- 
      chimera_sv_joined %>% 
      mutate(.up_distance = minus2negative(up_strand) * (Pos_up - up_pos),
             .down_distance = minus2negative(down_strand) * (down_pos - Pos_down)) %>% 
      filter(-sv_ambiguity_thres <= .up_distance &
               .up_distance <= sv_internal_exonend_distance_thres &
                -sv_ambiguity_thres <= .down_distance & 
               .down_distance <= sv_internal_exonend_distance_thres) %>% 
      select(-c(.up_distance, .down_distance))
  }
  
  sv_unmatch_id <- 
    anti_join(
      chimera_sv_candi,
      chimera_sv_match,
      by = colnames(chimera_sv_candi)
    ) %>% 
    `[[`("seq_id") %>% 
    unique()
  
  sv_partial_match_id <- sv_unmatch_id[sv_unmatch_id %in% chimera_sv_match$seq_id]
  sv_unmatch_id <- setdiff(sv_unmatch_id, sv_partial_match_id)
  sv_unmatch <- chimera_sv_candi %>% filter(seq_id %in% sv_unmatch_id)
  sv_partial_match <- chimera_sv_match %>% filter(seq_id %in% sv_partial_match_id)
  sv_partial_match <- 
    chimera_sv_candi %>% 
    filter(seq_id %in% sv_partial_match_id) %>% 
    left_join(sv_partial_match, by = colnames(chimera_sv_candi))
  
  sv_match <- chimera_sv_match %>% filter(!seq_id %in% sv_partial_match_id)
  
  list(
    sv_match_id = unique(sv_match$seq_id),
    sv_match = sv_match,
    sv_unmatch_id = sv_unmatch_id,
    sv_unmatch = sv_unmatch,
    sv_partial_match_id = sv_partial_match_id,
    sv_partial_match = sv_partial_match
  )
}

convert_paf_variant_multisample <- function(paf) {
  if (nrow(paf) == 0L) {
    dammy_res <- 
      tibble(
        seq_id = character(), chr = character(), type = character(),
        start = integer(), end = integer(), ref = character(), alt = character(), tx_strand = character(),
        sample = character(), tx_start = integer(), tx_end = integer(),
        .name_repair = "minimal"
      )
    return(dammy_res)
  }
  isocan::convert_paf_variant(paf, "") %>% 
    select(-sample) %>% 
    inner_join(paf %>% select(seq_id, sample), by = "seq_id")
}

retrieve_gtf_multisample <- function(paf, retreive_type = "exon", paf_variants = NULL) {
  if (nrow(paf) == 0L) {
    dammy_res <- 
      tibble(
        seq_id = character(), strand = character(), chr = character(), start = integer(),
        end = integer(), sample = character(), type = character(),
        .name_repair = "minimal"
      )
    return(dammy_res)
  }
  samples <- unique(paf$sample)
  map2_dfr(
    map(samples, ~ paf %>% filter(sample == .x)), 
    samples,
    isocan::retrieve_gtf,
    retreive_type = retreive_type,
    paf_variants = paf_variants
  )
}

get_id_noncanonical_notinref_notinsqanti <- function(seqid_idx_compressed_intron, ref_intron, sqanti_intron_additional, borbose = TRUE) {
  intron_no_ref <- 
    seqid_idx_compressed_intron %>% 
    anti_join(ref_intron, by = c("chr", "strand", "start", "end")) %>% 
    anti_join(sqanti_intron_additional, by = c("chr", "strand", "start", "end")) %>% 
    mutate(seq_id = str_remove(seq_id, "___[:digit:]+$"))
  
  intron_no_ref %>% 
    mutate(
      donor_seq =  get_seq_2_arrow_na(chr,strand,  start, start + 1L, genome),
      acceptor_seq = get_seq_2_arrow_na(chr, strand, end - 1L, end, genome)
    ) %>% 
    filter(!is_canonical(donor_seq, acceptor_seq)) %>% 
    `[[`("seq_id") %>% 
    unique()
}

add_exon_number <- function(chimera_gtf) {
  chimera_gtf %>% 
    mutate(start_mod = if_else(strand == "+", start, -start)) %>% 
    group_by(seq_id) %>% 
    arrange(seq_id, idx, start_mod) %>% 
    mutate(exon_number = row_number()) %>% 
    ungroup() %>% 
    select(-start_mod)
}

## main ----
candi_df_mod <- modify_candidate_df(candi_df)
modified_res <- modify_sv_junction(candi_df_mod) # only retain SV determined chimera in main results
modified_res_paf_2 <- modified_res[["paf"]] %>% combine_seqid_idx()
modified_res_partial_paf_2 <- modified_res[["partial_paf"]] %>% combine_seqid_idx()

candi_variant_list <- 
  map2(
    list(modified_res$sv$seq_id, modified_res$partial_sv$seq_id),
    list(modified_res_paf_2, modified_res_partial_paf_2),
    ~ {
      candi_df_mod %>% 
        filter(seq_id %in% .x) %>% 
        combine_seqid_idx() %>% 
        convert_paf_variant_multisample() %>% 
        select(-c(tx_start, tx_end)) %>% # modify tx_start/tx_end
        inner_join(.y %>% select(seq_id, tx_start = "start", tx_end = "end"),
                   by = "seq_id")
    }
  ) %>% 
  `names<-`(c("full", "partial"))

sj_overlap_id_list <- # additional dropout
  map(
    candi_variant_list, 
    ~ {
      .x %>% 
        filter(!(tx_start <= start & end <= tx_end) & type == "~") %>% 
        `$`("seq_id") %>% 
        str_remove("___[:digit:]+$")
    }
  ) %>% 
  `names<-`(c("full", "partial"))

# update
candi_variant_list <- 
  map2(
    candi_variant_list, sj_overlap_id_list, 
    ~ filter(..1, tx_start <= start & end <= tx_end & (!str_remove(seq_id, "___[:digit:]+$") %in% ..2))
  )

modified_res_paf_2 <-
  modified_res_paf_2 %>% 
  filter(!str_remove(seq_id, "___[:digit:]+$") %in% sj_overlap_id_list$full)
modified_res_partial_paf_2 <-
  modified_res_partial_paf_2 %>% 
  filter(!str_remove(seq_id, "___[:digit:]+$") %in% sj_overlap_id_list$partial)

candi_gtf_list <- 
  pmap(
    list(
      list(modified_res_paf_2, modified_res_paf_2, modified_res_partial_paf_2, modified_res_partial_paf_2),
      rep(c("exon", "intron"), 2L),
      list(candi_variant_list$full, candi_variant_list$full, candi_variant_list$partial, candi_variant_list$partial)
    ),
    retrieve_gtf_multisample
  ) %>% 
  `names<-`(c("exon_full", "intron_full", "exon_partial", "intron_partial"))

## prepare ref SJ ----
ref_gtf <- isocan::read_gtf(ref_gtf)
sqanti_junc <- 
  read_tsv(
    sqanti_junc, 
    col_types = cols(
      chrom = col_character(),
      strand = col_character(),
      genomic_start_coord = col_integer(),
      genomic_end_coord = col_integer(),
      .default = col_skip()
    )) %>% 
  rename(chr = "chrom", start = "genomic_start_coord", end = "genomic_end_coord")

ref_intron <- 
  ref_gtf %>% 
  filter(feature == "exon") %>% 
  group_by(transcript_id) %>% 
  mutate(intron_start = end + 1L,
         intron_end = lead(start) - 1L) %>% 
  ungroup() %>% 
  select(chr = "seqname", strand, start = "intron_start", end = "intron_end") %>% 
  drop_na() %>% 
  distinct()

sqanti_intron_additional <- 
  sqanti_junc %>% 
  distinct() %>% 
  anti_join(ref_intron, by = c("chr", "strand", "start", "end"))

## clean up ----
clean_up_obj <- c("candi_df", "ref_gtf", "sqanti_junc")

rm(list = clean_up_obj)

dammy_sv_df <- 
  tibble(
    Chr_1 = character(),
    Pos_1 = integer(),
    Dir_1 = character(),
    Chr_2 = character(),
    Pos_2 = integer(),
    Dir_2 = character(),
    Inserted_Seq = character(),
    Variant_Type = character(),
    Tumor_VAF = double(),
    Sample = character(),
    .name_repair = "minimal"
  )
read_genomon_sv <- function(path, sample_name) {
  if (is.na(path)) {
    dammy_sv_df
  } else {
    readr::read_tsv(
      path, 
      comment = "#", 
      col_types = cols(
        Chr_1 = col_character(),
        Pos_1 = col_integer(),
        Dir_1 = col_character(),
        Chr_2 = col_character(),
        Pos_2 = col_integer(),
        Dir_2 = col_character(),
        Inserted_Seq = col_character(),
        Variant_Type = col_character(),
        Tumor_VAF = col_double(),
        .default = col_skip()
      )
    ) %>% 
      mutate(Sample = sample_name)
  }
}

sv_df <- 
  map2_dfr(input_files[[sv_df_colname]], input_files$sample, read_genomon_sv) %>% 
  rename(Type = "Variant_Type") %>% 
  mutate(Inserted_Seq = str_remove(Inserted_Seq, "^---$"))

## main 2 ----
improper_sj_id_list <- 
  map(
    candi_gtf_list[c("intron_full", "intron_partial")],
    get_id_noncanonical_notinref_notinsqanti, # additional dropout
    ref_intron,
    sqanti_intron_additional
  )

chimera_res_gtf_list <- 
  map2(
    candi_gtf_list[c("exon_full", "exon_partial")],
    improper_sj_id_list,
    ~ {
      .x %>% 
        separate_seqid_idx() %>% 
        filter(!seq_id %in% .y) %>% # update
        add_exon_number()
    }
  ) %>% 
  `names<-`(c("full", "partial"))

chimera_res_variant_list <- 
  map2(
    candi_variant_list,
    improper_sj_id_list,
    ~ .x %>% separate_seqid_idx() %>% filter(!seq_id %in% .y) # update
  ) %>% 
  `names<-`(c("full", "partial"))

chimera_res_sv_list <- 
  map2(
    c("sv", "partial_sv"),
    list(chimera_res_gtf_list$full$seq_id, chimera_res_gtf_list$partial$seq_id),
    ~ modified_res[[.x]] %>% filter(seq_id %in% .y)
  ) %>% 
  `names<-`(c("full", "partial"))

# SV match
chimera_res_sv_match_list <- 
  map2(
    list(
      chimera_res_sv_list$full,
      chimera_res_sv_list$partial %>% filter(sv_type == "sj_not_confirmed"),
      chimera_res_sv_list$partial %>% filter(sv_type == "sj_confirmed")
    ),
    c("internal", "ambiguity", "internal"),
    judge_sv_match,
    sv_df, sv_ambiguity_thres, sv_internal_exonend_distance_thres
  ) %>% 
  `names<-`(c("full", "amb_partial", "int_partial"))

replace_df_mod <- function(df, x, y, replacement) if (nrow(df) == 0) df else `[<-`(df, x, y, replacement)

# IDs ----
sv_multi_match_id <-
  chimera_res_sv_match_list$full$sv_match %>% 
  distinct(seq_id, sample, up_chr, up_strand, up_pos, down_chr, down_strand, down_pos) %>% 
  count(seq_id) %>% 
  filter(n > 1L) %>% 
  `$`("seq_id")
sv_uni_match_id <- 
  setdiff(chimera_res_sv_match_list$full$sv_match_id, sv_multi_match_id)
sv_partial_match_id <-
  chimera_res_sv_match_list$full$sv_partial_match$seq_id %>% unique()

sv_multi_match_add_id <- 
  intersect(
    chimera_res_sv_match_list$amb_partial$sv_match_id,
    chimera_res_sv_match_list$int_partial$sv_match_id
  )
sv_partial_match_add_id <-
  list(chimera_res_sv_match_list$amb_partial$sv_match_id,
       chimera_res_sv_match_list$amb_partial$sv_partial_match_id,
       chimera_res_sv_match_list$int_partial$sv_match_id,
       chimera_res_sv_match_list$int_partial$sv_partial_match_id) %>% 
  purrr::reduce(union) %>% 
  unique() %>% 
  setdiff(sv_multi_match_add_id)

dropout_chimera_id <-
  list(
    no_sj_candidate = modified_res$no_candidate_id, # No canonical SV position candidates
    multi_sj_candidate = modified_res$multi_candidate_id, # Multi canonical SV position candidates
    sv_sj_overlap = c(sj_overlap_id_list$full, sj_overlap_id_list$partial), # SV overlaps with SJ
    improper_sj = improper_sj_id_list$full, # some non-SV SJ are not canonical or not in ref or sqanti-filtered PB.
    sv_unmatch = c(chimera_res_sv_match_list$full$sv_unmatch_id, # no SV matches with given SV list
                   intersect(chimera_res_sv_match_list$amb_partial$sv_unmatch_id,
                             chimera_res_sv_match_list$int_partial$sv_unmatch_id))
  )

# salvage cross-SV chimeras ----
chimera_sv_candi <- 
  candi_df_mod %>% 
  filter(seq_id %in% c(dropout_chimera_id$no_sj_candidate, dropout_chimera_id$multi_sj_candidate)) %>% 
  retrieve_sv_info(add_sv_number = TRUE)

sv_match_salvage <- judge_sv_match(chimera_sv_candi, "ambiguity", sv_df, sv_ambiguity_thres)

sv_match_salvage_paf_2_list <- 
  map(c("sv_match_id", "sv_partial_match_id"), ~ {
    candi_df_mod %>% 
      filter(seq_id %in% sv_match_salvage[[.x]]) %>% 
      combine_seqid_idx()
  }) %>% 
  `names<-`(c("full", "partial"))
  
sv_match_salvage_variant_list <- 
  map(sv_match_salvage_paf_2_list, convert_paf_variant_multisample) %>% 
  `names<-`(c("full", "partial"))

sv_match_salvage_gtf_list <- 
  map2(
    sv_match_salvage_paf_2_list,
    sv_match_salvage_variant_list,
    ~ retrieve_gtf_multisample(.x, "exon", .y)
  ) %>% 
  `names<-`(c("full", "partial"))

sv_match_salvage_intron_list <- 
  map2(
    sv_match_salvage_paf_2_list,
    sv_match_salvage_variant_list,
    ~ retrieve_gtf_multisample(.x, "intron", .y)
  ) %>% 
  `names<-`(c("full", "partial"))

improper_sj_id_match_salvage_list <- 
  map(
    sv_match_salvage_intron_list,
    get_id_noncanonical_notinref_notinsqanti,
    ref_intron,
    sqanti_intron_additional
  ) %>% 
  `names<-`(c("full", "partial"))

chimera_salvage_match_gtf_list <- 
  map2(
    sv_match_salvage_gtf_list,
    improper_sj_id_match_salvage_list,
    ~ {
      .x %>% 
        separate_seqid_idx() %>% 
        filter(!seq_id %in% .y) %>% # update
        add_exon_number()
    }
  ) %>% 
  `names<-`(c("full", "partial"))

chimera_salvage_match_variant_list <- 
  map2(
    sv_match_salvage_variant_list,
    improper_sj_id_match_salvage_list,
    ~ {
      .x %>% 
        separate_seqid_idx() %>% 
        filter(!seq_id %in% .y) # update
    }
  ) %>% 
  `names<-`(c("full", "partial"))

chimera_salvage_match_sv_list <- 
  map2(
    c("sv_match", "sv_partial_match"),
    improper_sj_id_match_salvage_list,
    ~ sv_match_salvage[[.x]] %>% filter(!seq_id %in% .y) %>% replace_df_mod( , "sv_type", "sj_not_confirmed")
    ) %>% 
  `names<-`(c("full", "partial"))

sv_multi_match_salvage_id <-
  chimera_salvage_match_sv_list$full %>% 
  distinct(seq_id, sample, up_chr, up_strand, up_pos, down_chr, down_strand, down_pos) %>% 
  count(seq_id) %>% 
  filter(n > 1L) %>% 
  `[[`("seq_id")
sv_uni_match_salvage_id <- 
  chimera_salvage_match_sv_list$full$seq_id %>% 
  unique() %>% 
  setdiff(sv_multi_match_salvage_id)
sv_partial_match_salvage_id <-
  chimera_salvage_match_sv_list$partial$seq_id %>% 
  unique()

# declare result object ----
result_unit <- vector(mode = "list", length = 4L) %>% `names<-`(c("id", "gtf", "sv", "variant"))

final_result <- 
  list(
    sv_uni_match = result_unit,
    sv_multi_match = result_unit,
    sv_partial_match = result_unit
  )

# uni_match
final_result$sv_uni_match$id <- c(sv_uni_match_id, sv_uni_match_salvage_id)
final_result$sv_uni_match$sv <- 
  bind_rows(
    chimera_res_sv_match_list$full$sv_match %>% filter(seq_id %in% sv_uni_match_id),
    chimera_salvage_match_sv_list$full %>% filter(seq_id %in% sv_uni_match_salvage_id)
  )
final_result$sv_uni_match$gtf <- 
  bind_rows(
    chimera_res_gtf_list$full %>% filter(seq_id %in% sv_uni_match_id),
    chimera_salvage_match_gtf_list$full %>% filter(seq_id %in% sv_uni_match_salvage_id)
  )
final_result$sv_uni_match$variant <- 
  bind_rows(
    chimera_res_variant_list$full %>% filter(seq_id %in% sv_uni_match_id),
    chimera_salvage_match_variant_list$full %>% filter(seq_id %in% sv_uni_match_salvage_id)
  )

# multi_match
final_result$sv_multi_match$id <- 
  c(sv_multi_match_id, sv_multi_match_add_id, sv_multi_match_salvage_id)
final_result$sv_multi_match$sv <- 
  bind_rows(
    chimera_res_sv_match_list$full$sv_match %>% filter(seq_id %in% sv_multi_match_id),
    chimera_res_sv_match_list$amb_partial$sv_match %>% filter(seq_id %in% sv_multi_match_add_id),
    chimera_res_sv_match_list$int_partial$sv_match %>% filter(seq_id %in% sv_multi_match_add_id),
    chimera_salvage_match_sv_list$full %>% filter(seq_id %in% sv_multi_match_salvage_id)
  )
final_result$sv_multi_match$gtf <-
  bind_rows(
    chimera_res_gtf_list$full %>% filter(seq_id %in% sv_multi_match_id),
    chimera_res_gtf_list$partial %>% filter(seq_id %in% sv_multi_match_add_id),
    chimera_salvage_match_gtf_list$full %>% filter(seq_id %in% sv_multi_match_salvage_id)
  )
final_result$sv_multi_match$variant <- 
  bind_rows(
    chimera_res_variant_list$full %>% filter(seq_id %in% sv_multi_match_id),
    chimera_res_variant_list$partial %>% filter(seq_id %in% sv_multi_match_add_id),
    chimera_salvage_match_variant_list$full %>% filter(seq_id %in% sv_multi_match_salvage_id)
  )

# partial_match
final_result$sv_partial_match$id <- 
  c(sv_partial_match_id, sv_partial_match_add_id, sv_partial_match_salvage_id)

final_result$sv_partial_match$sv <- 
  bind_rows(
    chimera_res_sv_match_list$full$sv_partial_match, 
    bind_rows(
      chimera_res_sv_match_list$amb_partial$sv_match,
      chimera_res_sv_match_list$amb_partial$sv_partial_match,
      chimera_res_sv_match_list$amb_partial$sv_unmatch,
      chimera_res_sv_match_list$int_partial$sv_match,
      chimera_res_sv_match_list$int_partial$sv_partial_match,
      chimera_res_sv_match_list$int_partial$sv_unmatch
    ) %>% 
      filter(seq_id %in% sv_partial_match_add_id),
    chimera_salvage_match_sv_list$partial
  ) %>% 
  distinct()
final_result$sv_partial_match$gtf <-
  bind_rows(
    chimera_res_gtf_list$full %>% filter(seq_id %in% sv_partial_match_id),
    chimera_res_gtf_list$partial %>% filter(seq_id %in% sv_partial_match_add_id),
    chimera_salvage_match_gtf_list$partial
  )
final_result$sv_partial_match$variant <-
  bind_rows(
    chimera_res_variant_list$full %>% filter(seq_id %in% sv_partial_match_id),
    chimera_res_variant_list$partial %>% filter(seq_id %in% sv_partial_match_add_id),
    chimera_salvage_match_variant_list$partial
  )


# popular style ----
change_pm <- function(x) factor(x, levels = c("+", "-")) %>% `levels<-`(c("-", "+")) %>% as.character()
final_result$sv_uni_match$sv <- final_result$sv_uni_match$sv %>% mutate(down_strand = change_pm(down_strand))
final_result$sv_multi_match$sv <- final_result$sv_multi_match$sv %>% mutate(down_strand = change_pm(down_strand))
final_result$sv_partial_match$sv <- final_result$sv_partial_match$sv %>% mutate(down_strand = change_pm(down_strand))

# output ----
saveRDS(final_result, chimera_output)

combine_seqid_idx <- function(df)  mutate(df, seq_id = str_c(seq_id, idx, sep = "___"), idx = NULL)
bind_rows(
  final_result$sv_uni_match$gtf,
  final_result$sv_multi_match$gtf,
  final_result$sv_partial_match$gtf
) %>%
    combine_seqid_idx() %>%
    select(seqname = "chr", strand, start, end, transcript_id = "seq_id", feature = "type") %>%
    mutate(gene_id = "gene_placeholder") %>%
    write_gtf(gene_fragment_gtf)
