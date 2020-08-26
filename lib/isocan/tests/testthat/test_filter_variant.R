# context("filter_variant")
# test_that("direct position",
#           {
#             genomon_1 <- testthis::read_testdata("genomon_1.rds")
#             variant_df <- testthis::read_testdata("variant_df.rds")
#             res_1 <- testthis::read_testdata("res_1.rds")
#             hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#             res <- filter_variant(variant_df, genomon_1, hg38)
#             expect_identical(res, res_1)
#           })
# test_that("alternative position",
#           {
#             genomon_2 <- testthis::read_testdata("genomon_2.rds")
#             variant_df <- testthis::read_testdata("variant_df.rds")
#             res_2 <- testthis::read_testdata("res_2.rds")
#             hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#             res <- filter_variant(variant_df, genomon_2, hg38)
#             expect_identical(res, res_2)
#           })
