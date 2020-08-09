# Editing a config file

MuSTA requires a config file which contains file paths and other optional configulations such as paths to external softwares.

The config file is written in the YAML format and includes at least 3 mandatory fields, `default`, `general`, and `long_read_hq`.

An example can be found at https://github.com/sinnhazime/MuSTA/test/config.yml.
You can easily write your own configulation by editing this file.


## YAML format

YAML is a sophisticated readable file format for representing structured data.

There are many websites describing YAML like (here)[https://rollout.io/blog/yaml-tutorial-everything-you-need-get-started/], so
we explain a brief overview with our example config file.


```
general:
  samples: [sample1, sample2, sample3, sample4]
  # Sample identifiers. Must not be overlapped each other.
  output_dir: test/test_musta
  reference_gtf: test/reference/test.ref.gtf
  genome_fasta: test/reference/chr16_18_head.fa

```

YAML utilizes (white-space) indentication for representing items included in other items.

In this example, **general** has 4 child items, **samples**, **output_dir**, **reference_gtf**, and **genome_fasta**.

Items can contain multiple values like the **samples** field.

Characters after "#" are treated as comment-out.

In addition to the basic YAML format, you can use some shortcuts (**dir**, **ext**, and **inherits**).
The details are described in the [Using shortcuts](#Using shortcuts) section.


## Mandatory fields

The minimum config file can be comprised of `default`, `general`, and `long_read_hq` fields.
Since `default` field needs not to be edited (as described below), only two 


### default

**default** is a special field for MuSTA, as its contents are inherited in any ather fields (i.e. the contents in **default** are treated as default values for other fields.) 
This feature is due to [{config} R package](https://github.com/rstudio/config), which MuSTA depends on.

We recommend to keep this field untouched in order to avoid unnecessary confusion.


### long_read_hq



## Optional fields

### short_read

### long_read_lq

### cluster_report


## Using shortcuts

### dir/ext

### inherits
