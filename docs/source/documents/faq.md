# FAQ

## Resuming MuSTA

MuSTA utilizes {drake} R package, which recodes the progress of the pipeline into `.drake` folder and enables resuming at a suspended point.

When MuSTA is suspended for any reasons, users can simply re-run it.

If you'd like to re-run the whole procedure, you can use the '--force' option (or delete the `.drake` directory located at your working directory).


## Using MuSTA in an HPC environment

Currently, MuSTA does not support HPC architecture, except for SHIROKANE super computer in Human Genome Center, Japan.

Users working with other HPC systems need to run MuSTA without batch submission until future updates.

