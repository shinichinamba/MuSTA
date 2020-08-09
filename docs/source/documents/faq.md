# FAQ

## Resuming MuSTA

MuSTA utilizes {drake} R package, which recodes the progress of the pipeline into `.drake` folder and enables resuming at a suspended point.

When MuSTA is suspended for any reasons, users can simply re-run it.

If you'd like to re-run the whole procedure, you can use the '--force' option (or delete the `.drake` directory located at your working directory).


## Using MuSTA in an HPC environment

Currently, MuSTA does not support HPC architecture, except for SHIROKANE super computer in Human Genome Center, Japan.

Users working with other HPC systems need to run MuSTA without batch submission until future updates.


## The status appearred in 'report/plan_post_run.pdf' is still 'outdated' even after the run has been completed.

If required time is written under each target, these targets seem to be successfully completed, but somehow {drake} doesn't know it.

You can check it by running MuSTA with the same option as your last run, and an additional '--dry-run' option.

It generates 'plan_pre_run.pdf' in your current directory, and it should be identical with 'report/plan_post_run.pdf' of your last run.