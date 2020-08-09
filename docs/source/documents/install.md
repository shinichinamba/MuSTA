# Installation

Our GitHub repogitory can be found [here](https://github.com/shinichinamba/MuSTA).

You can obtain MuSTA by typing

```
git clone https://github.com/shinichinamba/MuSTA.git
```


## R dependencies

MuSTA depends mainly on R and indirectly on python (via {argparse} R package), and should be run on unix OS because MuSTA utilizes some shell scripts.

As MuSTA always checks whether its dependent R packages are installed, you can use this function by running the code:

```
cd MuSTA
./MuSTA.R -v
```

If you see some messages showing that there are missing packages, start R and copy-and-paste the messages to the R console.

If you successfully see the version of MuSTA, now you have completed the installation!


## External softwares

The codes and R packages located at 'lib' directory can be run without any other requirement.

However, the MuSTA pipeline utilizes several external softwares, and you need to install them in order to use the main pipeline.

|  Softwares  |  URLs  |
| ----------- | ------ |
|  samtools   |  [http://www.htslib.org/](http://www.htslib.org/) |
|  Salmon     |  [https://salmon.readthedocs.io](https://salmon.readthedocs.io) |
|  minimap2   |  [https://github.com/lh3/minimap2](https://github.com/lh3/minimap2)  |
|  SQANTI     |  [https://github.com/ConesaLab/SQANTI](https://github.com/ConesaLab/SQANTI)  |
|  LoRDEC (optional)  |  [https://gite.lirmm.fr/lordec/lordec-releases/-/wikis/home](https://gite.lirmm.fr/lordec/lordec-releases/-/wikis/home)  |
|  seqkit (optional)  |  [https://bioinf.shenwei.me/seqkit/](https://bioinf.shenwei.me/seqkit/)  |



## Inform MuSTA of the paths to external softwares

MuSTA uses several external softwares listed above.

The locations of these softwares need to be informed to MuSTA either by exporting paths or by specifying appropriate optional arguments.


### Exporting paths

If you export all the paths so that MuSTA can simply use them with their basenames, you need not to write these paths to your config file.

MuSTA uses `bash` for calling external softwares, so the export should be done in the `bash` environment.


