You can build a [wally](https://github.com/tobiasrausch/wally) singularity container (SIF file) using

`sudo singularity build wally.sif wally.def`

Once you have built the container you can run analysis using

`singularity exec wally.sif wally --help`
