[![C/C++ CI](https://github.com/tobiasrausch/wally/workflows/C/C++%20CI/badge.svg)](https://github.com/tobiasrausch/wally/actions)
[![Docker CI](https://github.com/tobiasrausch/wally/workflows/Docker%20CI/badge.svg)](https://hub.docker.com/r/trausch/wally/)
[![GitHub license](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/tobiasrausch/wally/blob/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/tobiasrausch/wally.svg)](https://github.com/tobiasrausch/wally/releases)

# Wally

A genomic variant plotter.

# Installing Wally

The easiest way to get Wally is to download the statically linked binary from the [Wally github release page](https://github.com/tobiasrausch/wally/releases/). You can also build Wally from source using a recursive clone and make. 

`git clone --recursive https://github.com/tobiasrausch/wally.git`

`cd wally/`

`make all`


# Running Wally

Wally needs a sorted and indexed BAM/CRAM file and a reference genome. The output is a genomic alignment plot of the region specified on the command-line.

`wally region -r chrA:35-80 -g <genome> <input.bam>`

You can specify a plot name using another colon separator

`wally region -r chrA:35-80:myplot -g <genome> <input.bam>`

Most often you probably want to use a BED file with regions of interest and just execute wally in batch. The 4-th column specifies the plot name. For instance to plot variants from a VCF file you can use:

`bcftools query -f "%CHROM\t%POS\n" <input.vcf.gz> | awk '{print $1"\t"($2-50)"\t"($2+50)"\tid"NR;}' > regions.bed`

`wally region -R regions.bed -g <genome> <input.bam>`

# Gene annotations

Simple BED files are used to provide gene or other annotations. The BED file needs to be bgzipped and indexed via tabix (`tabix -p bed input.bed.gz`). The required columns are chromosome, start, end and an identifier which is displayed if there is sufficient space.

`wally region -b bed/gencode.hg19.bed.gz -r chr17:7573900-7574000 -g <genome> <input.bam>`

The chromosome names of the genome FASTA file, BAM file and BED annotation file need to match, "chr11" and "11" are not considered identical.

# Multiple BAM files

You can include multiple BAM files in a plot such as a tumor genome and a matched control in cancer genomics.

`wally region -r chr17:7573900-7574000 -g <genome> <tumor.bam> <control.bam>`

In fact, wally can be used to create "wallpapers" of genomes which gave rise to the tool name. For instance, the below command can be used to create a full genome view of 96 samples of a full plate of SARS-CoV-2 genomes.

`wally region -y 20480 -r NC_045512.2:1-30000 -g NC_045512.2.fa Plate*.bam`

# Missing features (work-in-progress)

* Paired-end support

* Abnormal paired-end coloring

* Split view for SVs


# License

Wally is distributed under the BSD 3-Clause license. Consult the accompanying [LICENSE](https://github.com/tobiasrausch/wally/blob/master/LICENSE) file for more details.
