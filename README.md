# Wally

A genomic variant plotter.

# Installing Wally

The easiest way to get Wally is to download a statically linked binary or the singularity container (SIF file) from the [Wally github release page](https://github.com/tobiasrausch/wally/releases/). You can also build Wally from source using a recursive clone and make. 

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

# Missing features (work-in-progess)

* Paired-end support

* Abnormal paired-end coloring

* Split view for SVs

* GTF support

* Multiple BAM files as input

# License

Wally is distributed under the BSD 3-Clause license. Consult the accompanying [LICENSE](https://github.com/tobiasrausch/wally/blob/master/LICENSE) file for more details.
