[![C/C++ CI](https://github.com/tobiasrausch/wally/workflows/C/C++%20CI/badge.svg)](https://github.com/tobiasrausch/wally/actions)
[![Docker CI](https://github.com/tobiasrausch/wally/workflows/Docker%20CI/badge.svg)](https://hub.docker.com/r/trausch/wally/)
[![GitHub license](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/tobiasrausch/wally/blob/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/tobiasrausch/wally.svg)](https://github.com/tobiasrausch/wally/releases)

# Wally

Plotting of aligned sequencing reads in BAM/CRAM format and visualization of genomic variants.

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

# Aligning multiple regions horizontally (split view)

Wally allows concatenating images of different regions horizontally using the `--split` option. This can be used, for instance, to zoom into a specific variant. On the command line the regions need to be separated by `,` without spaces. As an example, you can zoom into the N501Y variant of the alpha SARS-CoV-2 lineage using 

`wally region -s 3 -x 2048 -r NC_045512.2:22000-24000,NC_045512.2:23000-23100,NC_045512.2:23050-23070 -g NC_045512.2.fa <input.bam>`

You can split horizontally and vertically at the same time to view, for instance, a somatic inter-chromosomal translocation.

`wally region -r chrA:35-80,chrB:60-80 -g <genome> <tumor.bam> <control.bam>`

If you specify the regions in a BED file using the `-R` option then the split parameter operates row-wise, e.g., for `-s 3` row 1-3 of the BED file make up the first image, row 4-6 the second image, and so on.

# Paired-end view

With `-p` you can switch on the paired-end view. 

`wally region -p -r chrA:35-80 -g <genome> <input.bam>`

The paired-end coloring highlights candidate structural variants supported by read1 (R1) and read2 (R2). Below is a mapping of [delly's](https://github.com/dellytools/delly) structural variant types to [wally's](https://github.com/tobiasrausch/wally) paired-end coloring. For inter-chromosomal translocations I assumed for this illustration that R1 maps to chromosome A and R2 maps to chromosome B.

- ![#3fafaf](https://via.placeholder.com/15/3fafaf/000000?text=+) `  inversion-type paired-end, --R1-->   --R2-->, INV:3to3`
- ![#4e64d5](https://via.placeholder.com/15/4e64d5/000000?text=+) `  inversion-type paired-end, <--R1--   <--R2--, INV:5to5`
- ![#d53f3f](https://via.placeholder.com/15/d53f3f/000000?text=+) `   deletion-type paired-end, --R1-->   <--R2--, DEL:3to5`
- ![#3faf3f](https://via.placeholder.com/15/3faf3f/000000?text=+) `duplication-type paired-end, <--R1--   --R2-->, DUP:5to3`
- ![#fb9a99](https://via.placeholder.com/15/fb9a99/000000?text=+) `inter-chr paired-end, A:--R1--> B:--R2--> leads to --A--> <--B-- junction, BND:3to3`


For large and complex structural variants, wally supports split views (as explained above). For instance, for an inter-chromosomal translocation you probably want to use a 2-way horizontal split.

`wally region -up -x 2048 -s 2 -r chrA:35-80,chrB:60-80 -g <genome> <tumor.bam> <control.bam>`


# License

Wally is distributed under the BSD 3-Clause license. Consult the accompanying [LICENSE](https://github.com/tobiasrausch/wally/blob/master/LICENSE) file for more details.
