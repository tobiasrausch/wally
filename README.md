[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/wally/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/wally/badges/downloads.svg)](https://anaconda.org/bioconda/wally)
[![C/C++ CI](https://github.com/tobiasrausch/wally/workflows/C/C++%20CI/badge.svg)](https://github.com/tobiasrausch/wally/actions)
[![Docker CI](https://github.com/tobiasrausch/wally/workflows/Docker%20CI/badge.svg)](https://hub.docker.com/r/trausch/wally/)
[![GitHub license](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/tobiasrausch/wally/blob/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/tobiasrausch/wally.svg)](https://github.com/tobiasrausch/wally/releases)

# Wally

Plotting of aligned sequencing reads, assembled contigs or pan-genome graphs in BAM/CRAM format and visualization of genomic variants.

## Installing Wally

Wally is available as a [Bioconda package](https://anaconda.org/bioconda/wally), as a pre-compiled static binary from the [release page](https://github.com/tobiasrausch/wally/releases/), as a singularity container [SIF file](https://github.com/tobiasrausch/wally/releases/) or as a minimal [Docker container](https://hub.docker.com/r/trausch/wally/). You can also build Wally from source using a recursive clone and make. 

`git clone --recursive https://github.com/tobiasrausch/wally.git`

`cd wally/`

`make all`


## Running Wally

Wally uses subcommands for different visualization modes.

## Subcommand `region`: Visualization of alignment regions

Wally needs a sorted and indexed BAM/CRAM file and a reference genome. The output is a genomic alignment plot of the region specified on the command-line.

`wally region -r chrA:35-80 -g <genome> <input.bam>`

You can specify a plot name using another colon separator

`wally region -r chrA:35-80:myplot -g <genome> <input.bam>`

Most often you probably want to use a BED file with regions of interest and just execute wally in batch. The 4-th column specifies the plot name. For instance to plot variants from a VCF file you can use:

`bcftools query -f "%CHROM\t%POS\n" <input.vcf.gz> | awk '{print $1"\t"($2-50)"\t"($2+50)"\tid"NR;}' > regions.bed`

`wally region -R regions.bed -g <genome> <input.bam>`

### Gene annotations

Simple BED files are used to provide gene or other annotations. The BED file needs to be bgzipped and indexed via tabix (`tabix -p bed input.bed.gz`). The required columns are chromosome, start, end and an identifier which is displayed if there is sufficient space.

`wally region -b bed/gencode.hg19.bed.gz -r chr17:7573900-7574000 -g <genome> <input.bam>`

The chromosome names of the genome FASTA file, BAM file and BED annotation file need to match, "chr11" and "11" are not considered identical.

### Custom colors for BED file annotations

You can specify custom colors for each annotation record in the BED file using Hex color codes (column 5 of the BED file). For instance, an annotation file `anno.bed` with this content:

```
chr17	7350100	7350900	regRed	0xFF0000
chr17	7351100	7351900	regLime	0x00FF00
chr17	7352100	7352900	regBlue	0x0000FF
chr17	7353100	7353900	regTeal	0x008080
```

can be visualized using

`bgzip anno.bed`

`tabix -p bed anno.bed.gz`

`wally region -b anno.bed.gz -r chr17:7350000-7354000 -g <genome> <input.bam>`


### Multiple BAM files

You can include multiple BAM files in a plot such as a tumor genome and a matched control in cancer genomics.

`wally region -r chr17:7573900-7574000 -g <genome> <tumor.bam> <control.bam>`

In fact, wally can be used to create "wallpapers" of genomes which gave rise to the tool name. For instance, the below command can be used to create a full genome view of 96 samples of a full plate of SARS-CoV-2 genomes.

`wally region -y 20480 -r NC_045512.2:1-30000 -g NC_045512.2.fa Plate*.bam`

### Aligning multiple regions horizontally (split view)

Wally allows concatenating images of different regions horizontally using the `--split` option. This can be used, for instance, to zoom into a specific variant. On the command line the regions need to be separated by `,` without spaces. As an example, you can zoom into the N501Y variant of the alpha SARS-CoV-2 lineage using 

`wally region -s 3 -x 2048 -r NC_045512.2:22000-24000,NC_045512.2:23000-23100,NC_045512.2:23050-23070 -g NC_045512.2.fa <input.bam>`

You can split horizontally and vertically at the same time to view, for instance, a somatic inter-chromosomal translocation.

`wally region -s 2 -r chrA:35-80,chrB:60-80 -g <genome> <tumor.bam> <control.bam>`

If you specify the regions in a BED file using the `-R` option then the split parameter operates row-wise, e.g., for `-s 3` row 1-3 of the BED file make up the first image, row 4-6 the second image, and so on.

### Paired-end view

With `-p` you can switch on the paired-end view. 

`wally region -p -r chrA:35-80 -g <genome> <input.bam>`

The paired-end coloring highlights candidate structural variants supported by read1 (R1) and read2 (R2). Below is a mapping of [delly's](https://github.com/dellytools/delly) structural variant types to [wally's](https://github.com/tobiasrausch/wally) paired-end coloring. For inter-chromosomal translocations I assumed for this illustration that R1 maps to chromosome A and R2 maps to chromosome B.

- ![#3fafaf](https://via.placeholder.com/15/3fafaf/3fafaf.png) `  inversion-type paired-end, --R1-->   --R2-->, INV:3to3`
- ![#4e64d5](https://via.placeholder.com/15/4e64d5/4e64d5.png) `  inversion-type paired-end, <--R1--   <--R2--, INV:5to5`
- ![#d53f3f](https://via.placeholder.com/15/d53f3f/d53f3f.png) `   deletion-type paired-end, --R1-->   <--R2--, DEL:3to5`
- ![#3faf3f](https://via.placeholder.com/15/3faf3f/3faf3f.png) `duplication-type paired-end, <--R1--   --R2-->, DUP:5to3`
- ![#fb9a99](https://via.placeholder.com/15/fb9a99/fb9a99.png) `inter-chr paired-end, A:--R1--> B:--R2--> leads to --A--> <--B-- junction, BND:3to3`
- ![#fdbf6f](https://via.placeholder.com/15/fdbf6f/fdbf6f.png) `inter-chr paired-end, A:<--R1-- B:<--R2-- leads to <--A-- --B--> junction, BND:5to5`
- ![#ff8000](https://via.placeholder.com/15/ff8000/ff8000.png) `inter-chr paired-end, A:--R1--> B:<--R2-- leads to --A--> --B--> junction, BND:3to5`
- ![#cab2d6](https://via.placeholder.com/15/cab2d6/cab2d6.png) `inter-chr paired-end, A:<--R1-- B:--R2--> leads to <--A-- <--B-- junction, BND:5to3`

For large and complex structural variants, wally supports split views (as explained above). For instance, for an inter-chromosomal translocation you probably want to use a 2-way horizontal split with a larger image size.

`wally region -p -x 2048 -y 2048 -s 2 -r chrA:35-80,chrB:60-80 -g <genome> <tumor.bam> <control.bam>`

To visualize genomic breakpoints it's also helpful to highlight clipped reads `-c` and supplementary alignments `-u`.

`wally region -cup -x 2048 -y 2048 -s 2 -r chrA:35-80,chrB:60-80 -g <genome> <tumor.bam> <control.bam>`

## Subcommand `matches`: Visualization of chained alignment matches

For long read alignments or mappings of entire contigs, one often wants to explore all matches of a long input sequence with respect to the reference to visualize, for instance, complex rearrangements. For a single read or contig, you need a sorted and indexed BAM file and the reference genome.

`wally matches -r <read_name> -g <genome> <input.bam>`

Forward matches are colored in blue, reverse matches in orange and the black line traces all alignment matches along the input read sequence. You can also plot multiple reads to inspect, for instance, alignment concordance or heterogeneity using a list of reads as an input file.

`wally matches -R <read.lst> -g <genome> <input.bam>`

For many reads, this quickly leads to "wall(y)papers". With `-s` you get a separate plot for each read.

`wally matches -s -R <read.lst> -g <genome> <input.bam>`

You can also compare assemblies or align contigs to a reference genome using a workflow like this:

`minimap2 -ax asm5 ref.fa asm.fa | samtools sort -o assembly.bam -`

`samtools index assembly.bam`

`wally matches -r <contig_name> -g ref.fa assembly.bam`

## Subcommand `dotplot`: Pairwise sequence dotplots

All pairwise dotplots of a FASTA file can be generated using

`wally dotplot sequences.fa`

You can also extract reads directly from the BAM file and generate a reference-independent dotplot that includes the genomic mapping locations of read segments.

`wally dotplot -R reads -g hg38.fa input.bam`

For instance, to inspect haplotype differences or tumor heterogeneity at a given genomic locus (e.g., `chrA:1000-2000`) one can use:

`samtools view input.bam chrA:1000-2000 | cut -f 1 | sort | uniq > reads`

`wally dotplot -R reads -g hg38.fa input.bam`

Instead of plotting reads in pairs against each other, you can also calculate dot plots with respect to reference regions.

`wally dotplot -g hg38.fa -R reads -e chrA:35000-80000 input.bam`

`wally dotplot -g hg38.fa -e chrA:35000-80000 sequences.fa`

A minimum contig or read length can be specified using `-s`

`wally dotplot -s 10000 -R reads -g hg38.fa input.bam`

Lastly, you can flatten the genomic mappings into simple blocks using `-f`

`wally dotplot -f -s 10000 -R reads -g hg38.fa input.bam`

To visualize a self-alignment of a read with reference mappings use `-a`

`wally dotplot -a -r read_name -g hg38.fa input.bam`

## Web application

You can try wally on some selected data sets using the [wally web application](https://www.gear-genomics.com/wally/).

## Citation

Tobias Rausch, Rene Snajder, Adrien Leger, Milena Simovic, Mădălina Giurgiu, Laura Villacorta, Anton G. Henssen, Stefan Fröhling, Oliver Stegle, Ewan Birney, Marc Jan Bonder, Aurelie Ernst, Jan O. Korbel     
Long-read sequencing of diagnosis and post-therapy medulloblastoma reveals complex rearrangement patterns and epigenetic signatures     
Cell Genomics, 2023, 100281, [DOI: 10.1016/j.xgen.2023.100281](https://doi.org/10.1016/j.xgen.2023.100281)     

## License

Wally is distributed under the BSD 3-Clause license. Consult the accompanying [LICENSE](https://github.com/tobiasrausch/wally/blob/master/LICENSE) file for more details.

## Credits

Wally relies heavily on the [HTSlib](https://github.com/samtools/htslib) and [OpenCV](https://github.com/opencv/opencv). The visualization of genomic alignments was heavily inspired by [IGV](https://github.com/igvteam/igv).
