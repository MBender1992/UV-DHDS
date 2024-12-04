# https://cran.r-project.org/web/packages/vcfR/vignettes/intro_to_vcfR.html
library(vcfR)
pkg <- "pinfsc50"
vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = pkg)
dna_file <- system.file("extdata", "pinf_sc50.fasta", package = pkg)
gff_file <- system.file("extdata", "pinf_sc50.gff", package = pkg)

vcf <- read.vcfR( vcf_file, verbose = FALSE )
dna <- ape::read.dna(dna_file, format = "fasta")
gff <- read.table(gff_file, sep="\t", quote="")


chrom <- create.chromR(name='Supercontig', vcf=vcf, seq=dna, ann=gff)

plot(chrom)
chrom <- masker(chrom, min_QUAL = 1, min_DP = 300, max_DP = 700, min_MQ = 59.9,  max_MQ = 60.1)
plot(chrom)

chrom <- proc.chromR(chrom, verbose=TRUE)

plot(chrom)

chromoqc(chrom, dp.alpha=20)

## now the filtered and quality checked data can be stored as a new vcf file
