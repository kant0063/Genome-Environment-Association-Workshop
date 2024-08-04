# Getting Genomic Data Ready for Use for the GEA analysis

There are various data types that are commonly used for genomic analysis. The final format of the data will usually be .vcf or .hmp or .fam/.bim/.bed or .txt, these data require different levels of quality assurance and quality control. These data are also usually associated with different genotyping methods and different types of collections. 

# Common types of Genotyping
(1) Whole Genome shotgun sequencing
(2) Genotype-by-Seqeuncing or RAD-Sequencing
(3) DNA microarray (SNP-chip)
(4) Amplicon sequencing
(5) Sequence Capture array
(6) Copmarative Gneomic Hybridization (CGH) Array

# Data Cleaning pipeline
(1) Get .fastq file off of the sequncing machine
(2) Explore quality of .fastq file using FastQC
(3) Trim .fastq file
  - remove adapater sequence
  - remove poor quality bases
  - Common software include Sythe and trimmomatic
(4) Get reference genome
  - Reference genomes are housed at a number of different places
  - phytozome - https://phytozome-next.jgi.doe.gov/
  - species specific databases such as MazieGDB - https://www.maizegdb.org/ or Soybase - https://www.soybase.org/
  - The National Center for Biotechnology Information - https://www.ncbi.nlm.nih.gov/
  - The European Nucleotide Archive - https://www.ebi.ac.uk/ena/browser/
(5) Align to reference genome
  - Decide on an alignment software
  - Common software include BWA-MEM and Bowtie
(6) Generate SNP genotypes
  - Common tools include vcftools, GATK, freeBayes
(7) Output a .vcf file
