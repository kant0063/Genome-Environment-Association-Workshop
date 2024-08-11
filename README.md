# Genome Enviroment Association Workshop

This is a series of lessons accompanied by scripts to conduct Genome Environmental Association Analysis and Genomic Environmental Selection Analysis.

The observation unit for the data in this workshop is an `accession` with the location where it was collected as a `latitude` and `longitude`. The problem is that this information is disconnected from the data capture processes for the environmental data and the lab-based sequencing workflows for genotyping. An highly generalized workflow for GEA can be seen below.

```mermaid
graph LR
  Env1[Globally Gridded Climate Data -- GeoTIFFs or other raster format]
  Env2[Globally Gridded Soils Data -- GeoTIFFs or other raster format]
  Env3[Spatially Joined Globally Gridded Climate and Soils Data]
  Env1 --> Env3
  Env2 --> Env3

  Gen1[Accession Sequencing Data .fastq files]
  Gen2[Reference Genome]
  Gen3[Cleaned, Aligned, and SNP Called Data .vcf file]
  Gen4[Transformed Numeric SNPs .txt file]
  Gen1 --> Gen2 --> Gen3 --> Gen4

  Accession1[Field Sample Collection]
  Accession2[Accession Unique ID with Location csv or xlsx]
  Accession1 --> Accession2

  Merge1((Join using Accession ID))
  IntegratedData1[Integrated Genetic and Accession CSV]
  Gen4 --> Merge1
  Accession2 --> Merge1
  Merge1 --> IntegratedData1


  Merge2((Spatial Join by Accession Latitude and Longitude))
  IntegratedData2[Spatially Join Soil, \n Climate, Genotyping, and Accession Data]
  Env3 --> Merge2
  IntegratedData1 --> Merge2 --> IntegratedData2


```
# Steps
1) Getting accession geolocation information and mapping accession location
2) Accessing genomic data
3) Accessing climate and soil data 
4) Cleaning genomic data
5) Analysing population structure
6) Calculating LD decay
7) Performing the EAA (Environmental Association Analysis) & EAA (Environmental Association Analysis) Model Checking
8) Identifying candidate genomic regions and genes in LD with those regions
9) Performing the EGS (Environmental Genomic Selection)

