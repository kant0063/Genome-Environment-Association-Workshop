# Environmental Association Analysis
## Background
Many approaches have been developed to explore the importance of genetic local adaptation. Many of these have been based on allele frequency (e.g. Wright 1949; Lewontin and Krakauer 1973; Beaumont and Balding 2004; Beaumont 2005). Other methods have focused on exploring changes in allele frequency gradients (Yang et al., 2012), correlating allele frequencies with environmental variables (e.g. Gunther and Coop, 2013) or use mixed model approaches to identify an association between genetic variants and environmental variables (Eckert et al. 2010; Pyhäjärvi et al. 2013; Yoder et al. 2014).
## Mixed-Model Approach
Mixed-model association can be implemented in many different software packages. The general form of this model is 
Y = Qw + βx + Zu + e⁠, where y is a vector of bioclimatic or soil data, Q the incidence matrix of PCs (if modeling population structure), w the vector of fixed PC effects, X the SNP marker genotype matrix, β the fixed SNP effect being tested, Z the incidence matrix of individuals, u ~MVN(0,K$\theta^2)

SNPs and bioclimatic and biophysical variables. The following modelswereexplored:anaïvemodelwithnocontrolforpopulation structure, a model using the Q-matrix from STRUCTURE to control for population structure, a model using a kinship matrix (K-matrix) based on pairwise Manhattan distance amongindividuals at all SNPs, a model using both a K-matrix and a Q-matrix,and models that integrate latitude or latitude and longitude as covariates. 
## Model Testing
Quantile-Quantile (qq) plots were examined for each model, and the genomic inflation parameter lambda(L)(DevlinandRoeder1999),ametricofp-valueinflation, was calculated (Figure S3). Athreshold set at the 0.01% most extreme p-values was used to identify candidates for each environmental association, resulting in three significant marker associations per environmental variable. This threshold was chosen to focus the analysis on a minimum number of large effect variants, and to limit false positives. The qqman R package (Turner 2014) was used to plot the association results.

- Beaumont, M. A., 2005 Adaptation and speciation: what can Fst tell us? Trends Ecol. Evol. 20: 435–440. 
- Beaumont, M. A., D. J. Balding, 2004 Identifying adaptive genetic divergence among populations from genome scans. Mol. Ecol. 13: 969–980.
- Eckert, A. J., J. van Heerwaarden, J. L. Wegrzyn, C. D. Nelson, J. Ross-Ibarra et al., 2010 Patterns of population structure and environmental associations to aridity across the range of loblolly pine (Pinus taeda L., Pinaceae). Genetics 185: 969–982.
- Günther, T., Coop, G. 2013. Robust identification of local adaptation from allele frequencies. Genetics, 195(1), 205-220.
- Lewontin, R. C., and J. Krakauer, 1973 Distribution of gene frequency as a test of the theory of the selective neutrality of polymorphisms. Genetics 74: 175–195.
- Pyhäjärvi, T., M. B. Hufford, S. Mezmouk, and J. Ross-Ibarra, 2013 Complex patterns of local adaptation in teosinte. Genome Biol. Evol. 5: 1594–1609.
- Wright, S. 1949. The genetical structure of populations. Ann. Eugen. 15: 323–354.
- Yang, W. Y., J. Novembre, E. Eskin, and E. Halperin, 2012 A model-based approach for analysis of spatial structure in genetic data. Nat. Genet. 44: 725–731. 
- Yoder, J. B., J. Stanton-Geddes, P. Zhou, R. Briskine, N. D. Young et al., 2014 Genomic signature of adaptation to climate in Medicago truncatula. Genetics 196: 1263–1275.
