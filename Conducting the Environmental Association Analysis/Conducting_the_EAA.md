# Environmental Association Analysis
## Background
Many approaches have been developed to explore the importance of genetic local adaptation. Many of these have been based on allele frequency (e.g. Wright 1949; Lewontin and Krakauer 1973; Beaumont and Balding 2004; Beaumont 2005). Other methods have focused on exploring changes in allele frequency gradients (Yang et al., 2012), correlating allele frequencies with environmental variables (e.g. Gunther and Coop, 2013) or use mixed model approaches to identify an association between genetic variants and environmental variables (Eckert et al. 2010; Pyhäjärvi et al. 2013; Yoder et al. 2014).
## Mixed-Model Approach
Mixed-model association can be implemented in many different software packages. The general form of this model is 
Y = **Q**w + **β**x + **Z**u + e⁠, where y is a vector of bioclimatic or soil data, **Q** the incidence matrix of PCs (if modeling population structure), w the vector of fixed PC effects, X the SNP marker genotype matrix, **β** the fixed SNP effect being tested, **Z** the incidence matrix of individuals, **u** ~MVN(0,K $\sigma_u^2$), the vector of random polygenic background effects, and **e** ~MVN(0,K $\sigma_e^2$) the vector of random residuals. This is a flexible model structure that can provide multple different ways of optimizing the model. This can include covariates such as longitude, latitute, maturity group or altitude. This can include additonal metrics to control for population structure (e.g. Kinship matrix, structure assignemnt). 

## Model Testing
Best practices suggest that to ensure the model is not over or underfit by formally fitting the models with different types of covariates and comparing different versions. Models can be compared visually using QQ-plots. The genomic inflation parameter $\lambda$ (Van der most et la., 2017) can also be calculated for each model trait comparison. Addtionally a multiple comparision test correction (e.g. bonferonni) or a false discovery rate (FDR) should be used. 

- Beaumont, M. A., 2005 Adaptation and speciation: what can Fst tell us? Trends Ecol. Evol. 20: 435–440. 
- Beaumont, M. A., D. J. Balding, 2004 Identifying adaptive genetic divergence among populations from genome scans. Mol. Ecol. 13: 969–980.
- Eckert, A. J., J. van Heerwaarden, J. L. Wegrzyn, C. D. Nelson, J. Ross-Ibarra et al., 2010 Patterns of population structure and environmental associations to aridity across the range of loblolly pine (Pinus taeda L., Pinaceae). Genetics 185: 969–982.
- Günther, T., Coop, G. 2013. Robust identification of local adaptation from allele frequencies. Genetics, 195(1), 205-220.
- Lewontin, R. C., and J. Krakauer, 1973 Distribution of gene frequency as a test of the theory of the selective neutrality of polymorphisms. Genetics 74: 175–195.
- Pyhäjärvi, T., M. B. Hufford, S. Mezmouk, and J. Ross-Ibarra, 2013 Complex patterns of local adaptation in teosinte. Genome Biol. Evol. 5: 1594–1609.
- Van Der Most PJ, Küpers LK, Snieder H, Nolte I. QCEWAS: automated quality control of results of epigenome-wide association studies. Bioinformatics. 2017;33(8):1243–1245.
- Wright, S. 1949. The genetical structure of populations. Ann. Eugen. 15: 323–354.
- Yang, W. Y., J. Novembre, E. Eskin, and E. Halperin, 2012 A model-based approach for analysis of spatial structure in genetic data. Nat. Genet. 44: 725–731. 
- Yoder, J. B., J. Stanton-Geddes, P. Zhou, R. Briskine, N. D. Young et al., 2014 Genomic signature of adaptation to climate in Medicago truncatula. Genetics 196: 1263–1275.
