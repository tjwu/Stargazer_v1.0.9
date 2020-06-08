A note about Combined Annotation Dependent Depletion (CADD) scores:

CADD is a tool for scoring the deleteriousness of SNVs as well as indels in the human genome (PMID: 24487276, 30371827). There are two types of CADD scores: raw and scaled (PHRED). That is, variants at the 10th-% of raw CADD scores are assigned to CADD-10, top 1% to CADD-20, top 0.1% to CADD-30, etc. Stargazer currently uses scaled CADD scores from the GRCh37-v1.4 version. For details, please visit the CADD website (https://cadd.gs.washington.edu).


A note about variant impact in the SNP table:

Each SNV or indel is assigned variant impact (low, moderate, high). High impact means the variant is known to change enzyme function and/or cause abnormal phenotype. Additionally, a variant can be assigned as high impact without activity or phenotype data if it is predicted to be LoF (e.g. frameshift, stop gain, splice defect) or if the CADD score is >= 30 (i.e. among the 0.1% most deleterious variants). Moderate impact means the variant changes protein coding without knowing its effect in activity or phenotype. Finally, low impact means the variant does not affect protein coding (e.g. intron, 3’-UTR). Variants with low impact are most often used to “tag” the high- and moderate-impact variants.
