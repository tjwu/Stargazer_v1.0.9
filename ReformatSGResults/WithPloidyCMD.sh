#/stornext/hgsccl/next-gen/Dragen_Merge/SIA0000481/2019-10-31T163146_ILWGS_CLVALQC_SIA0000481_12431_2-FLOWCELL-HT3KVDSXX-HT3WYDSXX-HT5C2DSXX-HTVYJDSXX//BCM_NA17290_A2_SIA0000481.hard-filtered.vcf.gz

VCF=$1

#SN=`echo $VCF | cut -d/ -f6`

SN=`basename $VCF | cut -d. -f1 | awk -F'_' '{print $NF}'`
OrigN=`basename $VCF | cut -d. -f1 | awk -F'_' '{print $2}'`
Folder=`dirname $VCF`
WgsPloidy=`ls ${Folder}/*wgs_ploidy.csv`


python ReformatSGResults_withPloidy.py -dir /hgsccl/next-gen/Illumina/working_area/TJ_working/StargazerVCF_CombinedIntersectPreprocessedLiftover/${OrigN}_${SN} --outfile /hgsccl/next-gen/Illumina/working_area/TJ_working/StargazerVCF_CombinedIntersectPreprocessedLiftover/${OrigN}_${SN}/${OrigN}_${SN}_report.txt --ploidyfile $WgsPloidy --aouvar aou_variants.txt --starshared starall_shared.txt --starlook star_lookup.txt --starfunc star_function.txt

