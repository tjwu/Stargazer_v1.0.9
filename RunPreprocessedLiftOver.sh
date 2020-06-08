
#python3 ../../../../stargazer.py genotype -o NUDT15_SIA0000387 -t NUDT15 --vcf /stornext/hgsccl/next-gen/Dragen_Merge/SIA0000387/2019-10-31T163146_ILWGS_CLVALQC_SIA0000387_12475_2-FLOWCELL-HT3YHDSXX-HTHGGDSXX-HTHMKDSXX-HTVYJDSXX/BCM_NA18526_C3_SIA0000387.hard-filtered.vcf.gz --data wgs

#/stornext/hgsccl/next-gen/LiftOverResults/BCM_NA17290_A2_SIA0000481_liftover.hard-filtered.vcf.gz
#/stornext/hgsccl/next-gen/LiftOverResults/PreprocessingLiftOver/Liftover_Preprocessed_BCM_NA12878_NA12878.hard-filtered.vcf

VCF=$1
SN=`echo $VCF | awk -F'/' '{print $NF}' | cut -d. -f1 | awk -F'_' '{print $NF}'`
#OrigN=`echo $VCF | awk -F'/' '{print $NF}' | cut -d. -f1 | awk -F'_' '{print $3}'`
OrigN=`basename $VCF | cut -d. -f1 | awk -F'_' '{print $5}'`

for Gene in `cat /hgsccl_software/devel/TJ/Stargazer_v1.0.9/PGxList`
do
    if [ -d /hgsccl/next-gen/Illumina/working_area/TJ_working/StargazerVCF_CombinedIntersectPreprocessedLiftover/${OrigN}_${SN} ]
        then
        echo "" 
    else    
        mkdir /hgsccl/next-gen/Illumina/working_area/TJ_working/StargazerVCF_CombinedIntersectPreprocessedLiftover/${OrigN}_${SN}
    fi
    cd /hgsccl/next-gen/Illumina/working_area/TJ_working/StargazerVCF_CombinedIntersectPreprocessedLiftover/${OrigN}_${SN}
    if [ -e /hgsccl/next-gen/Illumina/working_area/TJ_working/StargazerVCF_CombinedIntersectPreprocessedLiftover/${OrigN}_${SN}/${Gene}_${OrigN}_${SN}.stargazer-genotype.txt ]
    then
        continue
    else
        #echo /hgsccl/next-gen/Illumina/working_area/TJ_working/StargazerVCF_PreprocessedLiftover/${OrigN}_${SN}/${Gene}_${OrigN}_${SN}.txt 
        echo "source /hgsccl_software/devel/TJ/Stargazer_v1.0.9/SetEnv.sh; 
        python3 /hgsccl_software/devel/TJ/Stargazer_v1.0.9/stargazer.py genotype -o ${Gene}_${OrigN}_${SN} -t $Gene --vcf $VCF --data wgs --output_dir /hgsccl/next-gen/Illumina/working_area/TJ_working/StargazerVCF_CombinedIntersectPreprocessedLiftover/${OrigN}_${SN} "  | msub -q hgsccl -l mem=2gb -N Stargazer_${Gene}_${OrigN}_${SN}
    fi
done
