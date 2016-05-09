#!/bin/bash

# Example call:
# bash /notebook/code/src/sequencing/neoepitope/predict_neoepitopes_from_VCF.sh \
# --input_path=${DIR_PRE}/Pre_prembo_ipi.passed.somatic.indels.dbsnp.snpeff.vcf \
# --output_folder=/notebook/human_sequencing_collison/P0207_pre_post_immunotherapy/pre \
# --root_output_filename=Pre_prembo_ipi.passed.somatic.indels.dbsnp.snpeff \
# --bedtools_bin_path=/opt/bedtools2/bin \
# --netmhc_path=/opt/netMHC-3.4 \
# --annovar_path=/opt/annovar \
# --gtf_path=/opt/annovar/humandb/hg19_refGene.txt \
# --generate_exon_bed_file=/notebook/code/src/sequencing/neoepitope/generate_exon_bed_file.py \
# --predict_neo_path=/notebook/code/src/sequencing/neoepitope/extract_neoepitopes_from_annovar_coding_change.py


# parse arguments
#----------------

for i in "$@"
do
case $i in
    -i=*|--input_path=*)
    VCF_PATH_RAW="${i#*=}"
    shift # past argument=value
    ;;
    -o=*|--output_folder=*)
    OUTPUT_PATH="${i#*=}"
    shift # past argument=value
    ;;
    -a=*|--annovar_path=*)
    ANNOVAR_PATH="${i#*=}"
    shift # past argument=value
    ;;
    -b=*|--bedtools_bin_path=*)
    BEDTOOLS_PATH="${i#*=}"
    shift # past argument=value
    ;;
    -n=*|--netmhc_path=*)
    NETMHC_PATH="${i#*=}"
    shift # past argument=value
    ;; 
    -g=*|--gtf_path=*)
    GTF_PATH="${i#*=}"
    shift # past argument=value
    ;;
    -p=*|--predict_neo_path=*)
    PREDICT_PATH="${i#*=}"
    shift # past argument=value
    ;;
    -e=*|--generate_exon_bed_file=*)
    EXON_BED_SCRIPT_PATH="${i#*=}"
    shift # past argument=value
    ;;    
    -r=*|--root_output_filename=*)
    ROOT_NAME="${i#*=}"
    shift # past argument=value
    ;;    
    --default)
    DEFAULT=YES
    shift # past argument with no value
    ;;
    *)
    # unknown option
    ;;
esac
done

if [ -f "${VCF_PATH_RAW}" ]
then
	echo "Parsing            ${VCF_PATH_RAW}"
else
	echo "ERROR: could not open input VCF file: ${VCF_PATH_RAW}"
	exit 1
fi
if [ -d "$OUTPUT_PATH" ]
then
	echo "Writing output:    ${OUTPUT_PATH}"
else
	echo "ERROR: could not find output folder ${OUTPUT_PATH}"
	exit 1
fi
if [ -d "$ANNOVAR_PATH" ]
then
	echo "ANNOVAR installed: ${ANNOVAR_PATH}"
else
	echo "ERROR: could not find output folder ${ANNOVAR_PATH}"
	exit 1
fi
if [ -f "${BEDTOOLS_PATH}/intersectBed" ]
then
    echo "Bedtools folder:   ${BEDTOOLS_PATH}"
else
    echo "ERROR: could not find bedtools folder at ${BEDTOOLS_PATH}"
	exit 1
fi
if [ -f "${NETMHC_PATH}/netMHC.py" ]
then
    echo "Netmhc folder:     ${NETMHC_PATH}"
    export NMHOME=${NETMHC_PATH}
else
    echo "ERROR: could not find netmhc at $NETHMC_PATH"
	exit 1
fi
if [ -f "${PREDICT_PATH}" ]
then
    echo "Neoepitope script: ${PREDICT_PATH}"
else
    echo "ERROR: could not find neoepitope prediction script at ${PREDICT_PATH}"
	exit 1
fi
if [ -f "${GTF_PATH}" ]
then
    echo "GTF for exons:     ${GTF_PATH}"
else
    echo "ERROR: could not find GTF at ${GTF_PATH}"
	exit 1
fi

export TMPDIR=/tmp

echo "------------------------------------------------------------"
   
# Restrict VCFs to lie within REFSEQ regions (although not necessarily within exons)
# create BED from the annovar GFF
# remove chr from the BED we created, remove that file
#
# Strip .vcf from input file to make subsequent filenames clean 
#

OUT_BASE=${OUTPUT_PATH}/${ROOT_NAME}
OUTPUT_BED_PATH_CHR=${OUTPUT_PATH}/hg19_refseq_chr.bed
OUTPUT_BED_PATH=${OUTPUT_PATH}/hg19_refseq.bed
OUTPUT_BED_EXON_PATH=${OUTPUT_PATH}/hg19_refseq_exons.bed

# improve to autodetect whether chr is in vcf or not
awk '{print($3"\t"$5"\t"$6"\t"$2)}' $GTF_PATH > ${OUTPUT_BED_PATH_CHR}
cp ${OUTPUT_BED_PATH_CHR} ${OUTPUT_BED_PATH}
#awk '{sub(/chr/,"")}; 1' ${OUTPUT_BED_PATH_CHR} > ${OUTPUT_BED_PATH}
#rm ${OUTPUT_BED_PATH_CHR}

${BEDTOOLS_PATH}/intersectBed -a ${VCF_PATH_RAW} -b ${OUTPUT_BED_PATH} \
  -u -header > ${OUT_BASE}_cleaned_inRS.vcf

# Generate exon BED file from the hg19_refseq gene file
#
python ${EXON_BED_SCRIPT_PATH} -i ${GTF_PATH} -o ${OUTPUT_BED_EXON_PATH}

# extract variants that intersect an exon (excluding introns)
#
${BEDTOOLS_PATH}/intersectBed -a ${OUT_BASE}_cleaned_inRS.vcf -b ${OUTPUT_BED_EXON_PATH} \
  -u -header > ${OUT_BASE}_cleaned_inRS_exons.vcf

# Annotate cleaned VCFs with ANNOVAR
#
${ANNOVAR_PATH}/table_annovar.pl ${OUT_BASE}_cleaned_inRS.vcf \
  ${ANNOVAR_PATH}/humandb/ -buildver hg19 -out ${OUT_BASE} \
  -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,snp138,ljb26_all \
  -operation g,r,r,f,f,f,f -nastring . -vcfinput 

# Predict protein products that result from germline or mutant DNA
# 
${ANNOVAR_PATH}/coding_change.pl --includesnp \
  ${OUT_BASE}.refGene.exonic_variant_function \
  ${GTF_PATH} ${ANNOVAR_PATH}/humandb/hg19_refGeneMrna.fa > ${OUT_BASE}_protein_changes.txt


python $PREDICT_PATH \
  -i ${OUT_BASE}_protein_changes.txt \
  -v ${OUT_BASE}.refGene.exonic_variant_function \
  -g ${GTF_PATH} \
  -n 9 \
  -w ${OUT_BASE}_neoantigen_wildtype.fa \
  -m ${OUT_BASE}_neoantigen_mutant.fa

python ${NETMHC_PATH}/netMHC.py -a HLA-A01:01 ${OUT_BASE}_neoantigen_wildtype.fa > ${OUT_BASE}_NetMHC_HLA-A01_WT.txt
python ${NETMHC_PATH}/netMHC.py -a HLA-A01:01 ${OUT_BASE}_neoantigen_mutant.fa > ${OUT_BASE}_NetMHC_HLA-A01_mut.txt
python ${NETMHC_PATH}/netMHC.py -a HLA-A02:01 ${OUT_BASE}_neoantigen_wildtype.fa > ${OUT_BASE}_NetMHC_HLA-A02_WT.txt
python ${NETMHC_PATH}/netMHC.py -a HLA-A02:01 ${OUT_BASE}_neoantigen_mutant.fa > ${OUT_BASE}_NetMHC_HLA-A02_mut.txt

# Consolidate all HLA results into a single file
#
timestamp(){
    date +"%Y-%m-%d_%H-%M-%S"
}
TIME=$(timestamp)
echo "#created:${TIME}" > ${OUT_BASE}_NetMHC.txt
echo "#raw:${VCF_PATH_RAW}" >> ${OUT_BASE}_NetMHC.txt
echo "#output:${OUTPUT_PATH}" >> ${OUT_BASE}_NetMHC.txt
echo "#annovar:${ANNOVAR_PATH}" >> ${OUT_BASE}_NetMHC.txt
echo "#bedtools:${BEDTOOLS_PATH}" >> ${OUT_BASE}_NetMHC.txt
echo "#netmhc:${NETMHC_PATH}" >> ${OUT_BASE}_NetMHC.txt
echo "#Neoepitope_script:${PREDICT_PATH}" >> ${OUT_BASE}_NetMHC.txt
echo "#GTF:${GTF_PATH}" >> ${OUT_BASE}_NetMHC.txt
tail -n +4 ${OUT_BASE}_NetMHC_HLA-A01_WT.txt >> ${OUT_BASE}_NetMHC.txt
tail -n +5 ${OUT_BASE}_NetMHC_HLA-A01_mut.txt >> ${OUT_BASE}_NetMHC.txt
tail -n +5 ${OUT_BASE}_NetMHC_HLA-A02_WT.txt >> ${OUT_BASE}_NetMHC.txt
tail -n +5 ${OUT_BASE}_NetMHC_HLA-A02_mut.txt >> ${OUT_BASE}_NetMHC.txt

echo "Final results: ${OUT_BASE}_NetMHC.txt"
