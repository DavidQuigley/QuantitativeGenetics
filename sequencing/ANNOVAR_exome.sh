# /notebook/code/src/sequencing/annotation/ANNOVAR_exome.sh \
#   -i=/raw/human_lines_cholangio/2016_05_05_exomes/GATK.C69.Recalibration.vcf \
#   -o=/datasets/human_lines_cholangio/2016_05_05_exomes/annovar \
#   -a=/opt/annovar \
#   -b=/opt/bedtools2/bin \
#   -t=/datasets/human_sequence_reference/annovar/humandb
#    -e=/notebook/code/src/sequencing/neoepitope/generate_exon_bed_file.py
#   -r=cholangio
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
    -t=*|--output_folder=*)
    ANNOVAR_ANNOTATIONS_PATH="${i#*=}"
    shift # past argument=value
    ;;
    -b=*|--bedtools_bin_path=*)
    BEDTOOLS_PATH="${i#*=}"
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

echo "------------------------------------------------------------"


OUT_BASE=${OUTPUT_PATH}/${ROOT_NAME}
OUTPUT_BED_PATH_CHR=${OUTPUT_PATH}/hg19_refseq_chr.bed
OUTPUT_BED_PATH=${OUTPUT_PATH}/hg19_refseq.bed
OUTPUT_BED_EXON_PATH=${OUTPUT_PATH}/hg19_refseq_exons.bed
GTF_PATH=${ANNOVAR_ANNOTATIONS_PATH}/hg19_refGene.txt

echo "Converting GTF to BED format, lacking 'chr' at beginning of line"
awk '{print($3"\t"$5"\t"$6"\t"$2)}' $GTF_PATH > ${OUTPUT_BED_PATH_CHR}
cp ${OUTPUT_BED_PATH_CHR} ${OUTPUT_BED_PATH}
awk '{sub(/chr/,"")}; 1' ${OUTPUT_BED_PATH_CHR} > ${OUTPUT_BED_PATH}
rm ${OUTPUT_BED_PATH_CHR}

echo "VCF to exons present in GTF"
echo "Bedtools intersectBed to restrict VCF to known exons"
${BEDTOOLS_PATH}/intersectBed -a ${VCF_PATH_RAW} -b ${OUTPUT_BED_PATH} \
  -u -header > ${OUT_BASE}_cleaned_inRS.vcf

# Generate exon BED file from annovar's hg19_refseq gene file
#
python ${EXON_BED_SCRIPT_PATH} -i ${GTF_PATH} -o ${OUTPUT_BED_EXON_PATH}

# extract variants that intersect an exon (excluding introns)
#
${BEDTOOLS_PATH}/intersectBed -a ${OUT_BASE}_cleaned_inRS.vcf \
                              -b ${OUTPUT_BED_EXON_PATH} \
                              -u \
                              -header > ${OUT_BASE}_cleaned_inRS_exons.vcf

echo "Beginning annotation"
# Annotate cleaned VCFs with ANNOVAR
#
${ANNOVAR_PATH}/table_annovar.pl ${OUT_BASE}_cleaned_inRS.vcf \
  ${ANNOVAR_ANNOTATIONS_PATH}/ --buildver hg19 -out ${OUT_BASE} \
  -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,snp138,ljb26_all,cosmic70,clinvar_20160302 \
  -operation g,r,r,f,f,f,f,f,f -nastring . -vcfinput 
