# https://www.broadinstitute.org/gatk/guide/best-practices?bpm=DNAseq
# ftp://ftp.broadinstitute.org/bundle/2.8/b37/
BUILD=Homo_sapiens.GRCh37.75
BUNDLE=/raw/human_sequence_reference/broad_bundle_2.8
REFDIR=/raw/human_sequence_reference
GTK=/opt/GenomeAnalysisTK.jar

cat $REFDIR/$BUILD.dna.chromosome.1.fa \
$REFDIR/$BUILD.dna.chromosome.2.fa \
$REFDIR/$BUILD.dna.chromosome.3.fa \
$REFDIR/$BUILD.dna.chromosome.4.fa \
$REFDIR/$BUILD.dna.chromosome.5.fa \
$REFDIR/$BUILD.dna.chromosome.6.fa \
$REFDIR/$BUILD.dna.chromosome.7.fa \
$REFDIR/$BUILD.dna.chromosome.8.fa \
$REFDIR/$BUILD.dna.chromosome.9.fa \
$REFDIR/$BUILD.dna.chromosome.10.fa \
$REFDIR/$BUILD.dna.chromosome.11.fa \
$REFDIR/$BUILD.dna.chromosome.12.fa \
$REFDIR/$BUILD.dna.chromosome.13.fa \
$REFDIR/$BUILD.dna.chromosome.14.fa \
$REFDIR/$BUILD.dna.chromosome.15.fa \
$REFDIR/$BUILD.dna.chromosome.16.fa \
$REFDIR/$BUILD.dna.chromosome.17.fa \
$REFDIR/$BUILD.dna.chromosome.18.fa \
$REFDIR/$REFDIR/$BUILD.dna.chromosome.19.fa \
$REFDIR/$BUILD.dna.chromosome.20.fa \
$REFDIR/$BUILD.dna.chromosome.21.fa \
$REFDIR/$BUILD.dna.chromosome.22.fa \
$REFDIR/$BUILD.dna.chromosome.X.fa \
$REFDIR/$BUILD.dna.chromosome.Y.fa \
$REFDIR/$BUILD.dna.chromosome.MT.fa > $BUILD.fa

bwa index -a bwtsw $REFDIR/$BUILD.fa
samtools faidx $REFDIR/$BUILD.fa 
java -jar /opt/picard/dist/picard.jar CreateSequenceDictionary \
    REFERENCE=$REFDIR/$BUILD.fa OUTPUT=$REFDIR/$BUILD.dict


java -Xmx4g -jar $GTK -T RealignerTargetCreator \
  -R $REFDIR/$BUILD.fa \
  -o $REFDIR/$BUILD.fa.intervals \
  -known $BUNDLE/Mills_and_1000G_gold_standard.indels.b37.vcf
