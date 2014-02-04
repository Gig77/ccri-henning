export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files

MERGE = ~/tools/samtools-0.1.19/samtools merge -f $@.part $^; mv $@.part $@; rm $^
LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

all: MHH_115 MHH_121 MHH_126 MHH_130 MHH_133 MHH_136 MHH_139

filtered-variants.tsv: filtered-variants/MHH_115.tsv filtered-variants/MHH_121.tsv filtered-variants/MHH_126.tsv filtered-variants/MHH_130.tsv filtered-variants/MHH_133.tsv filtered-variants/MHH_136.tsv filtered-variants/MHH_139.tsv 
	perl ~/henning/scripts/filter-variants.pl --header 1 > $@.part
	cat $^ >> $@.part
	mv $@.part $@

TEST: bam/TEST.bwa.merged.dup.bam bam/TEST.bwa.merged.dup.bam.bai picard/TEST.multiplemetrics picard/TEST.hs_metrics
TEST_BAM = bam/TEST_s_1.bwa.bam bam/TEST_s_2.bwa.bam
bam/TEST.bwa.merged.bam: $(TEST_BAM)
	$(MERGE)

MHH_115: varscan/MHH_115.varscan.snp.dbsnp.snpeff.vcf MHH_115_tumor MHH_115_normal
MHH_115_tumor: picard/MHH_115.tumor.multiplemetrics picard/MHH_115.tumor.hs_metrics 
MHH_115_normal: picard/MHH_115.normal.multiplemetrics picard/MHH_115.normal.hs_metrics
bam/MHH_115.tumor.bwa.merged.bam: bam/AML_015_s_1.bwa.bam bam/AML_015_s_2.bwa.bam bam/AML_015_s_3.bwa.bam bam/AML_015_s_4.bwa.bam bam/AML_015_s_5.bwa.bam bam/AML_015_s_6.bwa.bam bam/AML_015_s_7.bwa.bam bam/AML_015_s_8.bwa.bam
	$(MERGE)
bam/MHH_115.normal.bwa.merged.bam: bam/AML_016_s_1.bwa.bam bam/AML_016_s_2.bwa.bam bam/AML_016_s_3.bwa.bam bam/AML_016_s_4.bwa.bam bam/AML_016_s_5.bwa.bam bam/AML_016_s_6.bwa.bam bam/AML_016_s_7.bwa.bam bam/AML_016_s_8.bwa.bam
	$(MERGE)

MHH_121: varscan/MHH_121.varscan.snp.dbsnp.snpeff.vcf MHH_121_tumor MHH_121_normal
MHH_121_tumor: bam/MHH_121.tumor.bwa.merged.dup.bam picard/MHH_121.tumor.multiplemetrics picard/MHH_121.tumor.hs_metrics 
MHH_121_normal: bam/MHH_121.normal.bwa.merged.dup.bam picard/MHH_121.normal.multiplemetrics picard/MHH_121.normal.hs_metrics
bam/MHH_121.tumor.bwa.merged.bam: bam/AML_021_s_1.bwa.bam bam/AML_021_s_2.bwa.bam bam/AML_021_s_3.bwa.bam bam/AML_021_s_4.bwa.bam bam/AML_021_s_5.bwa.bam bam/AML_021_s_6.bwa.bam bam/AML_021_s_7.bwa.bam bam/AML_021_s_8.bwa.bam
	$(MERGE)
bam/MHH_121.normal.bwa.merged.bam: bam/AML_022_s_1.bwa.bam bam/AML_022_s_2.bwa.bam bam/AML_022_s_3.bwa.bam bam/AML_022_s_4.bwa.bam bam/AML_022_s_5.bwa.bam bam/AML_022_s_6.bwa.bam bam/AML_022_s_7.bwa.bam bam/AML_022_s_8.bwa.bam
	$(MERGE)

MHH_126: varscan/MHH_126.varscan.snp.dbsnp.snpeff.vcf MHH_126_tumor MHH_126_normal
MHH_126_tumor: picard/MHH_126.tumor.multiplemetrics picard/MHH_126.tumor.hs_metrics 
MHH_126_normal: picard/MHH_126.normal.multiplemetrics picard/MHH_126.normal.hs_metrics
bam/MHH_126.tumor.bwa.merged.bam: bam/AML_026_s_1.bwa.bam bam/AML_026_s_2.bwa.bam bam/AML_026_s_3.bwa.bam bam/AML_026_s_4.bwa.bam bam/AML_026_s_5.bwa.bam bam/AML_026_s_6.bwa.bam bam/AML_026_s_7.bwa.bam bam/AML_026_s_8.bwa.bam
	$(MERGE)
bam/MHH_126.normal.bwa.merged.bam: bam/AML_028_s_1.bwa.bam bam/AML_028_s_2.bwa.bam bam/AML_028_s_3.bwa.bam bam/AML_028_s_4.bwa.bam bam/AML_028_s_5.bwa.bam bam/AML_028_s_6.bwa.bam bam/AML_028_s_7.bwa.bam bam/AML_028_s_8.bwa.bam
	$(MERGE)

MHH_130: varscan/MHH_130.varscan.snp.dbsnp.snpeff.vcf MHH_130_tumor MHH_130_normal
MHH_130_tumor: picard/MHH_130.tumor.multiplemetrics picard/MHH_130.tumor.hs_metrics 
MHH_130_normal: picard/MHH_130.normal.multiplemetrics picard/MHH_130.normal.hs_metrics
bam/MHH_130.tumor.bwa.merged.bam: bam/AML_030_s_1.bwa.bam bam/AML_030_s_2.bwa.bam bam/AML_030_s_3.bwa.bam bam/AML_030_s_4.bwa.bam bam/AML_030_s_5.bwa.bam bam/AML_030_s_6.bwa.bam bam/AML_030_s_7.bwa.bam bam/AML_030_s_8.bwa.bam
	$(MERGE)
bam/MHH_130.normal.bwa.merged.bam: bam/AML_031_s_1.bwa.bam bam/AML_031_s_2.bwa.bam bam/AML_031_s_3.bwa.bam bam/AML_031_s_4.bwa.bam bam/AML_031_s_5.bwa.bam bam/AML_031_s_6.bwa.bam bam/AML_031_s_7.bwa.bam bam/AML_031_s_8.bwa.bam
	$(MERGE)

MHH_133: varscan/MHH_133.varscan.snp.dbsnp.snpeff.vcf MHH_133_tumor MHH_133_normal
MHH_133_tumor: picard/MHH_133.tumor.multiplemetrics picard/MHH_133.tumor.hs_metrics 
MHH_133_normal: picard/MHH_133.normal.multiplemetrics picard/MHH_133.normal.hs_metrics
bam/MHH_133.tumor.bwa.merged.bam: bam/AML_033_s_1.bwa.bam bam/AML_033_s_2.bwa.bam bam/AML_033_s_3.bwa.bam bam/AML_033_s_4.bwa.bam bam/AML_033_s_5.bwa.bam bam/AML_033_s_6.bwa.bam bam/AML_033_s_7.bwa.bam bam/AML_033_s_8.bwa.bam
	$(MERGE)
bam/MHH_133.normal.bwa.merged.bam: bam/AML_034_s_1.bwa.bam bam/AML_034_s_2.bwa.bam bam/AML_034_s_3.bwa.bam bam/AML_034_s_4.bwa.bam bam/AML_034_s_5.bwa.bam bam/AML_034_s_6.bwa.bam bam/AML_034_s_7.bwa.bam bam/AML_034_s_8.bwa.bam
	$(MERGE)

MHH_136: varscan/MHH_136.varscan.snp.dbsnp.snpeff.vcf MHH_136_tumor MHH_136_normal
MHH_136_tumor: picard/MHH_136.tumor.multiplemetrics picard/MHH_136.tumor.hs_metrics 
MHH_136_normal: picard/MHH_136.normal.multiplemetrics picard/MHH_136.normal.hs_metrics
bam/MHH_136.tumor.bwa.merged.bam: bam/AML_036_s_1.bwa.bam bam/AML_036_s_2.bwa.bam bam/AML_036_s_3.bwa.bam bam/AML_036_s_4.bwa.bam bam/AML_036_s_5.bwa.bam bam/AML_036_s_6.bwa.bam bam/AML_036_s_7.bwa.bam bam/AML_036_s_8.bwa.bam
	$(MERGE)
bam/MHH_136.normal.bwa.merged.bam: bam/AML_037_s_1.bwa.bam bam/AML_037_s_2.bwa.bam bam/AML_037_s_3.bwa.bam bam/AML_037_s_4.bwa.bam bam/AML_037_s_5.bwa.bam bam/AML_037_s_6.bwa.bam bam/AML_037_s_7.bwa.bam bam/AML_037_s_8.bwa.bam
	$(MERGE)

MHH_139: varscan/MHH_139.varscan.snp.dbsnp.snpeff.vcf MHH_139_tumor MHH_139_normal
MHH_139_tumor: picard/MHH_139.tumor.multiplemetrics picard/MHH_139.tumor.hs_metrics 
MHH_139_normal: picard/MHH_139.normal.multiplemetrics picard/MHH_139.normal.hs_metrics
bam/MHH_139.tumor.bwa.merged.bam: bam/AML_039_s_1.bwa.bam bam/AML_039_s_2.bwa.bam bam/AML_039_s_3.bwa.bam bam/AML_039_s_4.bwa.bam bam/AML_039_s_5.bwa.bam bam/AML_039_s_6.bwa.bam bam/AML_039_s_7.bwa.bam bam/AML_039_s_8.bwa.bam
	$(MERGE)
bam/MHH_139.normal.bwa.merged.bam: bam/AML_040_s_1.bwa.bam bam/AML_040_s_2.bwa.bam bam/AML_040_s_3.bwa.bam bam/AML_040_s_4.bwa.bam bam/AML_040_s_5.bwa.bam bam/AML_040_s_6.bwa.bam bam/AML_040_s_7.bwa.bam bam/AML_040_s_8.bwa.bam
	$(MERGE)
 
#-----------	
# PREPARE FASTQ INPUT
#-----------	
~/henning/data/fastq/TEST_s_1_1_sequence.txt.gz: ~/henning/data/fastq/AML_022_s_1_1_sequence.txt.gz ~/henning/data/fastq/AML_022_s_1_2_sequence.txt.gz ~/henning/data/fastq/AML_022_s_2_1_sequence.txt.gz ~/henning/data/fastq/AML_022_s_2_2_sequence.txt.gz
	~/tools/seqtk/seqtk sample -s100 ~/henning/data/fastq/AML_022_s_1_1_sequence.txt.gz 10000 | gzip > ~/henning/data/fastq/TEST_s_1_1_sequence.txt.gz
	~/tools/seqtk/seqtk sample -s100 ~/henning/data/fastq/AML_022_s_1_2_sequence.txt.gz 10000 | gzip > ~/henning/data/fastq/TEST_s_1_2_sequence.txt.gz
	~/tools/seqtk/seqtk sample -s324 ~/henning/data/fastq/AML_022_s_2_1_sequence.txt.gz 10000 | gzip > ~/henning/data/fastq/TEST_s_2_1_sequence.txt.gz
	~/tools/seqtk/seqtk sample -s324 ~/henning/data/fastq/AML_022_s_2_2_sequence.txt.gz 10000 | gzip > ~/henning/data/fastq/TEST_s_2_2_sequence.txt.gz

AML_015_input: AML_015_L2 AML_015_L3 AML_015_L4 AML_015_L5 AML_015_L6 AML_015_L7 AML_015_L8
AML_015_L%: ~/henning/data/fastq/B04GLABXX.%.AML_015.unmapped.bam
	#curl http://www.biomedical-sequencing.at/DNA/Pediatric_AML_Sequencing/AML_015/B04GLABXX.$*.AML_015.unmapped.bam -so ~/henning/data/fastq/B04GLABXX.$*.AML_015.unmapped.bam
	java -Xmx2g -jar ~/tools/picard-tools-1.101/SamToFastq.jar I=~/henning/data/fastq/B04GLABXX.$*.AML_015.unmapped.bam F=~/henning/data/fastq/AML_015_s_$*_1_sequence.txt F2=~/henning/data/fastq/AML_015_s_$*_2_sequence.txt
	gzip ~/henning/data/fastq/AML_015_s_$*_1_sequence.txt ~/henning/data/fastq/AML_015_s_$*_2_sequence.txt
AML_016_input: AML_016_L1 AML_016_L2 AML_016_L3 AML_016_L4 AML_016_L5 AML_016_L6 AML_016_L7 AML_016_L8
AML_016_L%: ~/henning/data/fastq/B04GLABXX.%.AML_016.unmapped.bam
	java -Xmx2g -jar ~/tools/picard-tools-1.101/SamToFastq.jar I=~/henning/data/fastq/B04GLABXX.$*.AML_016.unmapped.bam F=~/henning/data/fastq/AML_016_s_$*_1_sequence.txt F2=~/henning/data/fastq/AML_016_s_$*_2_sequence.txt
	gzip ~/henning/data/fastq/AML_016_s_$*_1_sequence.txt ~/henning/data/fastq/AML_016_s_$*_2_sequence.txt
AML_026_input: AML_026_L1 AML_026_L2 AML_026_L3 AML_026_L4 AML_026_L5 AML_026_L6 AML_026_L7 AML_026_L8
AML_026_L%: ~/henning/data/fastq/AB0B88ABXX.%.AML_026.unmapped.bam
	java -Xmx2g -jar ~/tools/picard-tools-1.101/SamToFastq.jar I=~/henning/data/fastq/AB0B88ABXX.$*.AML_026.unmapped.bam F=~/henning/data/fastq/AML_026_s_$*_1_sequence.txt F2=~/henning/data/fastq/AML_026_s_$*_2_sequence.txt
	gzip ~/henning/data/fastq/AML_026_s_$*_1_sequence.txt ~/henning/data/fastq/AML_026_s_$*_2_sequence.txt

~/henning/data/fastq/AML_028_s_%_1_sequence.txt.gz ~/henning/data/fastq/AML_028_s_%_2_sequence.txt.gz: ~/henning/data/fastq/AB0B88ABXX.%.AML_028.unmapped.bam
	java -Xmx2g -jar ~/tools/picard-tools-1.101/SamToFastq.jar I=~/henning/data/fastq/AB0B88ABXX.$*.AML_028.unmapped.bam F=~/henning/data/fastq/AML_028_s_$*_1_sequence.txt F2=~/henning/data/fastq/AML_028_s_$*_2_sequence.txt
	gzip ~/henning/data/fastq/AML_028_s_$*_1_sequence.txt ~/henning/data/fastq/AML_028_s_$*_2_sequence.txt
	rm ~/henning/data/fastq/AB0B88ABXX.$*.AML_028.unmapped.bam

~/henning/data/fastq/AML_030_s_%_1_sequence.txt.gz ~/henning/data/fastq/AML_030_s_%_2_sequence.txt.gz: ~/henning/data/fastq/AB0B88ABXX.%.AML_030.unmapped.bam
	java -Xmx2g -jar ~/tools/picard-tools-1.101/SamToFastq.jar I=~/henning/data/fastq/AB0B88ABXX.$*.AML_030.unmapped.bam F=~/henning/data/fastq/AML_030_s_$*_1_sequence.txt F2=~/henning/data/fastq/AML_030_s_$*_2_sequence.txt
	gzip ~/henning/data/fastq/AML_030_s_$*_1_sequence.txt ~/henning/data/fastq/AML_030_s_$*_2_sequence.txt
	rm ~/henning/data/fastq/AB0B88ABXX.$*.AML_030.unmapped.bam

~/henning/data/fastq/AML_031_s_%_1_sequence.txt.gz ~/henning/data/fastq/AML_031_s_%_2_sequence.txt.gz: ~/henning/data/fastq/AB0B88ABXX.%.AML_031.unmapped.bam
	java -Xmx2g -jar ~/tools/picard-tools-1.101/SamToFastq.jar I=~/henning/data/fastq/AB0B88ABXX.$*.AML_031.unmapped.bam F=~/henning/data/fastq/AML_031_s_$*_1_sequence.txt F2=~/henning/data/fastq/AML_031_s_$*_2_sequence.txt
	gzip ~/henning/data/fastq/AML_031_s_$*_1_sequence.txt ~/henning/data/fastq/AML_031_s_$*_2_sequence.txt
	rm ~/henning/data/fastq/AB0B88ABXX.$*.AML_031.unmapped.bam

~/henning/data/fastq/AML_033_s_%_1_sequence.txt.gz ~/henning/data/fastq/AML_033_s_%_2_sequence.txt.gz: ~/henning/data/fastq/B04GLABXX.%.AML_033.unmapped.bam
	java -Xmx2g -jar ~/tools/picard-tools-1.101/SamToFastq.jar I=~/henning/data/fastq/B04GLABXX.$*.AML_033.unmapped.bam F=~/henning/data/fastq/AML_033_s_$*_1_sequence.txt F2=~/henning/data/fastq/AML_033_s_$*_2_sequence.txt
	gzip ~/henning/data/fastq/AML_033_s_$*_1_sequence.txt ~/henning/data/fastq/AML_033_s_$*_2_sequence.txt
	rm ~/henning/data/fastq/B04GLABXX.$*.AML_033.unmapped.bam

~/henning/data/fastq/AML_034_s_%_1_sequence.txt.gz ~/henning/data/fastq/AML_034_s_%_2_sequence.txt.gz: ~/henning/data/fastq/B04GLABXX.%.AML_034.unmapped.bam
	java -Xmx2g -jar ~/tools/picard-tools-1.101/SamToFastq.jar I=~/henning/data/fastq/B04GLABXX.$*.AML_034.unmapped.bam F=~/henning/data/fastq/AML_034_s_$*_1_sequence.txt F2=~/henning/data/fastq/AML_034_s_$*_2_sequence.txt
	gzip ~/henning/data/fastq/AML_034_s_$*_1_sequence.txt ~/henning/data/fastq/AML_034_s_$*_2_sequence.txt
	rm ~/henning/data/fastq/B04GLABXX.$*.AML_034.unmapped.bam

~/henning/data/fastq/AML_036_s_%_1_sequence.txt.gz ~/henning/data/fastq/AML_036_s_%_2_sequence.txt.gz: ~/henning/data/fastq/AB0B88ABXX.%.AML_036.unmapped.bam
	java -Xmx2g -jar ~/tools/picard-tools-1.101/SamToFastq.jar I=~/henning/data/fastq/AB0B88ABXX.$*.AML_036.unmapped.bam F=~/henning/data/fastq/AML_036_s_$*_1_sequence.txt F2=~/henning/data/fastq/AML_036_s_$*_2_sequence.txt
	gzip ~/henning/data/fastq/AML_036_s_$*_1_sequence.txt ~/henning/data/fastq/AML_036_s_$*_2_sequence.txt
	rm ~/henning/data/fastq/AB0B88ABXX.$*.AML_036.unmapped.bam

~/henning/data/fastq/AML_037_s_%_1_sequence.txt.gz ~/henning/data/fastq/AML_037_s_%_2_sequence.txt.gz: ~/henning/data/fastq/AB0B88ABXX.%.AML_037.unmapped.bam
	java -Xmx2g -jar ~/tools/picard-tools-1.101/SamToFastq.jar I=~/henning/data/fastq/AB0B88ABXX.$*.AML_037.unmapped.bam F=~/henning/data/fastq/AML_037_s_$*_1_sequence.txt F2=~/henning/data/fastq/AML_037_s_$*_2_sequence.txt
	gzip ~/henning/data/fastq/AML_037_s_$*_1_sequence.txt ~/henning/data/fastq/AML_037_s_$*_2_sequence.txt
	rm ~/henning/data/fastq/AB0B88ABXX.$*.AML_037.unmapped.bam

#-----------	
# SAMTOOLS, ALIGNMENT
#-----------	
%.sai: %.gz
	~/tools/bwa-0.7.4/bwa aln -t 1 ~/generic/data/hg19/ucsc.hg19.fasta $^ 2>&1 1>$@.part | $(LOG)
	mv $@.part $@

bam/%.bwa.bam: ~/henning/data/fastq/%_1_sequence.txt.sai ~/henning/data/fastq/%_2_sequence.txt.sai ~/henning/data/fastq/%_1_sequence.txt.gz ~/henning/data/fastq/%_2_sequence.txt.gz
	sleep `expr $$RANDOM % 60`
	(~/tools/bwa-0.7.4/bwa sampe ~/generic/data/hg19/ucsc.hg19.fasta $^ \
		| ~/tools/samtools-0.1.19/samtools view -bS -) 2>&1 1>$@.part | $(LOG)
	mv $@.part $@
	~/tools/samtools-0.1.19/samtools sort -m 2G $@ $@.sorted 2>&1 | $(LOG)
	mv $@.sorted.bam $@
	rm ~/henning/data/fastq/$*_1_sequence.txt.sai ~/henning/data/fastq/$*_2_sequence.txt.sai

bam/%.bam.bai: bam/%.bam
	rm -f $@
	~/tools/samtools-0.1.19/samtools index $^ $@.part 2>&1 | $(LOG)
	mv $@.part $@

#-----------	
# PICARD: MARK DUPLICATES, METRICS
#-----------	
bam/%.bwa.merged.dup.bam: bam/%.bwa.merged.bam bam/%.bwa.merged.bam.bai
	java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar ~/tools/picard-tools-1.101/MarkDuplicates.jar \
		INPUT=$< \
		OUTPUT=$@.part \
		METRICS_FILE=picard/$*.mark_duplicates_metrics \
		VALIDATION_STRINGENCY=LENIENT \
		2>&1 | $(LOG)
	mv $@.part $@

picard/%.multiplemetrics: bam/%.bwa.merged.dup.bam bam/%.bwa.merged.dup.bam.bai
	java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar ~/tools/picard-tools-1.101/CollectMultipleMetrics.jar \
		INPUT=$< \
		OUTPUT=picard/$* \
		VALIDATION_STRINGENCY=LENIENT \
		PROGRAM=CollectAlignmentSummaryMetrics \
		PROGRAM=CollectInsertSizeMetrics \
		PROGRAM=QualityScoreDistribution \
		PROGRAM=MeanQualityByCycle \
		2>&1 | $(LOG)
	touch $@

picard/%.hs_metrics: bam/%.bwa.merged.dup.bam bam/%.bwa.merged.dup.bam.bai
	# truseq_exome_targeted_regions.hg19.bed.chr downloaded from http://supportres.illumina.com/documents/myillumina/5dfd7e70-c4a5-405a-8131-33f683414fb7/truseq_exome_targeted_regions.hg19.bed.chr.gz
	~/tools/samtools-0.1.19/samtools view -H bam/$*.bwa.merged.dup.bam 2>&1 1> picard/$*.truseq-for-picard.bed | $(LOG)
	gawk 'BEGIN { OFS="\t"} {print $$1,$$2,$$3,$$6,$$4 }' ~/generic/data/illumina/truseq_exome_targeted_regions.hg19.bed.chr >> picard/$*.truseq-for-picard.bed
	java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar ~/tools/picard-tools-1.101/CalculateHsMetrics.jar \
		BAIT_INTERVALS=picard/$*.truseq-for-picard.bed \
		TARGET_INTERVALS=picard/$*.truseq-for-picard.bed \
		INPUT=bam/$*.bwa.merged.dup.bam \
		OUTPUT=$@.part \
		REFERENCE_SEQUENCE=~/generic/data/hg19/ucsc.hg19.fasta \
		PER_TARGET_COVERAGE=picard/$*.hs_metrics.per_target_coverage.part \
		VALIDATION_STRINGENCY=LENIENT \
		2>&1 | $(LOG)
	mv picard/$*.hs_metrics.per_target_coverage.part picard/$*.hs_metrics.per_target_coverage
	rm picard/$*.truseq-for-picard.bed
	mv $@.part $@

#-----------	
# VARSCAN
#-----------	
varscan/%.varscan.snp.vcf: bam/%.normal.bwa.merged.dup.bam bam/%.tumor.bwa.merged.dup.bam bam/%.normal.bwa.merged.dup.bam.bai bam/%.tumor.bwa.merged.dup.bam.bai
	java -jar ~/tools/varscan-2.3.6/VarScan.v2.3.6.jar somatic \
		<(~/tools/samtools-0.1.19/samtools view -b -u -q 1 $(word 1,$^) | ~/tools/samtools-0.1.19/samtools mpileup -f ~/generic/data/hg19/ucsc.hg19.fasta -) \
		<(~/tools/samtools-0.1.19/samtools view -b -u -q 1 $(word 2,$^) | ~/tools/samtools-0.1.19/samtools mpileup -f ~/generic/data/hg19/ucsc.hg19.fasta -) \
		varscan/$*.varscan.part \
		--min-coverage 2 \
		--min-strands2 2 \
		--min-var-freq 0.2 \
		--normal-purity 1 \
		--tumor-purity 0.95 \
		--p-value 1 \
		--somatic-p-value 1 \
		--strand-filter 1 \
		--output-vcf 1 \
		2>&1 | $(LOG)
	mv varscan/$*.varscan.part.snp.vcf varscan/$*.varscan.snp.vcf
	mv varscan/$*.varscan.part.indel.vcf varscan/$*.varscan.indel.vcf

#-----------	
# SNPEFF
#-----------	
varscan/%.varscan.snp.dbsnp.snpeff.vcf: varscan/%.varscan.snp.dbsnp.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -Xmx2g -jar snpEff.jar -v -lof hg19 -stats $(PWD)/snpeff/$*.snpeff.summary.html $(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

varscan/%.varscan.snp.dbsnp.vcf: varscan/%.varscan.snp.vcf ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate -v ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf $(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	test -s $@.part
	mv $@.part $@

#-----------	
# FINAL LIST
#-----------	
filtered-variants/%.tsv: varscan/%.varscan.snp.dbsnp.snpeff.vcf ~/henning/scripts/filter-variants.pl
	perl ~/henning/scripts/filter-variants.pl \
		$< snp \
		--patient $* \
		--rmsk-file ~/generic/data/hg19/hg19.rmsk.txt.gz \
		--simpleRepeat-file ~/generic/data/hg19/hg19.simpleRepeat.txt.gz \
		--segdup-file ~/generic/data/hg19/hg19.genomicSuperDups.txt.gz \
		--blacklist-file ~/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz \
		--g1k-accessible ~/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz \
		--cosmic-mutation-file ~/generic/data/cosmic/v67/CosmicMutantExport_v67_241013.tsv \
		2>&1 1> $@.part | $(LOG)
	mv $@.part $@
