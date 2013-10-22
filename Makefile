export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work

SAMPLES = TEST

#bwa: $(foreach S, $(SAMPLES), bam/$S_s_5.bwa.bam)
#bwa: bam/TEST.bwa.merged.bam
all: picard/AML_022.multiplemetrics picard/AML_021.multiplemetrics

# created subsampled fastq for test purposes
~/henning/data/fastq/TEST_s_1_1_sequence.txt.gz: ~/henning/data/fastq/AML_022_s_1_1_sequence.txt.gz ~/henning/data/fastq/AML_022_s_1_2_sequence.txt.gz ~/henning/data/fastq/AML_022_s_2_1_sequence.txt.gz ~/henning/data/fastq/AML_022_s_2_2_sequence.txt.gz
	~/tools/seqtk/seqtk sample -s100 ~/henning/data/fastq/AML_022_s_1_1_sequence.txt.gz 10000 | gzip > ~/henning/data/fastq/TEST_s_1_1_sequence.txt.gz
	~/tools/seqtk/seqtk sample -s100 ~/henning/data/fastq/AML_022_s_1_2_sequence.txt.gz 10000 | gzip > ~/henning/data/fastq/TEST_s_1_2_sequence.txt.gz
	~/tools/seqtk/seqtk sample -s324 ~/henning/data/fastq/AML_022_s_2_1_sequence.txt.gz 10000 | gzip > ~/henning/data/fastq/TEST_s_2_1_sequence.txt.gz
	~/tools/seqtk/seqtk sample -s324 ~/henning/data/fastq/AML_022_s_2_2_sequence.txt.gz 10000 | gzip > ~/henning/data/fastq/TEST_s_2_2_sequence.txt.gz
	
bam/%.bwa.bam: ~/henning/data/fastq/%_1_sequence.txt.sai ~/henning/data/fastq/%_2_sequence.txt.sai ~/henning/data/fastq/%_1_sequence.txt.gz ~/henning/data/fastq/%_2_sequence.txt.gz
	~/tools/bwa-0.7.4/bwa sampe ~/generic/data/hg19/ucsc.hg19.fasta $^ \
		| ~/tools/samtools-0.1.19/samtools view -bS - \
		2>&1 1> $@.part | tee -a make.log
	mv $@.part $@
	~/tools/samtools-0.1.19/samtools sort $@ $@.sorted \
		2>&1 | tee -a make.log
	mv $@.sorted.bam $@
	~/tools/samtools-0.1.19/samtools index $@ \
		2>&1 | tee -a make.log
	ls $@.bai  # force make to quit if file does not exist

#.PRECIOUS: %.sai
%.sai: %.gz
	~/tools/bwa-0.7.4/bwa aln \
		-t 3 \
		~/generic/data/hg19/ucsc.hg19.fasta $^ \
		2>&1 1> $@.part | tee -a make.log
	mv $@.part $@
	
bam/%.bwa.merged.bam: bam/%_s_1.bwa.bam bam/%_s_2.bwa.bam bam/%_s_3.bwa.bam
	~/tools/samtools-0.1.19/samtools merge -f $@.part $^ \
		2>&1 | tee -a make.log
	rm $(foreach B, $^, $B.bai)
	mv $@.part $@

bam/%.dup.bam: bam/%.bam
	java -Xmx2g -jar ~/tools/picard-tools-1.101/MarkDuplicates.jar \
		INPUT=$< \
		OUTPUT=$@.part \
		METRICS_FILE=picard/$*.mark_duplicates_metrics \
		VALIDATION_STRINGENCY=LENIENT \
		2>&1 | tee -a make.log
	mv $@.part $@

picard/%.multiplemetrics: bam/%.bwa.merged.dup.bam
	java -Xmx2g -jar ~/tools/picard-tools-1.101/CollectMultipleMetrics.jar \
		INPUT=$< \
		OUTPUT=picard/$* \
		VALIDATION_STRINGENCY=LENIENT \
		PROGRAM=CollectAlignmentSummaryMetrics \
		PROGRAM=CollectInsertSizeMetrics \
		PROGRAM=QualityScoreDistribution \
		PROGRAM=MeanQualityByCycle \
		2>&1 | tee -a make.log
	touch $@

varscan/%.varscan.somatic: bam/%.bwa.merged.normal.bam bam/%.bwa.merged.tumor.bam
	java -jar ~/tools/varscan-2.3.6/VarScan.v2.3.6.jar somatic \
		<(~/tools/samtools-0.1.19/samtools view -b -u -q 1 $(word 1,$^) | ~/tools/samtools-0.1.19/samtools mpileup -f ~/generic/data/hg19/ucsc.hg19.fasta -) \
		<(~/tools/samtools-0.1.19/samtools view -b -u -q 1 $(word 2,$^) | ~/tools/samtools-0.1.19/samtools mpileup -f ~/generic/data/hg19/ucsc.hg19.fasta -) \
		$@ 2>&1 | tee -a make.log