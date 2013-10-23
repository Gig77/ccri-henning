export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work

MERGE = ~/tools/samtools-0.1.19/samtools merge -f $@.part $^ 2>&1 | tee -a make.log; mv $@.part $@;

all: TEST AML_021 AML_022 AML_039 AML_040
INTERMEDIATE_BAM = $(TEST_BAM) $(AML_021_BAM) $(AML_022_BAM) $(AML_039_BAM) $(AML_040_BAM)

TEST: bam/TEST.bwa.merged.dup.bam bam/TEST.bwa.merged.dup.bam.bai picard/TEST.multiplemetrics picard/TEST.hs_metrics
TEST_BAM = bam/TEST_s_1.bwa.bam bam/TEST_s_2.bwa.bam
bam/TEST.bwa.merged.bam: $(TEST_BAM)
	$(MERGE)

AML_021: bam/AML_021.bwa.merged.dup.bam bam/AML_021.bwa.merged.dup.bam.bai bam/AML_021.bwa.merged.bam.bai picard/AML_021.multiplemetrics picard/AML_021.hs_metrics
AML_021_BAM = bam/AML_021_s_1.bwa.bam bam/AML_021_s_2.bwa.bam bam/AML_021_s_3.bwa.bam bam/AML_021_s_4.bwa.bam bam/AML_021_s_5.bwa.bam bam/AML_021_s_6.bwa.bam bam/AML_021_s_7.bwa.bam bam/AML_021_s_8.bwa.bam
bam/AML_021.bwa.merged.bam: $(AML_021_BAM)
	$(MERGE)

AML_022: bam/AML_022.bwa.merged.dup.bam bam/AML_022.bwa.merged.dup.bam.bai bam/AML_022.bwa.merged.bam.bai picard/AML_022.multiplemetrics picard/AML_022.hs_metrics
AML_022_BAM = bam/AML_022_s_1.bwa.bam bam/AML_022_s_2.bwa.bam bam/AML_022_s_3.bwa.bam bam/AML_022_s_4.bwa.bam bam/AML_022_s_5.bwa.bam bam/AML_022_s_6.bwa.bam bam/AML_022_s_7.bwa.bam bam/AML_022_s_8.bwa.bam
bam/AML_022.bwa.merged.bam: $(AML_022_BAM)
	$(MERGE)

AML_039: bam/AML_039.bwa.merged.dup.bam bam/AML_039.bwa.merged.dup.bam.bai bam/AML_039.bwa.merged.bam.bai picard/AML_039.multiplemetrics picard/AML_039.hs_metrics
AML_039_BAM = bam/AML_039_s_1.bwa.bam bam/AML_039_s_2.bwa.bam bam/AML_039_s_3.bwa.bam bam/AML_039_s_4.bwa.bam bam/AML_039_s_5.bwa.bam bam/AML_039_s_6.bwa.bam bam/AML_039_s_7.bwa.bam bam/AML_039_s_8.bwa.bam
bam/AML_039.bwa.merged.bam: $(AML_039_BAM)
	$(MERGE)

AML_040: bam/AML_040.bwa.merged.dup.bam bam/AML_040.bwa.merged.dup.bam.bai bam/AML_040.bwa.merged.bam.bai picard/AML_040.multiplemetrics picard/AML_040.hs_metrics
AML_040_BAM = bam/AML_040_s_1.bwa.bam bam/AML_040_s_2.bwa.bam bam/AML_040_s_3.bwa.bam bam/AML_040_s_4.bwa.bam bam/AML_040_s_5.bwa.bam bam/AML_040_s_6.bwa.bam bam/AML_040_s_7.bwa.bam bam/AML_040_s_8.bwa.bam
bam/AML_040.bwa.merged.bam: $(AML_040_BAM)
	$(MERGE)

# created subsampled fastq for test purposes
~/henning/data/fastq/TEST_s_1_1_sequence.txt.gz: ~/henning/data/fastq/AML_022_s_1_1_sequence.txt.gz ~/henning/data/fastq/AML_022_s_1_2_sequence.txt.gz ~/henning/data/fastq/AML_022_s_2_1_sequence.txt.gz ~/henning/data/fastq/AML_022_s_2_2_sequence.txt.gz
	~/tools/seqtk/seqtk sample -s100 ~/henning/data/fastq/AML_022_s_1_1_sequence.txt.gz 10000 | gzip > ~/henning/data/fastq/TEST_s_1_1_sequence.txt.gz
	~/tools/seqtk/seqtk sample -s100 ~/henning/data/fastq/AML_022_s_1_2_sequence.txt.gz 10000 | gzip > ~/henning/data/fastq/TEST_s_1_2_sequence.txt.gz
	~/tools/seqtk/seqtk sample -s324 ~/henning/data/fastq/AML_022_s_2_1_sequence.txt.gz 10000 | gzip > ~/henning/data/fastq/TEST_s_2_1_sequence.txt.gz
	~/tools/seqtk/seqtk sample -s324 ~/henning/data/fastq/AML_022_s_2_2_sequence.txt.gz 10000 | gzip > ~/henning/data/fastq/TEST_s_2_2_sequence.txt.gz
	
.INTERMEDIATE: $(INTERMEDIATE_BAM) 
bam/%.bwa.bam: ~/henning/data/fastq/%_1_sequence.txt.sai ~/henning/data/fastq/%_2_sequence.txt.sai ~/henning/data/fastq/%_1_sequence.txt.gz ~/henning/data/fastq/%_2_sequence.txt.gz
	echo `date` [START] TARGET=$@ PRE=$^ | tee -a make.log
	~/tools/bwa-0.7.4/bwa sampe ~/generic/data/hg19/ucsc.hg19.fasta $^ \
		| ~/tools/samtools-0.1.19/samtools view -bS - \
		2>&1 1> $@.part | tee -a make.log
	mv $@.part $@
	~/tools/samtools-0.1.19/samtools sort $@ $@.sorted \
		2>&1 | tee -a make.log
	mv $@.sorted.bam $@
	echo `date` [END] TARGET=$@ PRE=$^ | tee -a make.log

%.sai: %.gz
	echo `date` [START] TARGET=$@ PRE=$^ | tee -a make.log
	~/tools/bwa-0.7.4/bwa aln \
		-t 3 \
		~/generic/data/hg19/ucsc.hg19.fasta $^ \
		2>&1 1> $@.part | tee -a make.log
	mv $@.part $@
	echo `date` [END] TARGET=$@ PRE=$^ | tee -a make.log

bam/%.bwa.merged.dup.bam: bam/%.bwa.merged.bam
	echo `date` [START] TARGET=$@ PRE=$^ | tee -a make.log
	java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar ~/tools/picard-tools-1.101/MarkDuplicates.jar \
		INPUT=$< \
		OUTPUT=$@.part \
		METRICS_FILE=picard/$*.mark_duplicates_metrics \
		VALIDATION_STRINGENCY=LENIENT \
		2>&1 | tee -a make.log
	mv $@.part $@
	echo `date` [END] TARGET=$@ PRE=$^ | tee -a make.log

bam/%.bam.bai: bam/%.bam
	echo `date` [START] TARGET=$@ PRE=$^ | tee -a make.log
	rm -f $@
	~/tools/samtools-0.1.19/samtools index $^
	echo `date` [END] TARGET=$@ PRE=$^ | tee -a make.log

picard/%.multiplemetrics: bam/%.bwa.merged.dup.bam
	echo `date` [START] TARGET=$@ PRE=$^ | tee -a make.log
	java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar ~/tools/picard-tools-1.101/CollectMultipleMetrics.jar \
		INPUT=$< \
		OUTPUT=picard/$* \
		VALIDATION_STRINGENCY=LENIENT \
		PROGRAM=CollectAlignmentSummaryMetrics \
		PROGRAM=CollectInsertSizeMetrics \
		PROGRAM=QualityScoreDistribution \
		PROGRAM=MeanQualityByCycle \
		2>&1 | tee -a make.log
	touch $@
	echo `date` [END] TARGET=$@ PRE=$^ | tee -a make.log

picard/%.hs_metrics: bam/%.bwa.merged.dup.bam
	echo `date` [START] TARGET=$@ PRE=$^ | tee -a make.log
	# truseq_exome_targeted_regions.hg19.bed.chr downloaded from http://supportres.illumina.com/documents/myillumina/5dfd7e70-c4a5-405a-8131-33f683414fb7/truseq_exome_targeted_regions.hg19.bed.chr.gz
	~/tools/samtools-0.1.19/samtools view -H bam/$*.bwa.merged.dup.bam > picard/$*.truseq-for-picard.bed
	gawk 'BEGIN { OFS="\t"} {print $$1,$$2,$$3,$$6,$$4 }' ~/generic/data/illumina/truseq_exome_targeted_regions.hg19.bed.chr >> picard/$*.truseq-for-picard.bed
	java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar ~/tools/picard-tools-1.101/CalculateHsMetrics.jar \
		BAIT_INTERVALS=picard/$*.truseq-for-picard.bed \
		TARGET_INTERVALS=picard/$*.truseq-for-picard.bed \
		INPUT=bam/$*.bwa.merged.dup.bam \
		OUTPUT=$@.part \
		REFERENCE_SEQUENCE=~/generic/data/hg19/ucsc.hg19.fasta \
		PER_TARGET_COVERAGE=picard/$*.hs_metrics.per_target_coverage.part \
		VALIDATION_STRINGENCY=LENIENT \
		2>&1 | tee -a make.log
	mv $@.part $@
	mv picard/$*.hs_metrics.per_target_coverage.part picard/$*.hs_metrics.per_target_coverage
	rm picard/$*.truseq-for-picard.bed
	echo `date` [END] TARGET=$@ PRE=$^ | tee -a make.log

varscan/MHH_121.varscan: bam/AML_022.bwa.merged.dup.bam bam/AML_021.bwa.merged.dup.bam
	echo `date` [START] TARGET=$@ PRE=$^ | tee -a make.log
	java -jar ~/tools/varscan-2.3.6/VarScan.v2.3.6.jar somatic \
		<(~/tools/samtools-0.1.19/samtools view -b -u -q 1 $(word 1,$^) | ~/tools/samtools-0.1.19/samtools mpileup -f ~/generic/data/hg19/ucsc.hg19.fasta -) \
		<(~/tools/samtools-0.1.19/samtools view -b -u -q 1 $(word 2,$^) | ~/tools/samtools-0.1.19/samtools mpileup -f ~/generic/data/hg19/ucsc.hg19.fasta -) \
		$@ \
		--min-coverage 2 \
		--min-strands2 2 \
		--min-var-freq 0.2 \
		--normal-purity 1 \
		--tumor-purity 0.95 \
		--p-value 1 \
		--somatic-p-value 1 \
		--strand-filter 1 \
		--output-vcf 1 \
		2>&1 | tee -a make.log
	touch $@
	echo `date` [END] TARGET=$@ PRE=$^ | tee -a make.log
	
varscan/%.varscan.snp.snpeff.vcf: varscan/%.varscan.snp.vcf
	echo `date` [START] TARGET=$@ PRE=$^ | tee -a make.log
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -Xmx2g -jar ~/tools/snpEff-3.3h/snpEff.jar -v -lof hg19 $(PWD)/$< > $(PWD)/$@.part)
	mv $@.part $@
	mv ~/tools/snpEff-3.3h/snpEff_summary.html snpeff/$*.snpEff_summary.html
	mv ~/tools/snpEff-3.3h/snpEff_genes.txt snpeff/$*.snpEff_genes.txt
	echo `date` [END] TARGET=$@ PRE=$^ | tee -a make.log

varscan/%.varscan.snp.snpeff.dbsnp.vcf: varscan/%.varscan.snp.snpeff.vcf ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf
	echo `date` [START] TARGET=$@ PRE=$^ | tee -a make.log
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate -v ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf $(PWD)/$< > $(PWD)/$@.part) 
	mv $@.part $@
	echo `date` [END] TARGET=$@ PRE=$^ | tee -a make.log

filtered-variants.tsv: varscan/MHH_121.varscan.snp.snpeff.dbsnp.vcf ~/henning/scripts/filter-variants.pl
	perl ~/henning/scripts/filter-variants.pl \
		varscan/MHH_121.varscan.snp.snpeff.dbsnp.vcf snp \
		--header 1 \
		--rmsk-file ~/generic/data/hg19/hg19.rmsk.txt.gz \
		--simpleRepeat-file ~/generic/data/hg19/hg19.simpleRepeat.txt.gz \
		--segdup-file ~/generic/data/hg19/hg19.genomicSuperDups.txt.gz \
		--blacklist-file ~/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz \
		--g1k-accessible ~/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz \
		--cosmic-mutation-file ~/generic/data/cosmic/v65/CosmicMutantExport_v65_280513.tsv \
		> $@.part
	mv $@.part $@
