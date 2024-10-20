REF := ref/contigs.fa
READS := $(wildcard reads/*_R1.fastq.gz)
SAMPLES := $(patsubst reads/%_R1.fastq.gz,%,$(READS))
BAMS := $(patsubst %,bam/%.bam,$(SAMPLES))
PROFILE_DIRS := $(patsubst %,bam/%,$(SAMPLES))

.PHONY: all

all: PROFILE.db

# Index reference if not already indexed
$(REF).bwt: $(REF)
	bwa index $(REF)

# Create BAM files
bam/%.bam: reads/%_R1.fastq.gz reads/%_R2.fastq.gz $(REF).bwt
	mkdir -p bam
	bwa mem -t 4 $(REF) $< $(word 2,$^) | samtools view -bS -F4 - | samtools sort -o $@ --write-index -
	touch $@.done

# Create CONTIGS.db
CONTIGS.db: $(REF)
	anvi-gen-contigs-database -f $(REF) -o CONTIGS.db
	anvi-run-hmms -c CONTIGS.db -T 4 
	anvi-run-ncbi-cogs -c CONTIGS.db -T 4
	anvi-scan-trnas -c CONTIGS.db 
	anvi-run-scg-taxonomy -c CONTIGS.db
	touch contigs.done

# Create profile directories
bam/%: bam/%.bam CONTIGS.db
	anvi-profile -i $< -c CONTIGS.db -T 4 -o $@
	touch $@.profile.done

# Merge profiles
PROFILE.db: $(PROFILE_DIRS) CONTIGS.db
	echo "anvi-merge $(addsuffix /PROFILE.db,$(PROFILE_DIRS)) -o ./MINITEST -c CONTIGS.db" > PROFILE.done
	anvi-merge $(addsuffix /PROFILE.db,$(PROFILE_DIRS)) -o ./MINITEST -c CONTIGS.db
	mv MINITEST/*.db .
	rm -rf MINITEST

# Clean rule
.PHONY: clean
clean:
	rm -f $(REF).bwt $(REF).pac $(REF).ann $(REF).amb $(REF).sa
	rm -rf bam
	rm -f CONTIGS.db
	rm -f PROFILE.db
	rm -f AUXILIARY-DATA.db
	rm *.tsv *.done 2>/dev/null || true
