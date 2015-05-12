# Compute coverage of regions masked and not masked by RepeatMasker.
# Based on a tutorial by Stephen Turner (http://gettinggeneticsdone.blogspot.com.au/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html)

# Compute the complement of all RepeatMaster tracks
# hg18
bedtools complement -i ../../repeat_masker/repeat_masker.hg18.bed.gz -g ../../repeat_masker/human.hg18.genome | bedtools merge -i - | gzip -c > ../../repeat_masker/repeat_masker.complement.hg18.bed.gz
# hg19
bedtools complement -i ../../repeat_masker/repeat_masker.hg19.bed.gz -g ../../repeat_masker/human.hg19.genome | bedtools merge -i - | gzip -c > ../../repeat_masker/repeat_masker.complement.hg19.bed.gz
# mm10
bedtools complement -i ../../repeat_masker/repeat_masker.mm10.bed.gz -g ../../repeat_masker/mouse.mm10.genome | bedtools merge -i - | gzip -c > ../../repeat_masker/repeat_masker.complement.mm10.bed.gz

# Compute coverage
# EPISCOPE data
cd EPISCOPE
parallel -j 12 "bedtools coverage -hist -abam {} -b ../../../repeat_masker/repeat_masker.hg19.bed.gz | grep ^all > {}.repeats.hist.all.txt" ::: *.bam
parallel -j 12 "bedtools coverage -hist -abam {} -b ../../../repeat_masker/repeat_masker.complement.hg19.bed.gz | grep ^all > {}.nonrepeats.hist.all.txt" ::: *.bam

# Lister data
cd ../Lister/
parallel -j 12 "bedtools coverage -hist -abam {} -b ../../../repeat_masker/repeat_masker.hg18.bed.gz | grep ^all > {}.repeats.hist.all.txt" ::: *.bam
parallel -j 12 "bedtools coverage -hist -abam {} -b ../../../repeat_masker/repeat_masker.complement.hg18.bed.gz | grep ^all > {}.nonrepeats.hist.all.txt" ::: *.bam

# Seisenberger data
cd ../Seisenberger/
parallel -j 3 "bedtools coverage -hist -abam {} -b ../../../repeat_masker/repeat_masker.mm10.bed.gz | grep ^all > {}.repeats.hist.all.txt" ::: *.bam
parallel -j 3 "bedtools coverage -hist -abam {} -b ../../../repeat_masker/repeat_masker.complement.mm10.bed.gz | grep ^all > {}.nonrepeats.hist.all.txt" ::: *.bam

# Ziller data
cd ../Ziller/
parallel -j 12 "bedtools coverage -hist -abam {} -b ../../../repeat_masker/repeat_masker.hg19.bed.gz | grep ^all > {}.repeats.hist.all.txt" ::: *.bam
parallel -j 12 "bedtools coverage -hist -abam {} -b ../../../repeat_masker/repeat_masker.complement.hg19.bed.gz | grep ^all > {}.nonrepeats.hist.all.txt" ::: *.bam
