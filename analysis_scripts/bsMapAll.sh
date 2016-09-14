#!/bin/bash

# Align the reads in parallel
cat fileManifest | parallel -j5 --colsep '\t' \
python /home/sbuckberry/working_data_01/bin/BSseeker2/bs_seeker2-align.py \
--aligner=bowtie2 --bt2--end-to-end --bt2-p 3 \
-i {1} \
-o {2} \
-p /usr/local/packages/bowtie2-2.2.5/ \
-d /scratchfs/sbuckberry/indexes/hs37d5_chrL.fa.gz_bowtie2/ \
-g /scratchfs/sbuckberry/indexes/hs37d5_chrL.fa.gz

# Call methylation in parallel
cat fileManifest | parallel -j5 --colsep '\t' \
python /home/sbuckberry/working_data_01/bin/BSseeker2/bs_seeker2-call_methylation.py \
-i {2} \
-o {3} \
-d /scratchfs/sbuckberry/indexes/hs37d5_chrL.fa.gz_bowtie2/

