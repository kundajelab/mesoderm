# This file prepares genomes for STAR v2.4 and RSEM v1.2.21.
# TO USE: Replace paths.

# Modified from https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/DAC/STAR_RSEM_prep.sh
# Commit 313830c7c10e8567091131c40bdec2b9477627e0

import os
from subprocess import call

# STAR genome directory
STARgenomeDir = "/srv/scratch/pangwei/mesoderm/bulk-rna/genome/"

# RSEM genome directory
RSEMgenomeDir = "/srv/scratch/pangwei/mesoderm/bulk-rna/genome/"

# fasta file(s)
fastaGenome = "/srv/scratch/pangwei/mesoderm/bulk-rna/genome/hg20.GRCh38/GrCh38_ERCC92.genome.fa"

# all-inclusive gtf file
gtf = "/srv/scratch/pangwei/mesoderm/bulk-rna/genome/hg20.GRCh38/GENCODE_ann/gencode.v22/gencode.v22_ERCC92.annotation.gtf"

assert os.path.isdir(STARgenomeDir) and os.path.isdir(RSEMgenomeDir)

# RSEM-prepare-reference 
RSEMcommand = "rsem-prepare-reference --gtf %s %s %s" % \
              (gtf, fastaGenome, os.path.join(RSEMgenomeDir, 'RSEMref'))
call(RSEMcommand, shell=True)

# STAR genome
STARcommand = ("/users/pangwei/STAR/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir %s "
               "--genomeFastaFiles %s --sjdbGTFfile %s "
               "--sjdbOverhang 100 --outFileNamePrefix %s --runThreadN 24") % \
               (STARgenomeDir, fastaGenome, gtf, STARgenomeDir)
call(STARcommand, shell=True)