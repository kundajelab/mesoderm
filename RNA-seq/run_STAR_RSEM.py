# This script runs STAR_RSEM.sh separately on each sample.

# TO USE: Replace paths in dataRoot, outputRoot, codeRoot, RSEMrefDir, and STARgenomeDir.
#         The latter two paths should be 

import os
from subprocess import call

# Set bulk_or_sc to either 'bulk' or 'sc'
# bulk_or_sc = 'bulk'
bulk_or_sc = 'sc'

if bulk_or_sc == 'bulk':
    dataRoot = "/srv/scratch/pangwei/mesoderm/data/bulk_RNAseq"
    outputRoot = "/srv/scratch/pangwei/mesoderm/bulk-rna/output"
elif bulk_or_sc == 'sc':
    dataRoot = "/srv/scratch/pangwei/mesoderm/data/single_cell"
    outputRoot = "/srv/scratch/pangwei/mesoderm/sc-rna/output"
    
codeRoot = "/users/pangwei/stemcells/mesoderm"
scriptName = "STAR_RSEM.sh"

RSEMrefDir = "/srv/scratch/pangwei/mesoderm/bulk-rna/genome/RSEMref"
STARgenomeDir = "/srv/scratch/pangwei/mesoderm/bulk-rna/genome"

dataType = "str_PE"
nThreadsSTAR = "8" 
nThreadsRSEM = "8" 

foldersToInclude = [folder for folder in os.listdir(dataRoot) if os.path.isdir(os.path.join(dataRoot, folder))]

for folderIdx, folder in enumerate(foldersToInclude):

    print('Processing %s...' % folder)

    folderPath = os.path.join(dataRoot, folder)
    assert os.path.isdir(folderPath), "%s does not exist" % folderPath

    outputPath = os.path.join(outputRoot, folder)

    # Check if Quant.genes.results already exists in outputPath. If so, skip.
    if os.path.isfile(os.path.join(outputPath, 'Quant.genes.results')):
        print('    %s has already been processed. Skipping.' % folder)
        continue
    
    os.chdir(folderPath)
    call("cp %s ." % os.path.join(codeRoot, scriptName), shell=True)
    call("chmod u+x %s" % scriptName, shell=True)

    # This is for the bulk-RNA file format
    if bulk_or_sc == 'bulk':
        call("./%s %s_R1.fastq_trimmed.gz %s_R2.fastq_trimmed.gz %s %s %s %s %s" % \
             (scriptName, folder, folder, STARgenomeDir, 
              RSEMrefDir, dataType, nThreadsSTAR, nThreadsRSEM), shell=True)
    # This is for the sc-RNA file format
    else:    
        call("./%s %s-trimmed-pair1.fastq.gz %s-trimmed-pair2.fastq.gz %s %s %s %s %s" % \
             (scriptName, folder, folder, STARgenomeDir, 
              RSEMrefDir, dataType, nThreadsSTAR, nThreadsRSEM), shell=True)

    # Output gets put in the working directory by default:
    # Aligned.sortedByCoord.out.bam                 # alignments, standard sorted BAM, agreed upon formatting
    # Log.final.out                                 # mapping statistics to be used for QC, text, STAR formatting
    # Quant.genes.results                           # RSEM gene quantifications, tab separated text, RSEM formatting
    # Quant.isoforms.results                        # RSEM transcript quantifications, tab separated text, RSEM formatting
    # Quant.pdf                                     # RSEM diagnostic plots
    # Signal.{Unique,UniqueMultiple}.strand{+,-}.bw # 4 bigWig files for stranded data
    # Signal.{Unique,UniqueMultiple}.unstranded.bw  # 2 bigWig files for unstranded data

    # Move output files to desired output ocation
    call("mkdir -p %s" % outputPath, shell=True)
    call("mv *.bam chrNL.txt sig.tmp SJ.out.tab Log* Signal* Quant* %s" % outputPath, shell=True)