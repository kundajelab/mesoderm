import os
from subprocess import call

# bulk_or_sc = 'bulk'
bulk_or_sc = 'sc'


# bulk-RNA
if bulk_or_sc == 'bulk':
    dataRoot = "/srv/scratch/pangwei/mesoderm/data/bulk_RNAseq"
    outputRoot = "/srv/scratch/pangwei/mesoderm/bulk-rna/output"
# sc-RNA
elif bulk_or_sc == 'sc':
    # dataRoot = "/mnt/lab_data/kundaje/users/pangwei/mesoderm/single_cell"
    # dataRoot = "/srv/scratch/pangwei/mesoderm/data/single_cell"
    # outputRoot = "/srv/scratch/pangwei/mesoderm/sc-rna/output"
    dataRoot = "/srv/scratch/pangwei/mesoderm/data/single_cell_unfiltered"
    outputRoot = "/srv/scratch/pangwei/mesoderm/sc-rna/output_unfiltered"

# common to both
codeRoot = "/users/pangwei/stemcells/mesoderm"
scriptName = "STAR_RSEM.sh"

# This uses the Quake lab hg19 files with ERCC
# STARgenomeDir = "/srv/scratch/pangwei/mesoderm/data/genomes"
# RSEMrefDir = "/srv/scratch/pangwei/mesoderm/data/genomes/RSEMref"

# This uses hg38 files with ERCC
RSEMrefDir = "/srv/scratch/pangwei/mesoderm/bulk-rna/genome/RSEMref"
STARgenomeDir = "/srv/scratch/pangwei/mesoderm/bulk-rna/genome"

dataType = "str_PE"
nThreadsSTAR = "8" 
nThreadsRSEM = "8" 

# bulk-RNA
# foldersToInclude = [
# 'GD01',
# 'GD02',
# 'GD03',
# 'GD04',
# 'GD05',
# 'GD07',
# 'GD08',
# 'GD12',
# 'GD13',
# 'GD14',
# 'GD16',
# 'GD17',
# 'GD19',
# 'GD20',
# 'GD24',
# 'GD25',
# 'GD26',
# 'GD27',
# 'RM1',
# 'RM2',
# 'RM3',
# 'RM4',
# 'RM5',
# 'RM6',
# 'RM8',
# 'RM9'
# ]

# sc-RNA
foldersToInclude = [folder for folder in os.listdir(dataRoot) if os.path.isdir(os.path.join(dataRoot, folder))]

# Set this to "parallelize" things
numJobs = 4
jobID = 3

for folderIdx, folder in enumerate(foldersToInclude):

    if (folderIdx % numJobs) != jobID:
        continue

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

    # ./STAR_RSEM.sh (read1) (read2 or "") (STARgenomeDir) (RSEMrefDir) (dataType) (nThreadsSTAR) (nThreadsRSEM) 
    # This is for the bulk-RNA file format
    # call("./%s %s_R1.fastq_trimmed.gz %s_R2.fastq_trimmed.gz %s %s %s %s %s" % \
    #      (scriptName, folder, folder, STARgenomeDir, 
    #       RSEMrefDir, dataType, nThreadsSTAR, nThreadsRSEM), shell=True)

    # This is for the sc-RNA file format
    call("./%s %s-trimmed-pair1.fastq.gz %s-trimmed-pair2.fastq.gz %s %s %s %s %s" % \
         (scriptName, folder, folder, STARgenomeDir, 
          RSEMrefDir, dataType, nThreadsSTAR, nThreadsRSEM), shell=True)

    # output: all in the working directory, fixed names
    # Aligned.sortedByCoord.out.bam                 # alignments, standard sorted BAM, agreed upon formatting
    # Log.final.out                                 # mapping statistics to be used for QC, text, STAR formatting
    # Quant.genes.results                           # RSEM gene quantifications, tab separated text, RSEM formatting
    # Quant.isoforms.results                        # RSEM transcript quantifications, tab separated text, RSEM formatting
    # Quant.pdf                                     # RSEM diagnostic plots
    # Signal.{Unique,UniqueMultiple}.strand{+,-}.bw # 4 bigWig files for stranded data
    # Signal.{Unique,UniqueMultiple}.unstranded.bw  # 2 bigWig files for unstranded data

    call("mkdir -p %s" % outputPath, shell=True)
    
    #print("mv Aligned.sortedByCoord.out.bam Log.final.out Quant.genes.results Quant.isoforms.results Quan.pdf Signal.*.bw %s" % outputPath)
    call("mv *.bam chrNL.txt sig.tmp SJ.out.tab Log* Signal* Quant* %s" % outputPath, shell=True)









