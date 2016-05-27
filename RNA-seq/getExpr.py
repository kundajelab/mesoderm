"""
Creates three files, {bulk/sc}_tpm.txt, {bulk/sc}_expr.txt and {bulk/sc}_geneLengths.txt, 
that contain a column for each sample containing gene expression (measured in TPM), gene expression
(measured in counts), and gene lengths, respectively. 
Each row is a gene. There is a header row (sample names) and header column (gene names).

minUniquePercentage is a percentage between 0-100. We use 70 for sc and 50 for bulk.
"""

from subprocess import call, check_output
import os
import argparse
import pandas as pd
import IPython

# Choose one
# bulk_or_sc = 'bulk'
bulk_or_sc = 'sc'

# bulk-RNA
if bulk_or_sc == 'bulk':
    outputRoot = "/srv/scratch/pangwei/mesoderm/bulk-rna/output"
    minUniquePercentage = 50
elif bulk_or_sc == 'sc':
    outputRoot = "/srv/scratch/pangwei/mesoderm/sc-rna/old-output"
    minUniquePercentage = 70

# Common to both
pipelinesRoot = "/users/pangwei/pipelines"
outputFolder = "/users/pangwei/stemcells/mesoderm"
minUniqueReads = 1000000
foldersToInclude = [folder for folder in os.listdir(outputRoot) if os.path.isdir(os.path.join(outputRoot, folder))]


# First, filter out the folders in foldersToInclude that haven't finished processing
# or that don't meet minUniqueReads and minUniquePercentage
filteredFolders = []

for folder in foldersToInclude:
    logPath = os.path.join(outputRoot, folder, 'Log.final.out')

    if not os.path.isfile(logPath): 
        continue

    uniqueReads = int(
        check_output(
            'cat %s | grep "Uniquely mapped reads number" | cut -f2' % logPath, 
            shell=True))

    print('Uniquely mapped reads for %s: %s' % (folder, uniqueReads))

    if uniqueReads < minUniqueReads:
        print('%s does not have enough uniquely mapped reads (%s).' % (folder, uniqueReads))
        continue

    uniquePercentage = float(
        check_output(
            'cat %s | grep "Uniquely mapped reads %%" | cut -f2' % logPath, 
            shell=True).split('%')[0])

    print('Uniquely mapped read %% for %s: %s' % (folder, uniquePercentage))

    if uniquePercentage < minUniquePercentage:
        print('%s does not have a high enough %% of uniquely mapped reads (%s).' % (folder, uniquePercentage))
        continue        

    filteredFolders.append(folder)

os.chdir(outputRoot)

# Get gene names
pasteCommand = "paste <(cut -f1 %s | tail -n +2)" % os.path.join(
    outputRoot,
    foldersToInclude[0],
    "Quant.genes.results")

tpmCommand = ""
countCommand = ""
lengthCommand = ""

# For each sample, append a column of the counts 
for folder in filteredFolders:

    resultsPath = os.path.join(
        outputRoot, 
        folder, 
        "Quant.genes.results")

    assert os.path.isfile(resultsPath), "No expression data for %s" % folder
    
    tpmCommand += " <(cut -f6 %s | tail -n +2)" % resultsPath            
    countCommand += " <(cut -f5 %s | tail -n +2)" % resultsPath
    lengthCommand += " <(cut -f4 %s | tail -n +2)" % resultsPath
    print("Pasting results from %s" % resultsPath)
    
# Concatenate folder names and add this as a header
call("cat <(echo -e '%s') <(%s) > %s" % (
        '\t'.join(filteredFolders), 
        pasteCommand + tpmCommand,
        os.path.join(outputFolder, "%s_tpm.txt" % bulk_or_sc)), 
    shell=True,
    executable="/bin/bash")

# call("cat <(echo -e '%s') <(%s) > %s" % (
        '\t'.join(filteredFolders), 
        pasteCommand + countCommand,
        os.path.join(outputFolder, "%s_expr.txt" % bulk_or_sc)), 
    shell=True,
    executable="/bin/bash")

# call("cat <(echo -e '%s') <(%s) > %s" % (
        '\t'.join(filteredFolders), 
        pasteCommand + lengthCommand,
        os.path.join(outputFolder, "%s_geneLengths.txt" % bulk_or_sc)), 
    shell=True,
    executable="/bin/bash")
