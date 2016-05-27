from subprocess import call, check_output
import os
import sys
import argparse
import pandas as pd
from multiprocessing import Pool
import inspect

PIPELINES_ROOT = "/users/pangwei/bds_atac"
CODE_ROOT = "/users/pangwei/stemcells/mesoderm"
DATA_ROOT = "/srv/scratch/pangwei/mesoderm/data/ATACseq"
OUTPUT_ROOT = "/srv/scratch/pangwei/mesoderm/bulk-atac/output-hg19-scidata"
GENOME_TYPE = "hg19"


def runFourTuple(inputTuple):
    """
    Helper function for ATACProcessor.mergePeaks. 
    This has to be declared at global scope for multiprocessing.Pool 
    to work correctly.
    """
    call("bash %s %s %s %s" % (inputTuple), shell=True)


def generateBDSScript():
    """
    Creates runATACpipeline.sh in the output directory.
    Combines replicates for celltypes, assuming that ATAC(N) and ATAC(N+1) are replicates.
    """
    scriptPath = os.path.join(OUTPUT_ROOT, 'runATACpipeline.sh')
    with open(scriptPath, 'w') as scriptFile:

        bdsPath = os.path.join(
            PIPELINES_ROOT,
            "atac.bds")

        for folder1 in dataFolders:
            
            if not folder1.startswith('ATAC'): continue
            ID = folder1.split('ATAC')[1]        
            assert(int(ID))
            ID = int(ID)
            if ID % 2 == 1:                    
                
                folder2 = 'ATAC%s' % (ID+1)             
                folder1_path = os.path.join(DATA_ROOT, folder1)
                folder2_path = os.path.join(DATA_ROOT, folder2)

                assert os.path.isdir(folder1_path)
                assert os.path.isdir(folder2_path)

                folder1_R1Path = os.path.join(folder1_path, "%s-trimmed-pair1.fastq.gz" % folder1)
                folder1_R2Path = os.path.join(folder1_path, "%s-trimmed-pair2.fastq.gz" % folder1)
                folder2_R1Path = os.path.join(folder2_path, "%s-trimmed-pair1.fastq.gz" % folder2)
                folder2_R2Path = os.path.join(folder2_path, "%s-trimmed-pair2.fastq.gz" % folder2)

                outputDir = os.path.join(OUTPUT_ROOT, "ATAC%s+%s" % (ID, ID+1))
                if not os.path.exists(outputDir):
                    os.makedirs(outputDir)

                scriptFile.write(("bds %s -out_dir %s -num_rep 2 "
                    "-fastq1_1 %s -fastq1_2 %s "
                    "-fastq2_1 %s -fastq2_2 %s "
                    "-subsample 17500000 -true_rep -species %s &") % 
                    (
                        bdsPath, outputDir, 
                        folder1_R1Path, folder1_R2Path, 
                        folder2_R1Path, folder2_R2Path, 
                        GENOME_TYPE))

                scriptFile.write("\n")
                
    print("Script written to %s" % scriptPath)    


def mergePeaksIDR():
    """
    Merge peaks across samples.

    This calls bedtools merge on all reproducible peaks (as measured by IDR)
    across all cell types.

    Outputs the following file in self.CODE_ROOT:
    1) mergedPeaks_idr.narrowPeak: 
        The genomic coordinates of each peak.
    """

    # Merge all peaks across samples
    catPeaksPath = os.path.join(OUTPUT_ROOT, 'catPeaks_idr.narrowPeak')
    mergedPeaksPath = os.path.join(OUTPUT_ROOT, 'mergedPeaks_idr.narrowPeak')
    
    mergeCommand = "zcat"
    print("Merging the following peaks:")
    for folder in outputFolders:
        print(folder)
        
        folderPath = os.path.join(OUTPUT_ROOT, folder)
        assert os.path.isdir(folderPath)

        narrowPeakPath = os.path.join(
            folderPath,
            "peak",
            "idr",
            "true_reps",
            "rep1-rep2",
            "rep1-rep2.IDR0.1.filt.narrowPeak.gz"
            )

        assert os.path.isfile(narrowPeakPath), "%s does not exist." % narrowPeakPath

        mergeCommand += " %s" % narrowPeakPath

    os.chdir(OUTPUT_ROOT)
    call("%s | sort -k1,1 -k2,2n > %s" % (mergeCommand, catPeaksPath), shell=True)
    call("bedtools merge -i %s > %s" % (catPeaksPath, mergedPeaksPath), shell=True)

    # Add peak names and dummy strand/score information to make it a proper 6-column BED file
    tempPath = mergedPeaksPath + ".temp"
    with open(mergedPeaksPath, "r") as r:
        with open(tempPath, "w") as w:
            for idx, line in enumerate(r):                        
                w.write("%s\tpeak%d\t0\t+\n" % (line.strip(), idx))

    call("cp %s %s" % (tempPath, mergedPeaksPath), shell=True)
    call("rm %s" % tempPath, shell=True)
    
    call("cp mergedPeaks_idr.narrowPeak %s" % (
        os.path.join(CODE_ROOT, 'mergedPeaks_idr.narrowPeak')), 
        shell=True)
    
    return


def getPValOnMergedPeaks():
    """
    Uses the universal list of reproducible peaks (merged across all cell types)
    to score the presence or absence (as measured by p-value) of each peak across 
    each cell type.

    For each cell type, we pool its two biological replicates together and call peaks (MACS2)
    on the pooled reads. (This is done automatically by the ATAC pipeline, not in this function.)

    To obtain a single measure of confidence at each peak P in the universal list for each cell-type 
    C, we took the highest -log10 p-value out of all peaks in the pooled replicates for C that
    intersected with P.

    Outputs the following files in self.CODE_ROOT:
    1) pvalOverlappingMergedPeaks_idr.txt: 
        A matrix containing the p-values as measured by the above method. 
        Every column is a cell type, and every row is a peak corresponding to rows of 
        mergedPeaks_idr.narrowPeak.
    2) pvalOverlappingMergedPeaks_idr.metadata: 
        The names of the cell types in the same order as the columns of the .txt file.
    """ 

    mergedPeaksPath = os.path.join(OUTPUT_ROOT, 'mergedPeaks_idr.narrowPeak')

    # Prepare parameters (filenames) for each sample
    parametersForMerge = []
    countsPaths = ""
    
    for folder in outputFolders:            
        folderPath = os.path.join(OUTPUT_ROOT, folder)
        tempPath = os.path.join(folderPath, "idr_%s.temp" % folder)
        assert os.path.isdir(folderPath)
    
        ID = int(folder.split('ATAC')[1].split('+')[0])            

        narrowPeakPath = os.path.join(
            folderPath,
            "peak",
            "macs2", 
            "pooled_rep",
            "ATAC%s-trimmed-pair1.trim.PE2SE.nodup.17.tn5_pooled.pf.narrowPeak.gz" % ID
            )
        assert os.path.isfile(narrowPeakPath), "%s does not have peaks called." % folder

        parametersForMerge.append(
            (os.path.join(CODE_ROOT, "getPValOnMergedPeaks.sh"), 
            narrowPeakPath,
            mergedPeaksPath,
            tempPath))

    # Count the number of fragments falling into each peak    
    p = Pool(len(outputFolders))
    p.map(runFourTuple, parametersForMerge)
    p.terminate()
        
    # For each sample, only take the highest IDR in each merged peak
    outputPaths = ""
    for folder in outputFolders:
        
        print(folder)
        folderPath = os.path.join(OUTPUT_ROOT, folder)
        tempPath = os.path.join(folderPath, "idr_%s.temp" % folder)
        outputPath = os.path.join(folderPath, "idr_%s.txt" % folder)

        df = pd.read_csv(tempPath, sep='\t', header=None)
        grouped = df.groupby(df[0]).max()

        # The groupby/max operation re-sorts the peaks lexicographically, so peak10 comes before peak2
        # We don't want this, so we re-sort by the integer bit, which is everything but the first 4 chars       
        for i in grouped.index:                
            grouped.loc[i, 'sortkey'] = i[4:]

        grouped['sortkey'] = grouped['sortkey'].astype('int')                
        grouped = grouped.sort_values(by='sortkey')

        grouped.to_csv(outputPath, header=False, index=False, columns=[1])
        call('rm %s' % tempPath, shell=True)

        outputPaths += " %s" % outputPath

    # Combine all of this data into one big matrix
    os.chdir(OUTPUT_ROOT)
    call("paste%s > pvalOverlappingMergedPeaks_idr.txt" % outputPaths, shell=True)
    
    with open(os.path.join(OUTPUT_ROOT, 'pvalOverlappingMergedPeaks_idr.metadata'), 'w') as f:
        f.write('\n'.join(outputFolders))
    
    call("cp pvalOverlappingMergedPeaks_idr.metadata %s" % (
        os.path.join(CODE_ROOT, 'pvalOverlappingMergedPeaks_idr.metadata')), 
        shell=True)

    # Copy it over to the shared code folder because RStudio is only available on mitra.stanford.edu
    call("cp pvalOverlappingMergedPeaks_idr.txt %s" % (
        os.path.join(CODE_ROOT, 'pvalOverlappingMergedPeaks_idr.txt')), 
        shell=True)        

if __name__ == '__main__':

    dataFolders = [folder for folder in os.listdir(DATA_ROOT) if 
                   os.path.isdir(os.path.join(DATA_ROOT, folder))]
    outputFolders = [folder for folder in os.listdir(OUTPUT_ROOT) if 
                     (os.path.isdir(os.path.join(OUTPUT_ROOT, folder)) and
                     folder.startswith('ATAC') and 'TEMP' not in folder)]

    if not os.path.exists(OUTPUT_ROOT):
        os.makedirs(OUTPUT_ROOT)

    args = sys.argv
    fxn_args = args[2:]
    print('Calling %s with arguments' % args[1], args[2:])
    locals()[args[1]](*args[2:])
