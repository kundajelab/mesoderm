# Human mesoderm differentiation
This repository contains scripts for dataset processing and QC for the mesoderm differentiation project.
If you have questions, please contact Pang Wei Koh (<pangwei@cs.stanford.edu>).

## RNA-seq
bulkDataViz.Rmd and scDataViz.Rmd (together with their corresponding .html files) reproduce Figures 2
in the Scientific Data manuscript. scAverageCorrelation.r measures the correlation
of gene expression between the bulk population and the average of single cells.

To process the raw RNA-seq reads, first run STAR\_RSEM\_prep.py and then run\_STAR\_RSEM.py. 
This will call STAR (for alignment) and then RSEM (for quantification) on each sample.
Next, run getExpr.py to filter out samples with poor mapping statistics and construct 
a matrix of gene expression across samples from the remaining samples.

## ATAC-seq
The ATAC-seq folder contains one file, mesodermATAC.py.

After downloading the data and changing the paths defined at the top of the file appropriately, 
run `python mesodermATAC.py generateBDSScript` to generate a .sh file that will call the ATAqC pipeline
on the samples. Figure 3 in the Scientific Data manuscript was taken from the output of the pipeline.

Run `python mesodermATAC.py mergePeaksIDR` followed by `python mesodermATAC.py getPValOnMergedPeaks` after the 
pipeline has finished running to obtain a universal list of reproducible peaks across all cell types.

## Surface marker screening
surfaceMarkerViz.Rmd (and its corresponding .html file) reproduces Figure 4 in the Scientific Data manuscript.
