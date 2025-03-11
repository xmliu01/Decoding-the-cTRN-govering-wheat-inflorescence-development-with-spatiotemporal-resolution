# Decoding-the-cTRN-gevering-wheat-inflorescence-development-with-spatiotemporal-resolution
These are the codes used in the article entitled "Docoding the cTRN govering wheat inflorescence development with spatiotemporal resolution"

To run these codes, you need to install R packages Seurat(v4.3.0), ArchR (v1.0.3), tidyverse, Monocle3, sctransform, BSgenome

snRNA.r: code used for snRNA-seq data analysis.
snATAC.r: code used for snATAC-seq data analysis.
Tangram.py: code used for scStereo-seq data analysis.
Monocle3.r: code used for trajectory inference.

The 'Test' folder contains test data, test code, and test results.
For snRNA-seq:
'snRNA.TestData' folder contains input data for snRNA-seq analysis test.
snRNA.test.r: test code for snRNA-seq data anslysis.
'test.snRNA.celltype.specific.genes.txt' and 'test.snRNA.Seurat.rds' are the results for snRNA-seq data analysis test.

For snATAC-seq:
'test.snATAC.fragment.tsv.gz' and 'test.snATAC.fragment.tsv.gz.tbi' are the input data for snATAC-seq analysis test.
snATAC.test.r: test code for snATAC-seq data analysis.
'test.snATAC.output' folder contains the results for snATAC-seq data analysis test.

