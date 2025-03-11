# scRNAseq-of-sCOTR-iNK
single-cell RNA-seq of synthetic COTR-iNK cells derived from hPSCs
Title: Synthetic COTR complex empowers iNK cells with intracellular antigen-targeting capability

Summary:  
Raw data of droplet-based scRNA-seq （10 × Genomics） can be downloaded and processed by CellRanger software package (version 7.0.1), then subjected to Seurat (version 4.3.0) for further analysis.  

Raw data of scRNA-seq:  
The scRNA-seq data of UCB-NK (NK), ESC-iNK (iNK), EBV sCOTR-iNK have been deposited in the GSA public database. The accession number of ESC-iNK is HRA001609 (published in "Lateral plate mesoderm cell-based organoid system for NK cell regeneration from human pluripotent stem cells "([https://pubmed.ncbi.nlm.nih.gov/36344493/](https://pubmed.ncbi.nlm.nih.gov/36344493/)). The accession number of UCB-NK is HRA007978 ([https://www.biorxiv.org/content/10.1101/2024.07.30.605741v1](https://www.biorxiv.org/content/10.1101/2024.07.30.605741v1)). The accession number of EBV sCOTR-iNK is HRA010093 (unpublished in this work).

# scRNA seq analysis pipeline
## Merge the UCB-NK and iNK scRNA-seq data, and present the expression levels of CD3, CD8, and CD28 subunits
    01_geneExpression_UCBNK_and_iNK.R

## Identification of DEGS and GO terms enriched between (EBV sCOTR-iNK vesus iNK)
    02_DEanalysis_GOenrichment.R

## Gene set-enrichment analysis (GSEA)
    03_GSEA.R
