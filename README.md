# single_cell_Lineage_Tracing_iFlpscLineage
This repository contains a snippet of the analysis I performed to validate the scLT technology iFlpscLineage that can be used for clonal analysis coupled with transcriptomical data. This work was carried during my PhD thesis at the Molecular Genetics of Angiogenesis at CNIC.

## single cell Lineage Tracing (scLT)

single cell Lineage Tracing seeks to follow the division history of a cell in order to establish the relationship between progenitors and their progenies at a single cell resolution. Currently scLT methods have been developed to be used in tandem with transcriptomical analysis (DARLIN, GESTALT, Polylox among many others). I worked on developing a novel scLT model, iFlpscLineage.

## iFlpscLineage

iFlpscLineage uses DNA barcodes that, through recombinase activity, are able to produce unique Barcode configurations that can be traced through the cells progeny to the moment where recombination ocurred. Under recombinase activity, the fixed DNA sequence undergoes either excision or deletion events. What sets each event to occur is mediated by the orientation of the recombination sites. If they face the same direction, the flanked DNA sequence will be excised. Alternatively, if they face opposite directions, inversion of the flanked sequence occurs.iFlpscLineage is comprised of three separate barcode arrays, arranged in tandem within the Rosa 26 locus. Independence of each array is ensured by using different recombination sites in each array, thereby impeding inter-array recombinations.

![iFlpscLineage](/images/iFlpscLineage.png)

### Experimental design

In order to retrieve both the transcriptomical and scLT profile of each sequenced cell, we used the 10x Chromium platform for performing scRNAseq followed by two independent amplification steps. Illumina sequencing would give information about the cells transcriptome and we used Nanopore sequencing for full lenght transcript sequencing of the iFlpscLineage barcodes.

![iFlpscLineage_protocol](/images/iFlpscLineage_Protocol.png)

### Computational pipeline

Once the sequencing data was produced the following pipeline was developed to merge Illumina and Nanopore sequenced reads onto one data object. The merge was followed by various QC steps aiming to work with:

1. Cells that passed standard quality controls
2. Cells from which both transcriptome and Barcode were retrieved
3. Cells expressing unique barcodes. These barcodes must assure that no other cell reached the same barcode through paralell means (dubbed barcode homoplasy)

![iFlpscLineage_pipeline](/images/iFlpscLineage_pipeline.png)
