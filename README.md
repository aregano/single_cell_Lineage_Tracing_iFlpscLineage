# single_cell_Lineage_Tracing_iFlpscLineage
This repository contains a snippet of the analysis I performed to validate the scLT technology iFlpscLineage that can be used for clonal analysis coupled with transcriptomical data. This work was carried during my PhD thesis at the Molecular Genetics of Angiogenesis at CNIC.

## single cell Lineage Tracing (scLT)

single cell Lineage Tracing seeks to follow the division history of a cell in order to establish the relationship between progenitors and their progenies at a single cell resolution. Currently scLT methods have been developed to be used in tandem with transcriptomical analysis (DARLIN, GESTALT, Polylox among many others). I worked on developing a novel scLT model, iFlpscLineage.

## iFlpscLineage

iFlpscLineage uses DNA barcodes that, through recombinase activity, are able to produce unique Barcode configurations that can be traced through the cells progeny to the moment where recombination ocurred. Under recombinase activity, the fixed DNA sequence undergoes either excision or deletion events. What sets each event to occur is mediated by the orientation of the recombination sites. If they face the same direction, the flanked DNA sequence will be excised. Alternatively, if they face opposite directions, inversion of the flanked sequence occurs.iFlpscLineage is comprised of three separate barcode arrays, arranged in tandem within the Rosa 26 locus. Independence of each array is ensured by using different recombination sites in each array, thereby impeding inter-array recombinations.

![screenshot](/images/iFlpscLineage.png)


