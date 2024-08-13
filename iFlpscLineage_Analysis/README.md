# Analysis of iFlpscLineage

My thesis work revolved around validating iFlpscLineage as a tool from which to retrieve Lineage information and carry scLT analysis couple with transcriptomics. We developed two main approaches for validating the technology which were based on recombinase activity.

## Constitutive Recombination

In this model the recombinase is permanently active once expressed.

**Pros:** easier detection of recombined barcodes

**Cons:** reduced Barcode variability

## Partial Recombination

In this model the recombinase activity is regulated. This regulation of recombinase activity comes through the use of a chimeric FlpO-ERT2 recombinase, which can only translocate inside the nucleus and edit the NA Barcode sequences when coupled with Tamoxifen. Therefore, recombinase activity can be regulated via the introduction of Tamoxifen in the system

**Pros:** Barcode variability that ensures unique clones are retrieved

**Cons:** Harder detection of Barcodes in each cell and harder retrieval of full length transcript barcodes

**Note:** As only Partial recombination barcodes give enough variability that should be the only viable approach when doing a biological assay that aims at producing a dataset that couples scLT and transcriptomics data through scRNAseq.
