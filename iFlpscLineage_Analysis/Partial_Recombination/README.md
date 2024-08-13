# Partial Recombination

Once iFlpscLineage was validated through the Constitutive recombination model, where it gave us valuable information of the expression levels of the barcodes, we moved onto trying Partial recombination models to understand the range of Barcode variability iFlpscLineage could offer us.

## Experimental design

![Experimental design](../../images/Experimental_design_Partial.png)

Organs from different adult mice were processed and cryostored, for later scRNASeq processing. On the experimental day cryovials were pooled and multiplexed using CMOs so all samples could be loaded onto the same 10x Chromium port.

## Barcoding Strategy

Barcode variability was tested using two allele models. The first one, CAG FlpOERT2/iFlpMosaics/ iFlpscLineage, hereafter referred as iFlpscL iFlpMosaics, was induced a P5 via single 4-hydroxytamoxifen (4OHT) injection, which prompts the activation of the FlpOERT2 driver that independently activates iFlpMosaics for fluorescent cell labelling, alongside iFlpscLineage. The second model, Cdh5-CreERT2/iSureHadCre/iFlpscLineage, referred as iFlpscL iSure-HadCre, underwent induction one week before processing via two tamoxifen injections, which in turn causes endothelial cell specific activation of CreERT2, initiating a genetic cascade with the activation of iSureHadCre that gives both permanent Tomato reporter expression and transient FlpO expression, that allows for partial recombination of iFlpscLineage.

![Barcoding](../../images/Partial.png)
