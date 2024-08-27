# Pathmatrix simulation of iFlpscLineage
Due to the sheer amount of possible barcodes assemblies that iFlpscLineage can produce (6.12 billion/array), we first designed a in sillico simulation approach to get the most common barcodes and store them in the form
of a pathmatrix, where rows represent each barcode and columns represent each recombination step. Values would be the probability of generation of each barcode under that amount of recombination steps.
With this pathmatrix it can be established the Pgen of the barcodes retrieved in each experiment. However, it would be ideal to have a full pathmatrix that registers the probability of generation of each barcode present
in the system. Nevertheless, if iFlpscLineage’s barcode generation probabilities were to be successfully produced in silico, the matrix with all the 6,12 billion possible barcode frequencies based on 11 recombination steps
would not be able to be manipulated with ease with a regular desktop computer as I have estimated its size to be of 534GB.
