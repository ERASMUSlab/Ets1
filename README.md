# Ets1 bioinformatics analysis user codes
The paper focuses on how chromatin architecture dynamics during osteoblast differentiation affect osteogenic genes, and on the mechanisms by which chromatin architecture dynamics occur in osteoblasts, and is published in AAAAA.

Osteogenesis is governed by the coordinated interplay between transcriptional networks and epigenetic regulation. Beyond linear gene control, emerging work highlights three-dimensional (3D) chromatin architecture—particularly the segregation into transcriptionally active A and repressive B compartments—as a key regulatory layer that shapes lineage-specific gene expression. Although compartment reorganization accompanies differentiation, the upstream mechanisms that drive this remodeling during osteoblast specification remain insufficiently defined. Lineage commitment proceeds through discrete transition points at which pioneer transcription factors access compact chromatin, displace nucleosomes, and recruit chromatin-modifying enzymes. Studies of Pax7, MyoD, and ETV2 demonstrate that pioneers regulate both chromatin accessibility and 3D genome organization, establishing an epigenetically primed landscape. Epigenetic priming, including H3K27ac pre-marking of enhancers, helps determine responsiveness to osteogenic cues such as BMP2 and Wnt3a, with committed cells (e.g., MC3T3-E1) showing robust induction. We hypothesize that Ets1 acts as a pioneer in early osteoblast differentiation, binding nucleosome-occluded regions, recruiting p300, depositing H3K27ac, and promoting B-to-A shifts. Integrating single-cell ATAC-seq, Hi-C compartment analysis, and scPRINTER footprinting in MC3T3-E1, we find that Ets1 precedes differentiation, increases local accessibility, primes enhancers with p300, and drives compartment-level reprogramming, positioning Ets1 as an upstream coordinator of early osteogenic fate.

![Alt text](./FIG/P1.png "Ets1")

# Input files
![Alt text](./FIG/P2.png "Ets1")

* Bulk RNA-seq at DAY0, 4, 7, 10 with 4 replicates
* ChIP-seq at GM, OM about H3K27ac, H3K4me1, H3K4me3
* Cut&Run-seq at GM, OM about Ets1, p300, H3K27ac
* Bulk ATAC-seq at GM, OM with 4 replicates
* scATAC-seq about osteoblast linage cell
* Micro-C at GM, OM with 3 replicates

# Processed datasets

| | link | Format | 
| --- | --- | --- |
| Bulk RNA-seq expression counts | [link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | txt | 

# Main Figures

## Figure 1

|  | Figures | link | 
| --- | --- | --- | 
| Principal component analysis (PCA) of RNA-seq profiles at four time points between day 0, 4, 7 and 10. Each dot represents an individual biological replicate. | Fig. 1b | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
          

# Citation

If you use this codes in your work, please cite: https://


# Contact information

* Author: 
* Affiliation: Department of Molecular Genetics, School of Dentistry and Dental Research Institute, Seoul National University, Seoul 08826, Republic of Korea
* Email: carpediemwj@snu.ac.kr, kitae@snu.ac.kr
