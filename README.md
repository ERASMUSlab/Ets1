# Ets1 bioinformatics analysis user codes
The paper focuses on how chromatin architecture dynamics during osteoblast differentiation affect osteogenic genes, and on the mechanisms by which chromatin architecture dynamics occur in osteoblasts, and is published in AAAAA.

Osteogenesis is governed by the coordinated interplay between transcriptional networks and epigenetic regulation. Beyond linear gene control, emerging work highlights three-dimensional (3D) chromatin architecture—particularly the segregation into transcriptionally active A and repressive B compartments—as a key regulatory layer that shapes lineage-specific gene expression. Although compartment reorganization accompanies differentiation, the upstream mechanisms that drive this remodeling during osteoblast specification remain insufficiently defined. Lineage commitment proceeds through discrete transition points at which pioneer transcription factors access compact chromatin, displace nucleosomes, and recruit chromatin-modifying enzymes. Studies of Pax7, MyoD, and ETV2 demonstrate that pioneers regulate both chromatin accessibility and 3D genome organization, establishing an epigenetically primed landscape. Epigenetic priming, including H3K27ac pre-marking of enhancers, helps determine responsiveness to osteogenic cues such as BMP2 and Wnt3a, with committed cells (e.g., MC3T3-E1) showing robust induction. We hypothesize that Ets1 acts as a pioneer in early osteoblast differentiation, binding nucleosome-occluded regions, recruiting p300, depositing H3K27ac, and promoting B-to-A shifts. Integrating single-cell ATAC-seq, Hi-C compartment analysis, and scPRINTER footprinting in MC3T3-E1, we find that Ets1 precedes differentiation, increases local accessibility, primes enhancers with p300, and drives compartment-level reprogramming, positioning Ets1 as an upstream coordinator of early osteogenic fate.

![Alt text](./FIG/P1.png "Ets1")

# Input files
![Alt text](./FIG/P2.png "Ets1")

* Bulk RNA-seq at DAY0, 4, 7, 10 with 4 replicates
> total RNA was extracted from MC3T3-E1 cells cultured under osteogenic induction conditions using the RNeasy Mini Kit (Qiagen). Only high-quality RNA samples with RNA integrity number (RIN) values above 7.0, as determined by TapeStation analysis (Agilent), were selected for library preparation. Libraries were constructed using the TruSeq Stranded mRNA Sample Prep Kit (Illumina), following poly-A mRNA enrichment with oligo-dT magnetic beads. Purified mRNAs were thermally fragmented and reverse transcribed to generate first-strand cDNA, followed by second-strand synthesis, end repair, A-tailing, adaptor ligation, and PCR amplification. Library quality was assessed by TapeStation and quantified by qPCR prior to sequencing. Paired-end sequencing was performed on either the NovaSeq or NextSeq 2000 platform (Illumina), depending on the differentiation time point.
* ChIP-seq at GM, OM about H3K27ac, H3K4me1, H3K4me3
> ChIP assays were performed using the SimpleChIP® Enzymatic Chromatin IP Kit (Cat# 9003; Cell Signaling Technology, USA) according to the manufacturer’s instructions. Briefly, MC3T3-E1 cells were cross-linked, washed, harvested, and stored at –80°C until use. For chromatin preparation, nuclei were isolated and digested with micrococcal nuclease to obtain DNA fragments of approximately 150–900 bp, followed by mild sonication. After confirming chromatin quality, chromatin was immunoprecipitated with antibodies against H3K4me3 (Cat# ab8580, Abcam, , UK), H3K27ac (Cat# ab4729, Abcam, , UK), and H3K4me1 (Cat# ab8895, Abcam, , UK) using Protein G magnetic beads. After sequential washes with low-salt and high-salt buffers, antibody-bound chromatin was eluted and cross-links were reversed. DNA was purified using spin columns provided in the kit and stored at –20°C until analysis. ChIP DNA was subjected to library preparation using the NEBNext Ultra II DNA Library Prep Kit for Illumina (Cat# E7645, New England Biolabs, USA). The concentration and size distribution of the libraries were evaluated using a TapeStation system (Agilent Technologies) and a Qubit Fluorometer (Thermo Fisher Scientific). Libraries were subsequently sequenced using 150 bp paired-end reads (Theragen Bio Institute, Suwon, Korea).
* Cut&Run-seq at GM, OM about Ets1, p300, H3K27ac
> CUT&amp;RUN were performed using the ChIC/CUT&amp;RUN Assay Kit (Active Motif, Cat# 53180) according to the manufacturer’s instructions. Briefly, MC3T3-E1 cells were prepared under two conditions: undifferentiated cells and cells differentiated for 4 days under osteogenic conditions. Cells (5 × 10⁵ per condition) were harvested, bound to Concanavalin A-coated magnetic beads and incubated overnight at 4°C with primary antibodies against ETS1 (Cat# 14069, Cell Signaling Technology, Inc., , USA) and p300 (Cat# 05-257, Sigma-Aldrich, USA). After three washes to remove unbound antibodies, pAG-MNase was added, and chromatin digestion was performed. Released DNA fragments were purified and used for library construction with the NEBNext Ultra II DNA Library Prep Kit for Illumina (Cat# E7645, New England Biolabs,USA). Library quality was assessed using a TapeStation system (Agilent Technologies) and a Qubit Fluorometer (Thermo Fisher Scientific). Libraries were sequenced using 150 bp paired-end reads on either a NovaSeq 6000 platform (Theragen Bio Institute, Suwon, Korea) or a NextSeq 2000 platform (Dental- Multiomics Center, Seoul National University, Seoul, Korea).
* Bulk ATAC-seq at GM, OM with 4 replicates
> ATAC-seq was performed following a well-established protocol (Buenrostro et al., 2013), incorporating slight modifications to accommodate experimental conditions. In briefly, MC3T3-E1 cells (50,000 cells per sample) were harvested and lysed to isolate nuclei, followed by transposition using the Illumina TD enzyme (Cat# 20034197, illumine). Transposed DNA was purified using the DNA Clean &amp; Concentrator-5 kit (Cat#D4004, Zymo Research). Libraries were PCR-amplified using NEBNext High-Fidelity 2× PCR Master Mix with Nextera index primers. Amplification was monitored by qPCR using SYBR Green (cat#RR420A, Takara, Japan) to determine the optimal cycle number. Final libraries were purified. Paired- end sequencing (2 × 150 bp) was performed on a HiSeq platform (Illumina) by Macrogen, Inc. (Seoul, South Korea).
* scATAC-seq about osteoblast linage cell
* Micro-C at GM, OM with 3 replicates
> Micro-C libraries were generated using the Dovetail™ Micro-C Kit (Dovetail Genomics, Scotts Valley, CA, USA) according to the manufacturer&#39;s protocol. Briefly, MC3T3-E1 cells were subjected to sequential cross-linking with disuccinimidyl glutarate (DSG) and formaldehyde to preserve chromatin architecture. Following cross-linking, nuclei were isolated and chromatin was digested with micrococcal nuclease (MNase) to achieve predominantly fragments. Proximity ligation of digested chromatin fragments was then performed, and cross-links were subsequently reversed. Purified DNA was used for library preparation, including end repair, adaptor ligation, and limited-cycle PCR amplification. Final libraries were size- selected to enrich for fragments between 350 bp and 1 kb. Library quality and concentration were assessed using a TapeStation system (Agilent Technologies) and a Qubit Fluorometer (Thermo Fisher Scientific). Sequencing was carried out at Macrogen (Seoul, Republic of Korea) using the Illumina HiSeq X Ten platform to generate 150 bp paired-end reads.

# Processed datasets

| | Figures | link | Format | 
| --- | --- | --- | --- |
| Bulk RNA-seq expression counts | Fig. 1 | [link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | txt | 

# Main Figures

## Figure 1

|  | Figures | link | 
| --- | --- | --- | 
| Principal component analysis (PCA) of RNA-seq profiles at four time points between day 0, 4, 7 and 10. Each dot represents an individual biological replicate. | Fig. 1b | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| Upset plot visualizing differential gene expression profiles between consecutive time points: undifferentiated condition Day 0 to Day 4 after osteogenic induction (stage1), Day 4 to Day 7 (stage2), and Day 7 to Day 10 (stage3). DEGs were defined as those with adjusted p-value < 0.0001 and fold change 1.5. | Fig. 1c | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| Overrepresentation enrichment analysis (ORA) using Gene Ontology (GO) terms was performed on DEGs identified between undifferentiated stage1, stage2 and stage3. | Fig. 1d | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| Gene set enrichment analysis (GSEA) results visualizing positive enrichment (red) and negative enrichment (blue) during the differentiation phase from stage1. | Fig. 1e | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) |
| Scatter plot showing concordance between histone signals and gene expression changes between day 0 and day 4. Fold changes in ChIP-seq fragment enrichment at predefined regulatory elements were compared with gene expression changes. Genes were classified as matching (green) when regulatory element activity and expression changed in the same direction, or unmatching (red) otherwise. | Fig. 1i | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) |
          

# Citation

If you use this codes in your work, please cite: https://


# Contact information

* Author: 
* Affiliation: Department of Molecular Genetics, School of Dentistry and Dental Research Institute, Seoul National University, Seoul 08826, Republic of Korea
* Email: carpediemwj@snu.ac.kr, kitae@snu.ac.kr
