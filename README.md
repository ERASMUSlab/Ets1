<img width="468" height="19" alt="image" src="https://github.com/user-attachments/assets/04af26e7-f7ca-4a16-ad42-91d0726cac5b" /><img width="468" height="19" alt="image" src="https://github.com/user-attachments/assets/9af60928-b19d-4461-8e76-b309216efbbf" /># Ets1 bioinformatics analysis user codes
The paper focuses on how chromatin architecture dynamics during osteoblast differentiation affect osteogenic genes, and on the mechanisms by which chromatin architecture dynamics occur in osteoblasts, and is published in AAAAA.

![Alt text](./FIG/Ets1_1.png "Ets1")
> Osteogenesis is governed by the coordinated interplay between transcriptional networks and epigenetic regulation. Beyond linear gene control, emerging work highlights three-dimensional (3D) chromatin architecture—particularly the segregation into transcriptionally active A and repressive B compartments—as a key regulatory layer that shapes lineage-specific gene expression. Although compartment reorganization accompanies differentiation, the upstream mechanisms that drive this remodeling during osteoblast specification remain insufficiently defined. Lineage commitment proceeds through discrete transition points at which pioneer transcription factors access compact chromatin, displace nucleosomes, and recruit chromatin-modifying enzymes. Studies of Pax7, MyoD, and ETV2 demonstrate that pioneers regulate both chromatin accessibility and 3D genome organization, establishing an epigenetically primed landscape. Epigenetic priming, including H3K27ac pre-marking of enhancers, helps determine responsiveness to osteogenic cues such as BMP2 and Wnt3a, with committed cells (e.g., MC3T3-E1) showing robust induction. We hypothesize that Ets1 acts as a pioneer in early osteoblast differentiation, binding nucleosome-occluded regions, recruiting p300, depositing H3K27ac, and promoting B-to-A shifts. Integrating single-cell ATAC-seq, Hi-C compartment analysis, and scPRINTER footprinting in MC3T3-E1, we find that Ets1 precedes differentiation, increases local accessibility, primes enhancers with p300, and drives compartment-level reprogramming, positioning Ets1 as an upstream coordinator of early osteogenic fate.

# Input files
![Alt text](./FIG/Ets1_2.png "Ets1")

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

| | Figures | Format | link | 
| --- | --- | --- | --- |
| Bulk RNA-seq expression counts | Fig. 1 | txt | [link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) |


# Main Figures

## Figure 1

|  | Figures | Link | 
| --- | --- | --- | 
| Principal component analysis (PCA) of RNA-seq profiles at four time points between day 0, 4, 7 and 10. Each dot represents an individual biological replicate. | Fig. 1b | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| Upset plot visualizing differential gene expression profiles between consecutive time points: undifferentiated condition Day 0 to Day 4 after osteogenic induction (stage1), Day 4 to Day 7 (stage2), and Day 7 to Day 10 (stage3). DEGs were defined as those with adjusted p-value < 0.0001 and fold change 1.5. | Fig. 1c | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| Overrepresentation enrichment analysis (ORA) using Gene Ontology (GO) terms was performed on DEGs identified between undifferentiated stage1, stage2 and stage3. | Fig. 1d | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| Gene set enrichment analysis (GSEA) results visualizing positive enrichment (red) and negative enrichment (blue) during the differentiation phase from stage1. | Fig. 1e | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) |
| Scatter plot showing concordance between histone signals and gene expression changes at stage1. Fold changes in ChIP-seq fragment enrichment at predefined regulatory elements were compared with gene expression changes. Genes were classified as matching (green) when regulatory element activity and expression changed in the same direction, or unmatching (red) otherwise. | Fig. 1i | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) |


## Figure 2

| | Figures | Link | 
| --- | --- | --- | 
| Violin plots showing the enrichment of histone activation marks (H3K27ac, H3K4me1, H3K4me3) in A and B compartments at GM and OM. Blue represents the A compartment; yellow represents the B compartment. | Fig. 2b | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| Micro-C interaction matrix showing changes in chromatin interactions during differentiation. Red dots indicate regions with increased interactions at OM, while blue dots indicate regions with stronger interactions at GM. | Fig. 2c | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| Donut charts showing the proportions of A (blue) and B (yellow) compartments at GM and OM. | Fig. 2d | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) |
| Heatmap showing compartment dynamics at stage1, with genomic regions classified as Shift, Strength, Stable, or Subtle change. Regions were included if there was at least 50% overlap between compartments at the two time points (GM and OM). A signal increase of ≥1.5-fold at day 4 was defined as Strength, minimal change as Stable, a signal reversal with a change of ≤1.5-fold as Subtle change, and a reversal with a change of ≥1.5-fold as Shift. Blue indicates A compartments, and yellow indicates B compartments. | Fig. 2e | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) |
| Gene ontology (GO) overrepresentation analysis of genes located in regions undergoing compartment switching. Each dot represents the bin containing the gene’s regulatory region, arranged by compartment signal at each time point. GO terms associated with regions shifting from B to A compartments are shown in blue (F), while those shifting from A to B compartments are shown in yellow (G). | Fig. 2f-g | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) |


## Figure 3

| | Figures | Link | 
| --- | --- | --- | 
| Donut charts indicating the directional preference of compartment shift among osteogenic genes identified in Supplementary Figure 3B. Red indicates regions shifting toward the A compartment at stage1, while blue indicates regions shifting toward the B compartment. | Fig. 3a | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| The promoters of genes directed toward each compartment shown in Figure 3A were represented as box plots of H3K27ac levels. The left plot shows the changes in H3K27ac levels at stage 1, the middle plot shows the levels in GM, and the right plot shows the levels in OM. | Fig. 3b | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| Box plots showing the H3K27ac levels in GM across four compartment shift categories between day 0 and day 4: A-to-A (blue), B-to-B (yellow), A-to-B (light blue), and B-to-A (light brown). | Fig. 3c | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) |
| Hi-C interaction map and histone activation profiles for chromosome 5, which exhibited the highest frequency of A/B compartment shfifts. To investigate whether regions undergoing compartment shifts show distinct histone activation levels at GM, we visualized compartment structure using Hi-C matrices overlaid with average H3K27ac, H3K4me1, and H3K4me3 signal per bin (lower left). The top 10% and bottom 10% of active histone signals are marked in red and blue, respectively (upper right). The longest B-to-A shifted bin is highlighted in light brown and coincides with a high-activation bin, while the longest A-to-B shifted bin, shaded in sky blue, corresponds to a bin with low histone activation. | Fig. 3d | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| Hi-C interaction map and histone activation profiles for chromosome 15, which exhibited the second highest frequency of A/B compartment shifts. To examine whether regions undergoing compartment shifts display distinct histone activation levels in GM, we visualized the compartment structure using Hi-C matrices overlaid with H3K27ac signals (lower left). The top 10% and bottom 10% of histone activation signals are highlighted in red and blue, respectively (upper right). The bin containing the osteogenic gene Sp7 is shaded in light brown and coincides with a highly active region, whereas the bin containing the proliferation-related gene "" is shaded in sky blue and corresponds to a region with low histone activation. The right panel shows a magnified view of the H3K27ac, H3K4me1, and H3K4me3 levels in the bins containing the osteogenic gene Sp7 and the proliferation-related gene "". | Fig. 3e | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| Box plots showing histone activation levels (H3K27ac, H3K4me1, and H3K4me3) at GM across four compartment dynamics categories: A-to-A (blue), B-to-B (yellow), A-to-B (light blue), and B-to-A (light brown).| Fig. 3f | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 


## Figure 4

| | Figures | Link | 
| --- | --- | --- | 
| | Fig. 3a | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 


## Figure 5

| | Figures | Link | 
| --- | --- | --- | 
| The binding levels of p300 relative to Ets1 binding levels were represented as a line plot. The binding levels of Ets1 and p300 were divided into quartiles (Q1–Q10), and the degree of p300 enrichment was indicated by the intensity of green shading. | Fig. 5a | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| The directional changes of A/B compartments according to Ets1 and p300 binding levels were represented as a line plot. The bins of the entire chromosome were divided into quartiles (Q1–Q10) based on the binding levels of Ets1 and p300, and the compartment shift direction was plotted relative to 1. Values above 1 indicate a shift toward the A compartment, while values below 1 indicate a shift toward the B compartment. | Fig. 5b | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| Clustering of UMAP reveals 22 distinct compartment behavior clusters. | Fig. 5c | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| UMAP dimensionality reduction of 24,639 genomic bins based on compartment interaction changes from stage1 divided into quartiles (Q1–Q10) by compartment shift direction, coloring them to represent A/B compartment shift direction and extent. | Fig. 5d | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| The directional changes of A/B compartments in each cluster were represented as a line plot. The compartment shift direction was plotted relative to thresholds of 1.2 and 0.8. Values above 1.2 indicate clusters shifting toward the A compartment, values below 0.8 indicate clusters shifting toward the B compartment, and clusters with values between 0.8 and 1.2 were considered unchanged. | Fig. 5e | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| Based on compartment interaction changes at stage 1, bins were divided into quartiles (Q1–Q10) according to the direction of compartment shifts. The distribution ratios of Q10, representing the strongest shift toward the A compartment, and Q1, representing the strongest shift toward the B compartment, were then plotted as a line graph for each cluster. | Fig. 5f | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| Donut plots display Ets1 and p300 binding distributions across clusters. | Fig. 5G | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 


## Figure 7

| | Figures | Link | 
| --- | --- | --- | 
| Cell type annotation of primary osteoblast populations based on scATAC-seq chromatin accessibility profiles. | Fig. 7e | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| Pseudotime trajectory analysis ordering cells along a differentiation continuum from early progenitors (pseudotime = 0) to terminally differentiated osteoblasts (pseudotime > 10). | Fig. 7f | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 
| Gene activity score plots showing peak chromatin accessibility at bone marker gene loci Sp7 and Ibsp, predominantly in pre-osteoblast clusters. | Fig. 7g-h | [Link](http://147.47.56.90:8895/Ets1/4DN_official_202508/) | 


# Citation

If you use this codes in your work, please cite: https://


# Contact information

* Author: 
* Affiliation: Department of Molecular Genetics, School of Dentistry and Dental Research Institute, Seoul National University, Seoul 08826, Republic of Korea
* Email: carpediemwj@snu.ac.kr, kitae@snu.ac.kr
