# CilioGenics is an integrated method for predicting the ciliary genes 

To more accurately predict ciliary genes, CilioGenics combines several approaches, including as single-cell RNA sequencing, protein-protein interactions (PPIs), comparative genomics, transcription factor (TF)-network analysis, and text mining.

The current repository contains codes for each method, and representative codes are presented below.

#### **Score calculations**
 
To calculate score of each section (interactions, scRNA-seq, motifs, text mining), the corresponding R script can be run. To calculate all scores, or final CilioGenics score, `ciliogenics_scores.R` script can be run.

Please note that `interaction_scores.R` script will download a relatively large file (\~6.7 GB) to generate scores.

To clone the repo and generate ciliogenics scores, open terminal and run:

```         
git clone https://github.com/thekaplanlab/CilioGenics_Analysis
cd CilioGenics_Analysis

R ciliogenics_scores.R
```
## CilioGenics : Step-by-Step Tutorial 

**Step 1 –** 
Go to (https://ciliogenics.com/?page=Home)

When you go to the website, you will see 1: “Gene-Search” bar and 2: “Explore the gene list”.

![1](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/25293a61-d0cf-4798-853a-6d308c45ae50)

**Let’s go with option 1:**

If you search for a gene (e.g ARL13B)

![Adsız tasarım(5)](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/9f4985dd-5c09-4df3-bd58-62a6e86e03e3)

You will see a box for “Gene info” and another box for “CilioGenics scores for each category”

![FİGURE 3](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/37867f5b-a1fc-4223-895a-9fbde75f376c)

**Let’s go with the second option:**

Click “Explore the gene list” and you will see a gene list. You can calibrate the range of the score with the help of the scale bar under each title. 

![Adsız tasarım(6)](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/efb9d930-7622-474f-8797-4c707f95443f)

When you go to the “Explore data” on the left sidebar, you will see some heads that the other functions or properties belong to CilioGenics. 

![Adsız tasarım(7)](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/265c7489-ba28-4493-835d-415c6f75a81e)

**Click the “Phylogenetic Analysis”**

There will be a display of an interactive heatmap depicting each cluster and illustrating the distribution of genes between ciliated and non-ciliated cells.

An interactive heatmap illustrating the distinctions between ciliated and non-ciliated cells will be displayed. Additionally, you can view the list of genes in the cluster located below on the page.

![figure8](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/4c3672e5-3009-4250-b822-014ca50899cb)

![figure 9](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/efdcb645-c63f-40ad-b859-f33a4e807844)

**When you click the “Single cell cluster”, you will see a selection bar**

Choose one of the options, and you'll be presented with the UMAP plot corresponding to the selected tissue. Additionally, there is another selection bar adjacent to the UMAP plot that can generate the DotPlot featuring the chosen gene. Finally, at the bottom of the page, you can observe the results of the selected cluster's differential analysis.

![figure 10](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/cefe827b-71f1-4d6b-9a0a-42448e08ea26)

For the DotPlot: 

![figure 11](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/12113891-b7a5-4a3e-bccf-12abafd6702d)

![figure 123](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/3c51afdc-687f-4d30-a22d-f3d4a68e08c6)

For the differential expression analysis: 

![figure 12](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/97318228-89cb-409a-a48a-1cc7f4c54b98)

**Click the “Publications”:**

Users have the opportunity to examine a list of ciliary candidate genes collected from publications.

![figure 13](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/1dd129a0-c92e-457c-95d9-adb449cb8068)

![figure 14](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/aade82f3-084a-4606-a6c7-b7e11cca17b4)

**Click the “Motifs”:**

You will see the selection bar for Motifs, then, after the selection motifs you will see the motifs information and other gene information. 

![figure 15](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/ba26eebc-abd3-4e6f-b0f5-41fa4047316b)

![figure 16](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/50c540e7-666a-4812-8c12-669a677134ce)

**Click the “Protein Atlas”:**

Users can explore the compilation of genes in the Human Protein Atlas databases, obtained through the use of the keywords Cilia, Cilium, Centrosome, Flagella, and Flagellum. 

 ![figure 17](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/2a741027-b57f-4b6c-a200-dc244a2245cd)



