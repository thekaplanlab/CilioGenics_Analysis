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

![resim](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/23f05436-40ec-4e2f-bd31-d063ee19879f) 

 
**Let’s go with option 1:**

If you search for a gene (e.g ARL13B)

![resim](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/4bdb9646-0db6-4b21-bd3b-dd80ba50ea9a)

 
You will see a box for “Gene info” and another box for “CilioGenics scores for each category”

 ![resim](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/db1ddb25-f29b-4553-a576-387d7765108f) 

**Let’s go with the second option:**

Click “Explore the gene list” and you will see a gene list. You can calibrate the range of the score with the help of the scale bar under each title. 

 ![resim](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/6273bfdb-4a13-405a-8880-8b64a3903643)

When you go to the “Explore data” on the left sidebar, you will see some heads that the other functions or properties belong to CilioGenics. 

 ![resim](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/3913e577-8cc1-4d79-b35c-1ca20453bb84)

**Click the “Phylogenetic Analysis”**

You will see an interactive heatmap that shows the clusters between ciliated and non-ciliated cells. Also, you can see the genes in the cluster list below part of the page. 

![resim](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/c1b30a28-8c91-47a9-b5ad-ade6fb41663b)

![resim](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/fddf719c-b42f-4c9a-9ee6-fb1905f72037)

**When you click the “Single Cell Clusters”, you will see a selection bar**

Select one of the options and you will see the UMAP plot for the selected tissue. Also, there is another selection bar next to the UMAP plot which can generate the DotPlot with the selected gene. Lastly, at the bottom of the page, you can see the differential analysis belonging to the selected cluster. 

![resim](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/a5954088-fd81-4ce8-978f-ba4b8b672e88)

Fort he DotPlot: 

![resim](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/98806d21-8e74-45a3-a0c6-de1874dc8799)

![resim](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/6e4ae348-8198-430c-87ef-f564583173ad)

For the differential expression analysis: 

![resim](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/5bb25f04-131a-4d80-a681-0a3cb072e6b3)

**Click the “Publications”:**

You will see the list of the publications and the number of the genes below. 

![resim](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/9ea66ed5-11d0-44be-ba8a-bd5f9aa9df9c)

![resim](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/0d3f90ee-f3d7-4ebf-9b72-2628a176eb99)


**Click the “Motifs”:**

You will see the selection bar for Motifs, then, after the selection motifs you will see the motifs information and other gene information. 

![resim](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/e9fb0e5b-9ca5-4c96-87fe-9f4e64d85662)4

![resim](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/3df7f90f-edc3-43f1-9e98-692d686995bc)

**Click the “Protein Atlas” :**

You will see the list of genes in the Human Protein Atlas databases by selecting the keywords Cilia, Cilium, Centrosome, Flagella, and Flagellum, respectively. 
 
![resim](https://github.com/thekaplanlab/CilioGenics_Analysis/assets/126083033/76479cbd-cef1-479f-9a27-f0d39d2aa199)


