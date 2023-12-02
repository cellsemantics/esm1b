Data 
The E.coli LTEE data consisting of varaints information was taken from the paper - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5788700/. 
The data was sourced from - https://github.com/benjaminhgood/LTEE-metagenomic
The E.coli LTEE data consisting of fitness information was taken from the paper - https://royalsocietypublishing.org/doi/10.1098/rspb.2015.2292

Analysis Overview 
In this analysis, the E. coli LTEE experiment dataset was utilized.
A total of 22,857 missense variants were reported out of which 18059 were successfully retrived and matched with their respective LLR score

Steps 
1) Preprocessing of the Data
All the DNA and Protein sequences for each gene were taken from NCBI - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/017/985/GCF_000017985.1_ASM1798v1/

2) Esm1b model
The pretrained model was sourced from - https://github.com/ntranoslab/esm-variants

3) Analysis
Script1.ipynb and Script2.R were employed for the analysis, while Script3.ipynb was utilized for generating plots of analysis.
[Note : Details of the steps are provided in the Script folder]
