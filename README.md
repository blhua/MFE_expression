# MFE_expression
Python script to calculate minimum free energy (MFE) of RNA secondary structures at the start codons of every gene within a NCBI nucleotide sequence and compare to relative protein expression in HEK-293 cells. iBAQ mass spectrometry HEK-293 data was sourced from [Geiger et al. Mol Cell. 2012.](https://www.mcponline.org/article/S1535-9476(20)30500-4/fulltext)
RNA sequences are extracted from defined regions relative to the start codon: -30:+30, +1:+60, and +31:+90 and the MFE of each sequence is calculated using RNAfold

![STR3030_scatter](https://github.com/blhua/MFE_expression/assets/66856632/56f7109b-b739-45c1-a322-0492cfc0488c)
![STR3030_kde](https://github.com/blhua/MFE_expression/assets/66856632/3d2a2fc8-7634-4e2c-81cd-0bfe15604de5)

## Prerequisites

RNAfold must be installed from [ViennaRNA Package 2](https://www.tbi.univie.ac.at/RNA/)

## Usage
```
python3 MFE_expression.py [accession_id]
```
## Inputs
- accession_id: An NCBI Nucleotide database accession number
  
## Outputs
- `$(accession_id}`.gb: Genbank file downloaded from NCBI using the accession id provided
- **HEK_iBAQ_with_RNA_info_and_MFE.csv**: Comma-separated values of each gene product from the genbank file merged with MFE values
- **MFE_-30:+30_scatter.jpg**: JPG file of the scatter plot of average expression (iBAQ average) vs. MFE of the -30:+30 region for each gene
- **MFE_+1:+60_scatter.jpg**: JPG file of the scatter plot of average expression (iBAQ average) vs. MFE of the +1:+60 region for each gene
- **MFE_+31:+90_scatter.jpg**: JPG file of the scatter plot of average expression (iBAQ average) vs. MFE of the +31:+90 region for each gene
- **MFE_-30:+30_kdeplot.jpg**: JPG file of the KDE plot of average expression (iBAQ average) vs. MFE of the -30:+30 region for each gene
- **MFE_+1:+60_kdeplot.jpg**: JPG file of the KDE plot of average expression (iBAQ average) vs. MFE of the +1:+60 region for each gene
- **MFE_+31:+90_kdeplot.jpg**: JPG file of the KDE plot of average expression (iBAQ average) vs. MFE of the +31:+90 region for each gene
