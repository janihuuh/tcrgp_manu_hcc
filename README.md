# Analysis of HBV-specificites at the single-cell level in HCC-patients with TCRGP

In order to reproduce the analysis in the publication:

### First, clone this repository

```
git clone https://github.com/janihuuh/tcrgp_manu_hcc
```

### Download the scRNAseq and TCRab data
This is from GEO. We need the GSE98638_HCC.TCell.S5063.count.txt.gz file, rename it zheng_count.txt and add it into the data folder.

[(Link)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98638)

The TCR&alpha;&beta; data is already provided, downloaded from the Supplementary files from the original Cell publication (Zheng et al., *Cell* 2017)

### Then run the ``` analyze_perc085.R ``` in the same folder

```
Rscript src/analyze_perc085.R
```

## References
Jokinen et al., biorxiv 2019
Zheng et al., Cell 2017
Chang et al., Sci Immunol 2019
