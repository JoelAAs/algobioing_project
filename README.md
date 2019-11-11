# FDR-adjusting associations p-values grouped on BRITE hierachy  
# About
Project for Algorithmic Bioinformatics course taken by participants of [**Medbioinfo**](http://www.medbioinfo.se/)

This workflow will annotate associations provided with associated genes and [KEGG Funtional hierachies](https://www.genome.jp/kegg/kegg3b.html). Using the distribution of for each feature group it calculates q-values and output visualisation (manhattan plots, p-value distribution etc). 


## Prerequisites
+ `Python 3.7.3`
  + `Snakemake 5.7.4+`
  + `seaborn`
  + `pandas`
  + `matplotlib`

+ `R 3.5.2 +`
  + `biomaRt`
  + `tidyverse`
  + `magrittr`
  + `KEGGREST`
  + `reshape2`


## Setup
First clone this repo into desired location.
```
git clone https://github.com/JoelAAs/algobioing_project 
```

For python packages I'd reccomend installing them through Conda.

The R to install run:
```{r}
install.packages(c("tidyverse",  "magrittr", "reshape2", "BiocManager"))
BiocManager::install(c("biomaRt", "KEGGREST"))
``` 


## Running
Move association files to `input/` and add filename to `python test = []` to the `Snakefile`.  
Now just run `snakemake` like in your project folder. Observe that the header of the association file is assumed to be that of a plink `.assoc` file.

## Output
A HTML repport can be found at `results/{input_file}.html`.

