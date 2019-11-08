# Entropy Filtering for Multiple-sequence alignments
## About
Project for the Applied Bioinformatics course taken by participants of [**Medbioinfo**](http://www.medbioinfo.se/)

This workflow is for trimming noisy regions from MSA using measurements of information entropy of said position. 

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
git clone https://github.com/patruong/appliedBioinformaticsMSA                   
```

## Running
Move association files to `input/` and add filename to `python test = []` to the `Snakefile`
Now just run `snakemake` like in your project folder.

## Output
A HTML repport can be found at `results/{input_file}.html`.

