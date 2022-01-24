# Tbud

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4173415.svg)](https://doi.org/10.5281/zenodo.4173415)

This repository contains all R code and datasets to reproduce the analysis used for the paper *Toward better phenology analyses based on bud traits and energy budget*.

A detailed description of the scope, method, and algorithms is given by the RMarkdown file `./R_src/Tbud_perspective.Rmd`, or the compiled pdf `./R_src/Tbud_perspective.pdf`. 

To get this repository, change into a suitable directory and clone by
```bash
cd <suitable-working-direcotry>
git clone https://github.com/mpeaucelle/Tbud.git .
```

To render the RMarkdown file, execute all code, and thereby fully reproduce all analysis, simply enter the following command in RStudio (after making sure the 'rmarkdown' package is installed):
```r
install.packages("rmarkdown", type = "source")
rmarkdown::render_site()
```

Please also note the data use policy, described in `./R_src/Tbud_perspective.Rmd` and `./LICENSE`. When using this code, please cite the paper:
Peaucelle, Penuelas and Verbeeck (2022), Toward better phenology analyses based on bud traits and energy budget, *under review*.
This code is also available on [Zenodo](https://doi.org/10.5281/zenodo.4173415), where the latest version corresponds to tag `v2.0.0 `.

Marc Peaucelle, 17.01.2022
