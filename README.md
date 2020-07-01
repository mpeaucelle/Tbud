# Tbud

Add DOI here

This repository contains all R code and datasets to reproduce the analysis used for the paper *Plants are warming faster than climate*.

A detailed description of the scope, method, and algorithms is given by the RMarkdown file `Tbud_perspective.Rmd`. 

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

Please also note the data use policy, described in `./Tbud_perspective.Rmd` and `./LICENSE`. When using this code, please cite the paper Peaucelle, Penuelas and Verbeeck (2020) *under review*.
This code is also available on [Zenodo](https://zenodo.org/XXX), where the latest version corresponds to tag `XXX `.

Marc Peaucelle, 01.07.2020