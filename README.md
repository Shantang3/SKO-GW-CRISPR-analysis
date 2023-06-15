# Sinle-gene KnockOut Genome0wide CRISPR Screen Anaylysis Pipeline
This pepilene is used for downstream data analysis and visualization after gRNA count and beta score calculation by MAGeCK. It is mainly derived from MAGeCKFlute R packages source code (https://rdrr.io/bioc/MAGeCKFlute/f/) and r code from the Nature Protocol paper (Paper: https://www.nature.com/articles/s41596-018-0113-7; script: http://cistrome.org/MAGeCKFlute/demo/Figures/script/). The purpose of using this pipeline instead of the ready-to-use package is to enable more customized control over the analysis and figures. 

## 1. Reads mapping by MAGeCKflute
A detailed demonstration of how to install and run MAGeCKflute can be found at: 
[https://sourceforge.net/p/mageck/wiki/Home/](https://sourceforge.net/p/mageck/wiki/Home/)

## 2.Quality control:
Mainly focus on the statistics of sequecing, count distribution and sample correlation.

## 3. Beta score cauculation
Use beta score as a quantitative value for essentiality measurement.

## 4. Visualization
