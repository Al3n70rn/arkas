#arkas[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.54752.svg)](http://dx.doi.org/10.5281/zenodo.54752)[![Travis-CI Build Status](https://travis-ci.org/RamsinghLab/arkas.svg?branch=master)](https://travis-ci.org/RamsinghLab/arkas)  
RNAseq analysis, from raw reads to pathways, typically in a few minutes.  
(Mostly by wrapping [Kallisto](http://pachterlab.github.io/kallisto/) and caching everything we possibly can.)  
Also a sandbox for wacky new ideas about RNA methylation,  
functional analysis of differential transcript abundance,  
and various other sorts of related malarkey (run the demo!).  


```R
suppressPackageStartupMessages(library(artemisData))
## this assumes Kallisto is in ~/bin/
demo("example", package="arkas")
```

![repeat expression](demo/example.png "Plot generated from example code")    
