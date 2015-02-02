## Reproducible Scripts for the Publication

Drost HG, Ó'Maoiléidigh DS, Gabel A, Bellstädt J, Ryan PT, Dekkers BJW, Bentsink L, Silva AT, Hilhorst H, Ligterink W, Wellmer F, Grosse I, and Quint M. (2015). __Molecular hourglass patterns in plants are disconnected from organogenesis and
body plan establishment__ (In Review)


## Performing Phylostratigraphy and Divergence Stratigraphy

The first step in performing phylotranscriptomic analyses is to
perform `Phylostratigraphy` and `Divergence Stratigraphy` to assign
each protein coding genes of a query organism an evolutionary age (Phylostratigraphy)
or degree of selection pressure (Divergence Stratigraphy). The resulting
`Phylostratigraphic Map` and `Divergence Map` of a query organism is then matched
with the corresponding transcriptome data covering the developmental process of interest.
A detailed introduction into phylotranscriptomics can be found [here](http://cran.r-project.org/web/packages/myTAI/vignettes/Introduction.html).

### Phylostratigraphy

Phylostratigraphy is a computational method to determine the evolutionary origin of protein coding genes
based on BLAST based homology searches.

A detailed description on how to perform [Phylostratigraphy](http://www.sciencedirect.com/science/article/pii/S0168952507002995) for _Arabidopsis thaliana_ as well as a reproducible `Perl Script` to perorm Phylostratigraphy can be found [here](https://github.com/HajkD/Active-maintenance-of-phylotranscriptomic-hourglasses#performing-phylostratigraphy).

Reproducible Scripts to obtain a `Phylostratigraphic Map` for `Arabidopsis thaliana` can be 
found in the section [Performing Phylostratigraphy](https://github.com/HajkD/Active-maintenance-of-phylotranscriptomic-hourglasses#performing-phylostratigraphy).

### Divergence Stratigraphy

Divergence Stratigraphy is a computational method to determine the degree of selection pressure acing on each
protein coding gene of a query organism against a reference organism.

A detailed tutorial on how to perform Divergence Stratigraphy for any pairwise genome comparison can
be found in the [Divergence Stratigraphy Vignette](https://github.com/HajkD/orthologr/blob/master/vignettes/divergence_stratigraphy.Rmd) of the [orthologr](https://github.com/HajkD/orthologr) package or in the [Advanced Phylotranscriptomics Analyses](http://cran.r-project.org/web/packages/myTAI/vignettes/Advanced.html) vignette of the [myTAI](http://cran.r-project.org/web/packages/myTAI/index.html) package.

Reproducible Scripts to obtain a `Divergence Map` for `Arabidopsis thaliana` versus `Arabidopsis lyrata` can be 
found in the section [Performing Divergence Stratigraphy](https://github.com/HajkD/Active-maintenance-of-phylotranscriptomic-hourglasses#performing-divergence-stratigraphy).

## Mapping Gene IDs

It is now assumed that the `Phylostratigraphic Map` and `Divergence Map` of interest and the corresponding gene expression data set are joined. For this purpose the `MatchMap()` function implemented in the [myTAI](http://cran.r-project.org/web/packages/myTAI/index.html) package can be used. See `?myTAI::MatchMap` for details. However, the `MatchMap()` function can only deal with identical gene ids present in the Phylo/Divergence-Maps and the corresponding gene expression set.


## Reading Datasets

The following script allows you to read `Supplementary datasets S1-S4`.

```r
# install.packages(gdata)
library("gdata")

## read PhyloExpressionSets
Ath.PhyloExpressionSet.Flowering <- read.xls("Supplementary dataset S1.xls",sheet = 1)
Ath.PhyloExpressionSet.wholeSeed <- read.xls("Supplementary dataset S2.xls",sheet = 1)

## read DivergenceExpressionSets
Ath_Aly.DivergenceExpressionSet.Flowering <- read.xls("Supplementary dataset S3.xls",sheet = 1)
Ath_Aly.DivergenceExpressionSet.wholeSeed <- read.xls("Supplementary dataset S4.xls",sheet = 1)



```

## Generating Figures

First load the following package into the work space:

```r
# install.packages("myTAI")
library(myTAI)
```


### Figure 2

```r

svg("Fig2.svg",16.9,5)
par(mfrow=c(1,2))

# plot the TAI of A. thaliana Flowering
PlotPattern( ExpressionSet =  Ath.PhyloExpressionSet.Flowering,
             permutations  = 10000, 
             type          = "l", 
             lwd           = 9, 
             col           = "black",
             xlab          = "Ontogeny",
             ylab          = "TAI", 
             main          = "",
             cex.lab       = 1.2,
             cex.axis      = 1.2,
             las           = 1 )
            
par(xpd = TRUE)
legend("topleft",legend = expression(bold("B")), bty = "n", cex = 1.5,inset = c(-0.08,-0.15))

# plot the Relative Expression Profiles of Phylostrata during A. thaliana Flowering
PlotRE( ExpressionSet = Ath.PhyloExpressionSet.Flowering,
        Groups        = list(c(1:12)),
        legendName    = "PS",
        lty           = 1, 
        lwd           = 5 )

par(xpd = TRUE)
legend("topleft",legend = expression(bold("C")), bty = "n", cex = 1.5,inset = c(-0.08,-0.15))


dev.off()

```


### Figure 3


```r
# plot the TAI of A. thaliana Germination
svg("Fig3.svg",16.9,5)
par(mfrow=c(1,2))

PlotPattern( ExpressionSet = Ath.PhyloExpressionSet.wholeSeed,
             permutations  = 10000, 
             type          = "l", 
             lwd           = 9, 
             col           = "black",
             xlab          = "Ontogeny",
             ylab          = "TAI", 
             main          = "",
             cex.lab       = 1.2,
             cex.axis      = 1.2,
             las           = 1 )

par(xpd = TRUE)
legend("topleft",legend = expression(bold("B")), bty = "n", cex = 1.5, inset = c(-0.08,-0.15))


# plot the Relative Expression Profiles of Phylostrata during A. thaliana Germination
PlotRE( ExpressionSet = Ath.PhyloExpressionSet.wholeSeed,
        Groups        = list(c(1:12)),
        legendName    = "PS",
        lty           = 1, 
        lwd           = 5 )

par(xpd = TRUE)
legend("topleft",legend = expression(bold("C")), bty = "n", cex = 1.5,inset = c(-0.08,-0.15))


dev.off()

```


### Figure 4


```r

svg("Fig4.svg",16.9,10)
par(mfrow=c(2,2))

# plot the TDI of A. thaliana Flowering
PlotPattern( ExpressionSet = Ath_Aly.DivergenceExpressionSet.Flowering,
             permutations  = 10000, 
             type          = "l", 
             lwd           = 9, 
             col           = "black",
             xlab          = "Ontogeny",
             ylab          = "TDI", 
             main          = "",
             cex.lab       = 1.2,
             cex.axis      = 1.2,
             las           = 1 )

par(xpd = TRUE)
legend("topleft",legend = expression(bold("A")), bty = "n", cex = 1.5, inset = c(-0.08,-0.15))


# plot the Relative Expression Profiles of Divergence Strata during A. thaliana Flowering
PlotRE( ExpressionSet = Ath_Aly.DivergenceExpressionSet.Flowering,
        Groups        = list(c(1:10)),
        legendName    = "DS",
        lty           = 1, 
        lwd           = 5 )

par(xpd = TRUE)
legend("topleft",legend = expression(bold("B")), bty = "n", cex = 1.5,inset = c(-0.08,-0.15))

# plot the TDI of A. thaliana Germination
PlotPattern( ExpressionSet = Ath_Aly.DivergenceExpressionSet.wholeSeed,
             permutations  = 10000, 
             type          = "l", 
             lwd           = 9, 
             col           = "black",
             xlab          = "Ontogeny",
             ylab          = "TDI", 
             main          = "",
             cex.lab       = 1.2,
             cex.axis      = 1.2,
             las           = 1 )

par(xpd = TRUE)
legend("topleft",legend = expression(bold("C")), bty = "n", cex = 1.5, inset = c(-0.08,-0.15))


# plot the Relative Expression Profiles of Divergence Strata during A. thaliana Germination
PlotRE( ExpressionSet = Ath_Aly.DivergenceExpressionSet.wholeSeed,
        Groups        = list(c(1:10)),
        legendName    = "DS",
        lty           = 1, 
        lwd           = 5 )

par(xpd = TRUE)
legend("topleft",legend = expression(bold("D")), bty = "n", cex = 1.5,inset = c(-0.08,-0.15))



dev.off()

```


## Supplementary Figures

### Supplementary Figure S1 







