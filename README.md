## Reproducible Scripts for the Publication

Drost HG, Ó'Maoiléidigh DS, Gabel A, Bellstädt J, Ryan PT, Dekkers BJW, Bentsink L, Silva AT, Hilhorst H, Ligterink W, Wellmer F, Grosse I, and Quint M. (2015). __Molecular hourglass patterns in plants are disconnected from organogenesis and
body plan establishment__


## Performing Phylostratigraphy and Divergence Stratigraphy

The first step in performing phylotranscriptomic analyses is to
perform `Phylostratigraphy` and `Divergence Stratigraphy` to assign
each protein coding genes of a query organism an evolutionary age (Phylostratigraphy)
or degree of selection pressure (Divergence Stratigraphy). The resulting
`Phylostratigraphic Map` and `Divergence Map` of a query organism is then matched
with the corresponding transcriptome data covering the developmental process of interest.
A detailed introduction into phylotranscriptomics can be found [here](http://cran.r-project.org/web/packages/myTAI/vignettes/Introduction.html).

### Phylostratigraphy

[Phylostratigraphy](http://www.sciencedirect.com/science/article/pii/S0168952507002995) is a computational method to determine the evolutionary origin of protein coding genes
based on BLAST homology searches.

For this study, we loaded the proteome of `Arabidopsis thaliana` from [Phytozome](http://www.phytozome.net/) as follows:


Open the [R](http://cran.r-project.org/) Command Line Interface and type:

```r
# download the proteome of A. thaliana
download.file( url      = "ftp://ftp.jgi-psf.org/pub/compgen/phytozome/v9.0/Athaliana/annotation/Athaliana_167_protein.fa.gz", 
               destfile = "Athaliana_167_protein.fa.gz" )
```

The next step is to perfom [Phylostratigraphy](http://www.sciencedirect.com/science/article/pii/S0168952507002995) using the `A. thaliana` proteome as query. The following steps have to be done to retrieve a phylostratigraphic map for `A. thaliana`:
    
1) Make sure that BLAST (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/2.2.21/) is installed on your machine.

2) Download the sequence database <a href="http://msbi.ipb-halle.de/download/phyloBlastDB_Drost_Gabel_Grosse_Quint.fa.tbz">phyloBlastDB_Drost_Gabel_Grosse_Quint.fa</a> used for BLAST searches and unpack it (`tar xfvj phyloBlastDB_Drost_Gabel_Grosse_Quint.fa.tbz`).

3) Make sure that the header of your FASTA-files (e.g. Athaliana_167_protein.fa) fullfills the following specification:<br />
  <code>>GeneID | [organism_name] | [taxonomy]</code><br />
  Notice, the taxonomy begins after the node "Cellular organisms" e.g.
```{terminal}
>NP_146894.1 | [Aeropyrum pernix] | [Archaea; Crenarchaeota; Thermoprotei; Desulfurococcales; Desulfurococcaceae; Aeropyrum]
or
>YP_001514406.1 | [Acaryochloris marina MBIC11017] | [Bacteria; Cyanobacteria; Oscillatoriophycideae; Chroococcales; Acaryochloris; Acaryochloris marina]
or
>ATCG00500.1|PACid:19637947 | [Arabidopsis thaliana] | [Eukaryota; Viridiplantae; Streptophyta; Streptophytina; Embryophyta; Tracheophyta; Euphyllophyta; Spermatophyta; Magnoliophyta; eudicotyledons; core eudicotyledons; rosids; malvids; Brassicales; Brassicaceae; Camelineae; Arabidopsis]
```

4) Use the following command to start the [Perl script](https://github.com/HajkD/Active-maintenance-of-phylotranscriptomic-hourglasses/blob/master/createPsMap.pl)

```terminal
perl createPsMap.pl -i Athaliana_167_protein_with_new_Header.fa -d phyloBlastDB_Drost_Gabel_Grosse_Quint.fa -p BLAST_Athaliana 
                    -r athaliana_blast_results -t 30 -a 64             
Arguments:
-i,--input          input file name in FASTA format
-d,--database       BLAST sequence database name
-p,--prefix         Prefix for generated mysql-files containing BLAST results
-r,--resultTable    mysql table name
-t,--threshold      threshold for sequence length (Default 30 amino acids)
-a                  threads for BLAST searches
-e,--evalue         e-value threshold for BLAST 
```

### Divergence Stratigraphy

[Divergence Stratigraphy](http://mbe.oxfordjournals.org/content/early/2015/01/27/molbev.msv012.abstract) is a computational method to determine the degree of selection pressure acting on each
protein coding gene of a query organism against a reference organism.

A detailed tutorial on how to perform Divergence Stratigraphy for any pairwise genome comparison can
be found in the [Divergence Stratigraphy Vignette](https://github.com/HajkD/orthologr/blob/master/vignettes/divergence_stratigraphy.Rmd) of the [orthologr](https://github.com/HajkD/orthologr) package or in the [Advanced Phylotranscriptomics Analyses](http://cran.r-project.org/web/packages/myTAI/vignettes/Advanced.html) vignette of the [myTAI](http://cran.r-project.org/web/packages/myTAI/index.html) package.

Following steps are performed to obtain a standard divergence map for `A. thaliana` versus `A. lyrata`:

1) Orthology Inference using BLAST reciprocal best hit ("RBH") based on blastp

2) Pairwise global amino acid alignments of orthologous genes using the [Needleman-Wunsch algorithm](http://www.sciencedirect.com/science/article/pii/0022283670900574)

3) Codon alignments of orthologous genes using [PAL2NAL](http://www.bork.embl.de/pal2nal/)

4) dNdS estimation using [Comeron's method (1995)](http://link.springer.com/article/10.1007/BF00173196)

5) Assigning estimated dNdS values to divergence strata (deciles of all dNdS values)

When using the `divergence_stratigraphy()` function implemented in `orthologr` it is assumed that you have BLAST installed on your machine.


To obtain a `Divergence Map` for `A. thaliana` versus `A. lyrata` the following commands need to be passed to the [R](http://cran.r-project.org/) Command Line Interface:


a) Installing the [orthologr](https://github.com/HajkD/orthologr) package: 

```r

# install package 'orthologr' from: https://github.com/HajkD/orthologr
install.packages("devtools") # note for wondows installation see https://github.com/HajkD/orthologr for details
devtools::install_github("HajkD/orthologr", build_vignettes = TRUE, dependencies = TRUE)
library(orthologr)

# install Bioconductor base packages
source("http://bioconductor.org/biocLite.R")
biocLite()

# install package: Biostrings
biocLite("Biostrings")

# install package: S4Vectors
source("http://bioconductor.org/biocLite.R")
biocLite("S4Vectors")

```

b) Downloading the CDS files of `A. thaliana` and `A. lyrata`:

```r
# download the CDS file of A. thaliana
download.file( url      = "ftp://ftp.jgi-psf.org/pub/compgen/phytozome/v9.0/Athaliana/annotation/Athaliana_167_cds.fa.gz", 
               destfile = "Athaliana_167_cds.fa.gz" )
               
# download the CDS file of A. lyrata
download.file( url      = "ftp://ftp.jgi-psf.org/pub/compgen/phytozome/v9.0/Alyrata/annotation/Alyrata_107_cds.fa.gz", 
               destfile = "Alyrata_107_cds.fa.gz" )
```

c) Compute the `Divergence Map` of `A. thaliana` versus `A. lyrata`:

```r
library(orthologr)

# compute the divergence map of A. thaliana vs. A. lyrata
Ath_vs_Aly_DM <- divergence_stratigraphy(
                         query_file      = "Athaliana_167_cds.fa",
                         subject_file    = "Alyrata_107_cds.fa",
                         eval            = "1E-5", 
                         ortho_detection = "RBH",
                         comp_cores      = 1, 
                         quiet           = TRUE, 
                         clean_folders   = TRUE )
   

```


## Mapping Gene IDs

It is now assumed that the `Phylostratigraphic Map` and `Divergence Map` of interest and the corresponding gene expression data set are joined. For this purpose the `MatchMap()` function implemented in the [myTAI](http://cran.r-project.org/web/packages/myTAI/index.html) package can be used. See `?myTAI::MatchMap` for details.


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

First install and load the [myTAI](http://cran.r-project.org/web/packages/myTAI/index.html) package:

```r
install.packages("myTAI")
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





