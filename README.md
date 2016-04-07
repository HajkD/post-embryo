## Reproducible Scripts for the Publication

Drost HG, Bellstädt J, Ó'Maoiléidigh DS, Silva AT, Gabel A, Weinholdt C, Ryan PT, Dekkers BJW, Bentsink L, Hilhorst H, Ligterink W, Wellmer F, Grosse I, and Quint M. (2016). __Post-embryonic hourglass patterns mark ontogenetic transitions in plant development__. _Mol Biol Evol_ [doi:10.1093/molbev/msw039](http://mbe.oxfordjournals.org/content/early/2016/03/15/molbev.msw039).

Preprint: [Post-embryonic hourglass patterns mark ontogenetic transitions in plant development](http://biorxiv.org/content/early/2015/12/28/035527).

MBE Journal: [Post-embryonic hourglass patterns mark ontogenetic transitions in plant development](http://mbe.oxfordjournals.org/content/early/2016/03/15/molbev.msw039).

## Performing Phylostratigraphy

__A new version of the Phylostratigraphy algorithm implemented by Alexander Gabel can be found at:__ [https://github.com/AlexGa/Phylostratigraphy](https://github.com/AlexGa/Phylostratigraphy).

__We recommend to use the new version for performing custom phylostratigraphy due to easier applicability to any genome of interest. The older version implemented by Alexander is useful to repreduce the results shown in this study, but was not implemented to be applicable to any genome of interest.__ 


The first step in performing phylotranscriptomic analyses is to
perform `Phylostratigraphy` to assign
each protein coding genes of a query organism an evolutionary age. The resulting
`Phylostratigraphic Map` of a query organism is then matched
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
    
1) Make sure that [BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/2.2.21/) is installed on your machine.

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

4) Use the following command to start the Perl createPsMap.pl script (implemented (cc) by Alexander Gabel)
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


## Reading Datasets

The following script allows users to read [Supplementary Dataset 1](http://www.biorxiv.org/highwire/filestream/9526/field_highwire_adjunct_files/1/035527-2.xlsx).

```r
install.packages("readxl")
library(readxl)

### read PhyloExpressionSets
# PhyloExpressionSet: Germination
Ath.PhyloExpressionSet.Germination <- read_excel("035527-2.xlsx",sheet = 1)
# PhyloExpressionSet: Flowering
Ath.PhyloExpressionSet.Flowering <- read_excel("035527-2.xlsx",sheet = 2)
# PhyloExpressionSet: Flower Development
Ath.PhyloExpressionSet.FlowerDevelopment <- read_excel("035527-2.xlsx",sheet = 3)
```

### Users can find the raw datasets here:

- [Germination: GSE65394](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?token=onyxsyycjtaxxux&acc=GSE65394)
- [Flowering: PRJNA311774](http://www.ncbi.nlm.nih.gov/bioproject/311774)
- [Flower Development: GSE64581](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64581)

## Mapping Gene IDs

It is now assumed that the `Phylostratigraphic Map` of interest and the corresponding gene expression data set are joined. For this purpose the `MatchMap()` function implemented in the [myTAI](http://cran.r-project.org/web/packages/myTAI/index.html) package can be used. See `?myTAI::MatchMap` for details.


## Generating Figures

For matters of scientific reproducibility we developed the [myTAI](http://cran.r-project.org/web/packages/myTAI/index.html) R package to allow users to perform and reproduce all analyses presented in this study in an easy to use, well documented, and intuitive way. Please read the [Introduction to myTAI](https://github.com/HajkD/myTAI/blob/master/vignettes/Introduction.Rmd) to see the full potential of this software package.

First install and load the [myTAI](https://github.com/HajkD/myTAI) package:

```r
install.packages("myTAI")
library(myTAI)
```


### Figure 2

```r
# Visualize the Transcriptome Age Index of A. thaliana Germination
 PlotPattern(Ath.PhyloExpressionSet.Germination,
            TestStatistic = "FlatLineTest",
            permutations  = 10000,
            type          = "l",
            lwd           = 9,
            cex           = 2,
            cex.lab       = 1.2,
            cex.axis      = 1.1,
            xlab          = "Development [stages]",
            ylab          = "Transcriptome Age Index",
            shaded.area   = FALSE,
            p.value       = TRUE,
            y.ticks       = 5)


# The corresponding `Reductive Hourglass Test` statistically quantifies that
# the observed hourglass pattern does follow an high-low-high pattern of transcriptome age
ReductiveHourglassTest(Ath.PhyloExpressionSet.wholeSeed,
                       modules       = list(early = 1:2, mid = 3:5, late = 6:7),
                       plotHistogram = TRUE,
                       permutations  = 10000)

# Visualize the Relative Expression Profiles of Phylostrata during A. thaliana Germination
 PlotRE( ExpressionSet = Ath.PhyloExpressionSet.Germination,
         Groups        = list(1:3,4:12),
         legendName    = "PS",
         lty           = 1,
         cex           = 1.5,
         cex.lab       = 1.3,
         cex.axis      = 1.3,
         lwd           = 7, 
         xlab          = "Development [stages]")
```


### Figure 3


```r
# Visualize the Transcriptome Age Index of A. thaliana Flowering
PlotPattern(Ath.PhyloExpressionSet.Flowering,
            TestStatistic = "FlatLineTest",
            permutations  = 10000,
            type          = "l",
            lwd           = 9,
            cex           = 2,
            cex.lab       = 1.2,
            cex.axis      = 1.1,
            xlab          = "Development [stages]",
            ylab          = "Transcriptome Age Index",
            shaded.area   = FALSE,
            p.value       = TRUE,
            y.ticks       = 5)

# The corresponding `Reductive Hourglass Test` statistically quantifies that
# the observed hourglass pattern does follow an high-low-high pattern of transcriptome age
ReductiveHourglassTest(Ath.PhyloExpressionSet.Flowering,
                       modules       = list(early = 1:3, mid = 4:6, late = 7:9),
                       plotHistogram = TRUE,
                       permutations  = 10000)


# Visualize the Relative Expression Profiles of Phylostrata during A. thaliana Flowering
 PlotRE( ExpressionSet = Ath.PhyloExpressionSet.Flowering,
         Groups        = list(1:3,4:12),
         legendName    = "PS",
         lty           = 1,
         cex           = 1.5,
         cex.lab       = 1.3,
         cex.axis      = 1.3,
         lwd           = 7, 
         xlab          = "Development [stages]")
```


### Figure 4

```r
# Visualize the Transcriptome Age Index of A. thaliana Flower Development
PlotPattern(Ath.PhyloExpressionSet.FlowerDevelopment,
            TestStatistic = "FlatLineTest",
            permutations  = 10000,
            type          = "l",
            lwd           = 9,
            cex           = 2,
            cex.lab       = 1.2,
            cex.axis      = 1.1,
            xlab          = "Development [stages]",
            ylab          = "Transcriptome Age Index",
            shaded.area   = FALSE,
            p.value       = TRUE,
            y.ticks       = 5)

# Visualize the Relative Expression Profiles of Phylostrata during A. thaliana Flower Development
 PlotRE( ExpressionSet = Ath.PhyloExpressionSet.FlowerDevelopment,
         Groups        = list(1:3,4:12),
         legendName    = "PS",
         lty           = 1,
         cex           = 1.5,
         cex.lab       = 1.3,
         cex.axis      = 1.3,
         lwd           = 7, 
         xlab          = "Development [stages]")
```


## Supplementary Figures

### Suppl. Figure S3a

```r
par(mfrow = c(1,2))
# log2 transformed expression levels
PlotPattern(tf(Ath.PhyloExpressionSet.Germination,log2),
            TestStatistic = "FlatLineTest",
            permutations  = 10000,
            type          = "l",
            lwd           = 9,
            cex           = 2,
            cex.lab       = 1.5,
            cex.axis      = 1.5,
            main          = "log2 transformed",
            xlab          = "Development [stages]",
            ylab          = "Transcriptome Age Index",
            shaded.area   = FALSE,
            p.value       = TRUE)

# sqrt transformed expression levels
PlotPattern(tf(Ath.PhyloExpressionSet.Germination,sqrt),
            TestStatistic = "FlatLineTest",
            permutations  = 10000,
            type          = "l",
            lwd           = 9,
            cex           = 2,
            cex.lab       = 1.5,
            cex.axis      = 1.5,
            main          = "sqrt transformed",
            xlab          = "Development [stages]",
            ylab          = "Transcriptome Age Index",
            shaded.area   = FALSE,
            p.value       = TRUE)
```

### Suppl. Figure S3b

```r
par(mfrow = c(1,2))
# log2 transformed expression levels
PlotPattern(tf(Ath.PhyloExpressionSet.Flowering,log2),
            TestStatistic = "FlatLineTest",
            permutations  = 10000,
            type          = "l",
            lwd           = 9,
            cex           = 2,
            cex.lab       = 1.5,
            cex.axis      = 1.5,
            main          = "log2 transformed",
            xlab          = "Development [stages]",
            ylab          = "Transcriptome Age Index",
            shaded.area   = FALSE,
            p.value       = TRUE)

# sqrt transformed expression levels
PlotPattern(tf(Ath.PhyloExpressionSet.Flowering,sqrt),
            TestStatistic = "FlatLineTest",
            permutations  = 10000,
            type          = "l",
            lwd           = 9,
            cex           = 2,
            cex.lab       = 1.5,
            cex.axis      = 1.5,
            main          = "sqrt transformed",
            xlab          = "Development [stages]",
            ylab          = "Transcriptome Age Index",
            shaded.area   = FALSE,
            p.value       = TRUE)
```



### Suppl. Figure S5

```r
# reading bolting data
Bolting <- read.csv("BoltingData.csv",header = FALSE, sep = ";")
names(Bolting) <- c("nPlantsBolting","days")

# Visualize bolting
plot((Bolting[ , 1]/56) * 100,
     type     = "b", 
     lwd      = 7, 
     xaxt     = "n",
     xlab     = "days after shift to LD",
     ylab     = "% plants", 
     cex.lab  = 1.5, 
     cex.axis = 1.5)
     
axis(1,seq_len(length(Bolting[ , 1])), 
     labels   = Bolting[ , 2], 
     cex.lab  = 1.5, 
     cex.axis = 1.5) 
     
abline(h = 100, lty = 2, lwd = 2)
```
