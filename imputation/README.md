# Genetic imputation
This [AncestryMapper](https://cran.r-project.org/web/packages/AncestryMapper/vignettes/AncestryMapper2.0.html) package

## How-to
1) Install R for [Linux](https://cran.r-project.org/bin/linux/), [Mac OS](https://cran.r-project.org/bin/macosx/) or [Windows](https://cran.r-project.org/bin/windows/base/)
2) Get [plink](https://www.cog-genomics.org/plink/2.0/) from this website https://www.cog-genomics.org/plink/2.0/. *plink* will be used to convert e.g. VCF or 23andme format to tped format for *AncestryMaper*.
3) Set path to *plink* in `config.yml` file.
4) Replace reduced reference data set by the complete data
   - Download full reference data set from http://bit.ly/1OUstDP
   - Find the path to AncestryMapper package by typing `.libPaths()` in R terminal
   - Replace the reduced data `AncestryMaper/data/*.rda` by the downloaded reference data *(for testing purpose, this step can be skipped. However, the result accuracy may be affected!)*
   - Copy the provided `MinMaxFreq.rds` to the same folder above `AncestryMaper/data/`
5) Run `runAncestryMapper.py` script to calculate ancestry composition
```
python runAncestryMapper.py -i ../example.23andme.txt -t 23andme -o output
```
or
```
python runAncestryMapper.py -i ../example.vcf -t vcf -o output
```

## Input and Output
### Input
Input is SNP data in VCF or 23andme format (see [example files](https://github.com/trvinh/genomes-io-prj/tree/master/ancestry)).

### Output
The python script will give 3 output files:

- *output.txt* shows the distances of the input person to the **geography regions**.

- *output.amid* shows the distances of the individual to the **reference ethnic groups**.

- *output.pdf* is a plot representing the distances in the *output.amid* file.

Abbreviations using for the geo regions:
```
AFR Africa
AME America
AMER North, Central and South America
CSA Central South Asia
CSEA Central South East Asia
EUR Europa
ME Middle East
NEA North East Atlantic
OCE Ocenia
SEA South East asia
NaN Unknown
```

*Note: this method does not give the ancestry composition (e.g. 80% Asia, 19% Africa, 1% Moon). Its result is represented as **the distances between the person to the reference ethnic groups**. These distances are ranged between 0 and 1; but in practice, they are mostly between 0.6 and 0.9*
