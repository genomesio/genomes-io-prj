# Ancestry mapping using [AncestryMapper]() package

1) Install R for [Linux](https://cran.r-project.org/bin/linux/), [Mac OS](https://cran.r-project.org/bin/macosx/) or [Windows](https://cran.r-project.org/bin/windows/base/)
2) Get [plink](https://www.cog-genomics.org/plink/2.0/) from this website https://www.cog-genomics.org/plink/2.0/. *plink* will be used to convert e.g. VCF or 23andme format to tped format for *AncestryMaper*.
3) Set path to *plink* in `config.yml` file.
4) Replace reduced reference data set by the complete data
   - Download full reference data set from http://bit.ly/1OUstDP
   - Find the path to AncestryMapper package by typing `.libPaths()` in R terminal
   - Replace the reduced data `AncestryMaper/data/*.rda` by the downloaded reference data (for testing purpose, this step can be skipped. However, the result accuracy can be affected!)
   - Copy the provided `MinMaxFreq.rds` to the same folder above `AncestryMaper/data/`
5) Run `runAncestryMapper.py` script to calculate ancestry composition
```
python runAncestryMapper.py -i ../example.23andme.txt -t 23andme -o output
```
or
```
python runAncestryMapper.py -i ../example.vcf -t vcf -o output
```
