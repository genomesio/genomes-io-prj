# Genetic imputation based on impute.me
This imputation module is created based on the source code of [impute.me](https://www.impute.me/). The original souce code can be found at https://github.com/lassefolkersen/impute-me

## How-to
1. Install R for [Linux](https://cran.r-project.org/bin/linux/), [Mac OS](https://cran.r-project.org/bin/macosx/) or [Windows](https://cran.r-project.org/bin/windows/base/)
2. Clone this repository into your computer.
3. Run `setup.sh` script to download reference data from 1000genomes and other tools for the imputation ([plink1.9](https://www.cog-genomics.org/plink/1.9/), [impute2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html), [shapeit](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html) and [gtool](https://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool.html))
4. Check paths to the reference data, tools and imputation modules in the `config.yml` file.
5. Run `runImputation.py` script to do the genetic imputation and tests
```
python runImputation.py -i ../example.23andme.txt -t 23andme -o impute_output [-n job_id]
```

## Input and Output
### Input
Input is SNP data in VCF or 23andme format (see [example files](https://github.com/trvinh/genomes-io-prj/tree/master/ancestry)). *Note: the example.vcf will not work with this imputation module!!! It has SNPs for chromosome 1 only.*

VCF inputs need to be annotated with [rsIDs](https://www.ncbi.nlm.nih.gov/books/NBK44417/#Content.what_is_a_reference_snp_or__rs_i). You can annotate your genome using [this instruction](https://gist.github.com/trvinh/43a0e0724ee7330a45d0c7074f1c0e5f) if needed.

### Output
The `runImputation.py` will generate a random job ID each time it's called if `-n/--id` not given.

The **imputation result** will stored in the `imputation/job_id` folder.

Based on the output in the `imputation` folder, several **genetic tests** will be performed. Each test has its own reference SNP data and script saved in the corresponding subfolder within the `imputeTraits` folder. **Genetic testing results** are saved separately in the `output/job_id` folder, where all output in json format can be found in `output/job_id/json_out/`. See [example](https://github.com/trvinh/genomes-io-prj/blob/master/imputation/example/).

## Genetic tests

Currently there are several available tests

| No. | Test | Output | Note |
| --- | ---- | ------ | ---- |
| 1 | AllDiseases | GRS, trait, percentage, message | *GRS* = Z-score, *percentage* shows the percentage of the population that have a slower risk score than this person |
| 2 | autoimmuneDiseases | GRS, disease |   |
| 3 | BRCA | gene , user phenotype, normal phenotype, [clinVar](https://www.ncbi.nlm.nih.gov/clinvar/intro/), polyphen, sift, consequence type | gene is either BRCA1 or BRCA2, which produce a hereditary breast-ovarian cancer syndrome. Read more about [polyphen](http://genetics.bwh.harvard.edu/pph2/), [sift](https://www.ncbi.nlm.nih.gov/pubmed/19561590) |
| 4 | drugResponse | GRS, disease, drug , percentage |   |
| 5 | ethnicity | guessed super population, PCA coordinates, SNP count | a PCA plot can be generated using `pcaPlot.R` script in the `impute-me/ethnicity` folder. Usage: `Rscript pcaPlot.R job_id_data.json impute-me/ethnicity/2017-04-03_ethnicity_pca.rdata impute-me/ethnicity/2017-04-03_ethnicity_descriptions.txt` |
| 6 | intelligence | GRS, trait, percentage, message | same as `AllDiseases`  |
| 7 | precisionMedicine | GRS, disease, drug, percentage |   |
| 8 | rareDiseases | message, diseases of interest, all comparison (all_findings) (*) | (*) if user has a SNP that is present in the reference data, the entry in `all_findings` will have 5 sub-entries, including SNP ID, user genotype, risk allele, non-risk-allele and the corresponding disease. If user's genotype has risk allele, that disease will be mentioned in the `advice` and `diseases of interest` |
| 9 | ukbiobank | GRS, trait, percentage, message | same as `AllDiseases` |
| 10 | athletics |athletics SNPs, injuries and dietary |   |
| 11 | appearance | guessed height (gheight) Z score, gheight SNP count, gheight estimate, hair colour |   |
| 12 | leukemia | GRS, mean and SD for case, mean and SD for control, message and 2 plots for CLL and ALL |   |
| 13 | microbiome | SNP, user genotype, increasing allele, bacteria |   |
| 14 | nonsenser | SNP, user genotype, common allele, freq, SIFT, Mutation taster, polyphen2, type |   |


"A polygenic risk score is a value that gives a summary of a large number of different SNPs - each of which contribute a little to disease risk. The higher the value, the higher the risk of developing disease. Of course the interpretation of this risk depends a lot on other factors as well: How heritable the disease is. How much of this heritability we can explain with known SNPs. And not least, what would the risk of disease be for you otherwise, i.e. without taking the genetic component into account.

Because the polygenic risk score is only a risk-modifier, knowledge of these three other values are all required if you want to know what your overall risk is, i.e. what's the chance in percent. This calculator cannot provide that. But it can provide a view of the known genetic component of your disease risk, based on all the SNPs that we know are associated with the disease. This, we believe, makes it a better choice for complex diseases than the typical one-SNP-at-the time analysis typically seen in consumer genetics." (copied from impute.me)

The polygenic risk score is represented as the [Z-score](https://en.wikipedia.org/wiki/Standard_score), which is calculated by the mean and standard deviation of the population. If the score is lower than 20% of the population, it will be a *low score*. If it is greater than 90% of the population, it will be a *high score*. Otherwise, it will be a *fairly average score*.

An example of a *message* for the asthma test: *Ethnicity-corrected trait Z-score is -0.24. This genetic risk score is higher than 40% of the general population. This is a fairly average score. Result from the analysis of 59 SNPs from http://www.ncbi.nlm.nih.gov/pubmed/29273806 Demenais F et al (PMID 29273806), which were reported to be associated with asthma. This study reports a total sample size of 132486, as entered on date 2017-12-22.*

*__Note__: The numbers in the IDs for the traits (1), (4), (6), and (7) are the NCBI pubmed IDs (for example: acute_insulin_response_28490609 will be referenced to https://www.ncbi.nlm.nih.gov/pubmed/28490609).*
