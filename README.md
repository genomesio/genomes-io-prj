# Genetic testing

Modules to perform genetic testing for human genomes or [SNPs data](https://github.com/trvinh/genomes-io-prj/tree/master/example_input).
## Ancestry testing
Ancestry test is done based on 2 approaches: [admix](https://github.com/stevenliuyi/admix) and [AncestryMapper](https://cran.r-project.org/web/packages/AncestryMapper/vignettes/AncestryMapper2.0.html).

**admix** calculates the `ancestry composition` (e.g. 48.08% North European, 17.77% Mediterranean, 10.21% East African, etc.) for input genome based on different public [admixture models](https://en.wikipedia.org/wiki/Genetic_admixture). The algorithm is based on the maximum likelihood estimation.

**AncestryMapper** estimates the similarity of the input genome to several reference ethnic groups by computing the `genetic distances` based on the [Principal component analysis (PCA) procedure](https://en.wikipedia.org/wiki/Principal_component_analysis). This apporach will not give the ancestry composition as the output like *admix*, but a sorted list from the most similar ethnic group to the most distance one together with the calculated genetic distances.

## Genetic testing based on imputation

([Definition from wiki](https://en.wikipedia.org/wiki/Imputation_(genetics)): Imputation in genetics refers to the statistic reference of unobserved genotypes. It is achieved by using known haplotypes in a population, for instance from the HapMap or the 1000 Genomes Project in humans, thereby allowing to test for association between a trait of interest (e.g. a disease) and experimentally untyped genetic variants, but whose genotypes have been statistically inferred ("imputed"). Genotype imputation is usually performed on SNP, the most common kind of genetic variation.)

This imputation module is adapted and modified from the source code of [www.impute.me](https://github.com/lassefolkersen/impute-me). The imputation is carried out using the tool [impute2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html). After that, genetic tests are performed by calculating the genetic risk scores for the individual with different traits including several complex diseases, autoimmune diseases, breast cancer, drug response, precision medicine (personal relevant diseases), rare diseases, public traits from the UK biobank, as well as the ethnicity, genetic height and hair color, and traits related to intelligence.
