# Ancestry composition calculation using Admix

1) Install [Admix](https://github.com/stevenliuyi/admix)
```
pip install git+https://github.com/stevenliuyi/admix
```

2) Get [plink](https://www.cog-genomics.org/plink/2.0/) from this website https://www.cog-genomics.org/plink/2.0/. *plink* will be used to convert e.g. VCF format into 23andme format for *admix*.

3) Set path to *plink* in `config.yml` file.
4) Run `runAdmix.py` script to calculate ancestry composition
```
python runAdmix.py -i example.vcf -t vcf -o output.txt
```
or
```
python runAdmix.py -i example.23andme.txt -t 23andme -o output.txt
```
