# Ancestry composition calculation using Admix

## HowTo
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

## Input and Output
### Input
Input is SNP data in VCF or 23andme format (see [example files](https://github.com/trvinh/genomes-io-prj/blob/master/ancestry/example.23andme.txt)).

### Output
Output is a text file showing the ancestry composition predicted using all available reference data (admixture models) in *Admix*. For general usage, one should use the results based on Dodecad's global models (global13, global10 or world9) or HarappaWorld.

For example:
```
globe13
Siberian: 5.63%
Amerindian: 7.21%
West African: 0.00%
Palaeo African: 0.63%
Southwest Asian: 5.19%
East Asian: 0.00%
Mediterranean: 17.77%
Australasian: 0.00%
Artic: 0.97%
West Asian: 4.31%
North European: 48.08%
South Asian: 0.00%
East African: 10.21%

globe10
Ameriandian: 8.92%
West Asian: 4.41%
Australasian: 0.00%
Palaeo African: 2.30%
Neo African: 6.41%
Siberian: 7.09%
Southern: 11.11%
East Asian: 0.00%
Atlantic Baltic: 59.46%
South Asian: 0.31%

world9
Amerindian: 7.83%
East Asian: 0.00%
African: 8.04%
Atlantic Baltic: 55.40%
Australasian: 0.00%
Siberian: 8.49%
Caucasus Gedrosia: 11.05%
Southern: 9.20%
South Asian: 0.00%

HarappaWorld
South-Indian: 0.00%
Baloch: 14.43%
Caucasian: 0.00%
Northeast-Euro: 48.78%
Southeast-Asian: 1.33%
Siberian: 5.64%
Northeast-Asian: 0.00%
Papuan: 0.00%
American: 3.89%
Beringian: 1.92%
Mediterranean: 13.49%
Southwest-Asian: 0.53%
San: 0.00%
East-African: 8.21%
Pygmy: 1.78%
West-African: 0.00%
```

*Please be aware about the precision by using each model as well as the **Term of use** for the models!*
