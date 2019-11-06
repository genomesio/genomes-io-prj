# Ancestry composition calculation using Admix

## HowTo
1) Install [Admix](https://github.com/stevenliuyi/admix)
```
pip install git+https://github.com/stevenliuyi/admix
```

2) Get [plink](https://www.cog-genomics.org/plink/1.9/) from this website https://www.cog-genomics.org/plink/1.9/. *plink* will be used to convert e.g. VCF format into 23andme format for *admix*.

3) Set path to *plink* in `config.yml` file.
4) Run `runAdmix.py` script to calculate ancestry composition
```
python runAdmix.py -i example.vcf -t vcf -m K47 -o output
```
or
```
python runAdmix.py -i example.23andme.txt -t 23andme -m K47 -o output
```
or
```
python runAdmix.py --help
```
for parameter description :)

## Input and Output
### Input
Input is SNP data in VCF or 23andme format (see [example files](https://github.com/trvinh/genomes-io-prj/tree/master/ancestry)).

### Output
Output is a JSON file showing the ancestry composition predicted using all available reference data (admixture models) in *Admix*. For general usage, one should use the results based on Dodecad's global models (global13, global10 or world9) or HarappaWorld.

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

## Models used in admix tool
([source](https://github.com/stevenliuyi/admix/blob/master/README.md#models))
Admix supports many publicly available admixture models. All the calculator files are properties of their authors, and are not covered by the license of this program. Links are provided which contain more information for each model.

| model value | model name | source |
| ----- | --------- | ---- |
| `K7b` | Dodecad K7b | [Link](http://dodecad.blogspot.com/2012/01/k12b-and-k7b-calculators.html) |
| `K12b` | Dodecad K12b | [Link](http://dodecad.blogspot.com/2012/01/k12b-and-k7b-calculators.html) |
| `globe13` | Dodecad globe13 | [Link](http://dodecad.blogspot.com/2012/10/globe13-calculator.html) |
| `goble10` | Dodecad globe10 | [Link](http://dodecad.blogspot.com/2012/10/globe10-calculator.html) |
| `world9` | Dodecad world9 | [Link](http://dodecad.blogspot.com/2011/12/world9-calculator.html) |
| `Eurasia7` | Dodecad Eurasia7 | [Link](http://dodecad.blogspot.com/2011/10/eurasia7-calculator.html) |
| `Africa9` | Dodecad Africa9 | [Link](http://dodecad.blogspot.com/2011/09/africa9-calculator.html) |
| `weac2` | Dodecad weac (West Eurasian cline) 2 | [Link](http://dodecad.blogspot.com/2012/06/weac2-calculator.html) |
| `E11` | E11 | [Link](http://www.ranhaer.com/thread-32241-1-1.html) |
| `K36` | Eurogenes K36 | [Link](http://bga101.blogspot.com/2013/03/eurogenes-k36-at-gedmatch.html) |
| `EUtest13` | Eurogenes EUtest K13 | [Link](http://bga101.blogspot.com/2013/11/updated-eurogenes-k13-at-gedmatch.html) |
| `Jtest14` | Eurogenes Jtest K14 | [Link](http://bga101.blogspot.com/2012/09/eurogenes-ashkenazim-ancestry-test-files.html) |
| `HarappaWorld` | HarappaWorld | [Link](http://www.harappadna.org/2012/05/harappaworld-admixture/) |
| `TurkicK11` | Turkic K11 | [Link](http://www.anthrogenica.com/showthread.php?8817-Turkic-K11-Admixture-Calculator) |
| `KurdishK10` | Kurdish K10 | [Link](http://www.anthrogenica.com/showthread.php?8571-K10-Kurdish-Calculator-Version-1/page6) |
| `AncientNearEast13` | Ancient Near East K13 | [Link](http://www.anthrogenica.com/showthread.php?8193-ancient-DNA-in-the-Gedrosia-Near-East-Neolithic-K13) |
| `K7AMI` | Eurogenes K7 AMI | [Link](http://www.anthrogenica.com/showthread.php?4548-Upcoming-DIY-Eurogenes-K7-amp-K8-Calculator-amp-Oracles-for-tracking-E-Asian-amp-ASI) |
| `K8AMI` | Eurogenes K8 AMI | [Link](http://www.anthrogenica.com/showthread.php?4548-Upcoming-DIY-Eurogenes-K7-amp-K8-Calculator-amp-Oracles-for-tracking-E-Asian-amp-ASI) |
| `MDLPK27` | MDLP K27 | [Link](http://www.anthrogenica.com/showthread.php?4557-Post-MDLP-K27-Results) |
| `puntDNAL` | puntDNAL K12 Ancient World | [Link](http://www.anthrogenica.com/showthread.php?8034-PuntDNAL-K12-Ancient-World-Results) |
| `K47` | LM Genetics K47 | [Link](https://anthrogenica.com/showthread.php?12788-New-K30-K47-World-Calculator) |
| `K7M1` | Tolan K7M1 | [Link](http://gen3553.pagesperso-orange.fr/ADN/Calc.htm) |
| `K13M2` | Tolan K13M2 | [Link](http://gen3553.pagesperso-orange.fr/ADN/Calc.htm) |
| `K14M1` | Tolan K14M1 | [Link](http://gen3553.pagesperso-orange.fr/ADN/Calc.htm) |
| `K18M4` | Tolan K18M4 | [Link](http://gen3553.pagesperso-orange.fr/ADN/Calc.htm) |
| `K25R1` | Tolan K25R1 | [Link](http://gen3553.pagesperso-orange.fr/ADN/Calc.htm) |
| `MichalK25`| Michal World K25 | [Link](https://anthrogenica.com/showthread.php?13359-Michal-s-World-K25-calculator) |

*Please be aware about the precision by using each model as well as the **Term of use** for the models!*
