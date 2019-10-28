# Example output
Results of AncestryMapper using this [example input](https://github.com/trvinh/genomes-io-prj/blob/master/ancestry/example.23andme.txt)

* testvn.amid: a matrix contains the normalised euclidean distances (columns started with `I_`) and the the crude distances (columns started with `C_`)
* testvn.pdf: a plot showing how similar the user genome in comparison to the reference ethnic groups (populations). Data used for this plot are the normalised euclidean distances
* testvn_ethnic_dist.txt: the crude distances of the user genome to the reference ethnic groups. This table contains 3 columns including the name of the ethnic group, the corresponding distance and source of the reference data for that ethnic group.
```
ethnic_group	distance	source
Kinh	0.619509500078284	1KG
Han.South.China	0.619698050810253	1KG
Han.Dallas	0.622679660310717	HapMap
Dyak	0.623454868777498	Reich2011
Han	0.625168614122902	HGDP
Dai	0.625507551578514	1KG
...
```
* testvn_region_dist.txt: a table shows the distances of the user genome to the geographical regions. The distance to a specific region is calculated as the mean distance of the user genome to all ethnic groups within that region. The abbreviations of the regions are shown in the following table
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

The command to generate these files:

```
python runAncestryMapper.py -i example.23andme.txt -t 23andme -o testvn
```

