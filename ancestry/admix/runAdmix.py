#!/usr/bin/env python
import getopt
import sys
import subprocess
import os
import yaml
import json

def load_config(config_file):
	with open(config_file, 'r') as stream:
		try:
			return yaml.safe_load(stream)
		except yaml.YAMLError as exc:
			print(exc)

def convert_to_json(inFile):
	outDict = { "model" : "composition" }
	key = ""
	value = ""
	with open(inFile) as fp:
		for line in fp:
			str = line.rstrip()
			if len(str) > 0:
				if str.find(':') < 0:
					key = str
				else:
					if len(key) > 0:
						if key in outDict:
							outDict[key] += "; " + str
						else:
							outDict[key] = str
	del outDict['model']
	return(outDict)

def subprocess_cmd(commands):
	for cmd in commands:
		subprocess.call(cmd, shell = True)

def check_vcf(vcfFile):
	fileExt = vcfFile.split(".")[-1]
	flag = 0
	if fileExt == 'gz':
		extractCmd = 'bgzip -dc %s > temp.vcf' % vcfFile
		subprocess.call([extractCmd], shell = True)
		vcfInput = 'temp.vcf'
		flag = 1
	else:
		vcfInput = vcfFile
	# check if vcf file already annotated
	cmd = 'cut -f 3 %s | sort -u | grep "^rs"' % vcfInput
	proc = subprocess.Popen([cmd], stdout = subprocess.PIPE, shell = True)
	(out, err) = proc.communicate()
	if len(out) == 0:
		if fileExt == 'gz':
			subprocess.call(['rm', 'temp.vcf'])
		sys.exit('Input file need to be annotated!')

def main(argv):
	input = ''
	type = '23andme'
	output = 'output'
	model = 'K47'
	try:
		opts, args = getopt.getopt(argv,"i:t:m:o:h",["in", "type", "model", "out", "help"])
	except getopt.GetoptError:
		print('runAdmix.py -i input.vcf -t vcf -m K47 -o output')
		sys.exit(2)

	for opt,arg in opts:
		if opt in ('-h','--help'):
			print('runAdmix.py -i <input file> -t <input type> -o <output file>')
			print('Input types (default = 23andme):')
			print('\tvcf (require plink for recoding!)')
			print('\t23andme')
			print('Model: please choose one of the following options (default = K47)')
			print('\tall, K7b, K12b, globe13, globe10, world9, Eurasia7, Africa9, weac2, E11, K36,')
			print('\tEUtest13, Jtest14, HarappaWorld, TurkicK11, KurdishK10, AncientNearEast13, K7AMI,')
			print('\tK8AMI, MDLPK27, puntDNAL, K47, K7M1, K13M2, K14M1, K18M4, K25R1, MichalK25')
			sys.exit()
		elif opt in ('-i','--in'):
			input = arg
		elif opt in ('-t','--type'):
			type = arg
		elif opt in ('-o','--out'):
			output = arg
		elif opt in ('-m','--model'):
			model = arg

	cfg = load_config('config.yml')

	allModels = ["all", "K7b", "K12b", "globe13", "globe10", "world9", "Eurasia7", "Africa9", "weac2", "E11", "K36", "EUtest13", "Jtest14", "HarappaWorld", "TurkicK11", "KurdishK10", "AncientNearEast13", "K7AMI", "K8AMI", "MDLPK27", "puntDNAL", "K47", "K7M1", "K13M2", "K14M1", "K18M4", "K25R1", "MichalK25"]
	if model not in allModels:
		print('Model: Invalid model given. Please choose one of the following options')
		print('\tall, K7b, K12b, globe13, globe10, world9, Eurasia7, Africa9, weac2, E11, K36,')
		print('\tEUtest13, Jtest14, HarappaWorld, TurkicK11, KurdishK10, AncientNearEast13, K7AMI,')
		print('\tK8AMI, MDLPK27, puntDNAL, K47, K7M1, K13M2, K14M1, K18M4, K25R1, MichalK25')
		exit()

	if not os.path.exists('tmp'):
		os.makedirs('tmp')
	if type == "vcf":
		print('Check VCF input...')
		check_vcf(input)
		print('Convert VCF to 23andme...')
		fileExt = input.split(".")[-1]
		if fileExt == 'gz':
			removeContig = 'bgzip -dc %s | sed "s/^chr//g" | grep "^[#,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,(MT)]" > tmp/tmp.vcf' % (input)
		else:
			removeContig = 'sed "s/^chr//g" %s | grep "^[#,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,(MT)]" %s > tmp/tmp.vcf' % (input)
		plinkCMD = '%s --vcf tmp/tmp.vcf --snps-only --recode 23 --out tmp/tmp.23andme' % (cfg['plink'])
		removeVCF = 'rm tmp/tmp.vcf'
		subprocess_cmd((removeContig, plinkCMD, removeVCF))
	elif type == "23andme":
		print('Copy input to tmp folder...')
		cpCMD = 'cp %s tmp/tmp.23andme.txt' % (input)
		subprocess.call(cpCMD, shell = True)
	else:
		print('Invalid input types given. Accepted input types:')
		print('\tvcf')
		print('\t23andme')
		exit()

	# remove unidentified SNPs
	removeSNP = 'awk \'$2 != \".\"\' tmp/tmp.23andme.txt > tmp/tmp.23andme.txt.mod'
	replaceCMD = 'mv tmp/tmp.23andme.txt.mod tmp/tmp.23andme.txt'
	subprocess_cmd((removeSNP, replaceCMD))

	# run admix
	if model == "all":
		admixCMD = 'admix -f %s -v 23andme -o %s' % ("tmp/tmp.23andme.txt", output)
	else:
		admixCMD = 'admix -f %s -v 23andme -m %s -o %s' % ("tmp/tmp.23andme.txt", model, output)
	subprocess.call(admixCMD, shell = True)

	# print output
	jsonOut = convert_to_json(output)
	jsonFile = output + ".json"
	with open(jsonFile, 'w') as fp:
		json.dump(jsonOut, fp)
	mvCMD = 'mv %s tmp' % (output)
	rmCMD = 'rm -rf tmp'
	subprocess_cmd((mvCMD, rmCMD))
	print('Finished! Check output %s.json' % (output))

if __name__ == "__main__":
	if len(sys.argv[1:]) == 0:
		print('Missing input!')
		print('runAdmix.py -i input.vcf -t vcf -m K74 -o output.txt')
		sys.exit(2)
	else:
		main(sys.argv[1:])
