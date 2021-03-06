#!/usr/bin/env python
import getopt
import sys
import subprocess
import os
import yaml

def load_config(config_file):
	with open(config_file, 'r') as stream:
		try:
			return yaml.safe_load(stream)
		except yaml.YAMLError as exc:
			print(exc)
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
	type = ''
	output = 'output'
	try:
		opts, args = getopt.getopt(argv,"i:t:o:h",["in", "type", "out", "help"])
	except getopt.GetoptError:
		print('runAncestryMapper.py -i input.vcf -t vcf -o output.json')
		sys.exit(2)

	for opt,arg in opts:
		if opt in ('-h','--help'):
			print('runAncestryMapper.py -i <input file> -t <input type> -o <name of output files>')
			print('Input types: ')
			print('\tvcf (require plink for recoding!)')
			print('\t23andme')
			sys.exit()
		elif opt in ('-i','--in'):
			input = arg
		elif opt in ('-t','--type'):
			type = arg
		elif opt in ('-o','--out'):
			output = arg

	cfg = load_config('config.yml')
	if not os.path.exists('tmp'):
		os.makedirs('tmp')

	if type == "vcf":
		print('Check VCF input...')
		check_vcf(input)
		print('Convert VCF to tped...')
		fileExt = input.split(".")[-1]
        if fileExt == 'gz':
			removeContig = 'bgzip -dc %s | sed "s/^chr//g" | grep "^[#,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,(MT)]" > tmp/tmp.vcf' % (input)
		else:
			removeContig = 'sed "s/^chr//g" %s | grep "^[#,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,(MT)]" %s > tmp/tmp.vcf' % (input)
		plinkCMD = '%s --vcf tmp/tmp.vcf --recode transpose --out tmp/input' % (cfg['plink'])
		removeVCF = 'rm tmp/tmp.vcf'
		subprocess_cmd((removeContig, plinkCMD, removeVCF))
	elif type == "23andme":
		print('Convert 23andme to tped...')
		plinkCMD = '%s --23file %s --recode transpose --out tmp/input' % (cfg['plink'], input)
		subprocess.run(plinkCMD, shell = True)
	else:
		print('Accepted input types: ')
		print('\tvcf (require plink for recoding!)')
		print('\t23andme')
		sys.exit()

	# remove NA SNP IDs and create list of snps
	removeSNP = 'awk \'$2 != \".\"\' tmp/input.tped > tmp/input.tpedMod'
	replaceTped = 'mv tmp/input.tpedMod tmp/input.tped'
	snpList = 'awk -F \' \' \'{print $2}\' tmp/input.tped > tmp/input.snplist'
	subprocess_cmd((removeSNP, replaceTped, snpList))

	# ancestry mapper
	mapperCMD = 'Rscript mapper.R tmp %s' % (output)
	removeTMP = 'rm -rf tmp'
	removeReport = 'rm report_*'
	renameAmid = 'mv AMid*.amid %s.amid' % (output)
	subprocess_cmd((mapperCMD, removeTMP, removeReport, renameAmid))

	print('Finished! Check output %s.json and %s.amid' % (output, output))

if __name__ == "__main__":
	if len(sys.argv[1:]) == 0:
		print('Missing input!')
		print('runAncestryMapper.py -i input.vcf -t vcf -o output')
		sys.exit(2)
	else:
		main(sys.argv[1:])
