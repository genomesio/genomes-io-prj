#!/usr/bin/env python
import getopt
import sys
import subprocess
import os
import yaml
import random
import string

def randomStringDigits(stringLength = 6):
    lettersAndDigits = string.ascii_letters + string.digits
    return ''.join(random.choice(lettersAndDigits) for i in range(stringLength))

def load_config(config_file):
    with open(config_file, 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

def subprocess_cmd(commands):
	for cmd in commands:
		subprocess.call(cmd, shell = True)

def main(argv):
	input = ''
	type = ''
	mode = 'full'
	try:
		opts, args = getopt.getopt(argv,"i:t:m:h",["in", "type", "mode", "help"])
	except getopt.GetoptError:
		print('runImputation.py -i input.txt -t 23andme -m full')
		sys.exit(2)

	for opt,arg in opts:
		if opt in ('-h','--help'):
			print('runImputation.py -i <input file> -t <input type> -m <run mode>')
			print('Input types: ')
			print('\tvcf (require plink for recoding!)')
			print('\t23andme')
			print('Run mode:')
			print('\tfull (imputation and genetic testing)')
			print('\ttest (only genetic testing)')
			sys.exit()
		elif opt in ('-i','--in'):
			input = arg
		elif opt in ('-t','--type'):
			type = arg
		elif opt in ('-m','--mode'):
			mode = arg

	cfg = load_config('config.yml')

	if mode == 'full':
		# create random string used as ID
		uniqueID = 'id_' + randomStringDigits(9)

		# [convert format and] copy input to imputation folder
		if not os.path.exists('imputation'):
			os.makedirs('imputation')

		rawInput = uniqueID + '_raw_data'

		if type == 'vcf':
			print('Convert VCF to 23andme...')
			plinkCMD = '%s --vcf %s --snps-only --recode 23 --out imputation/%s' % (cfg['plink'], input, rawInput)
			subprocess.run(plinkCMD, shell = True)
		elif type == '23andme':
			cpCMD = 'cp %s imputation/%s\.txt' % (input, rawInput)
			subprocess.run(cpCMD, shell = True)
		else:
			print('Input types not given or invalid:')
			print('\tvcf (require plink for recoding!)')
			print('\t23andme')
			sys.exit()

		# prepare folders
		runDir = os.getcwd() + '/imputation'
		if not os.path.exists('output'):
			os.makedirs('output')
		outputDir = os.getcwd() + '/output/' + uniqueID
		if not os.path.exists(outputDir):
			os.makedirs(outputDir)
	elif mode == 'test':
		# get existing job ID
		outputDirTmp = os.getcwd() + '/output'
		uniqueID = ''
		subFolder = os.listdir(outputDirTmp)
		if len(subFolder) < 1:
			print('No output folder found! Please check again or run with full mode!')
			sys.exit()
		elif len(subFolder) > 1:
			print('More than one output folders found! Please check again!')
			sys.exit()
		else:
			uniqueID = subFolder[0]
		runDir = os.getcwd() + '/imputation'
		outputDir = os.getcwd() + '/output/' + uniqueID
	else:
		print('Run mode not given or invalid:')
		print('\tfull (imputation and genetic testing)')
		print('\ttest (only genetic testing)')
		sys.exit()

	# do imputation and testing
	imputeCMD = 'Rscript imputation_a2z.R %s %s %s %s %s %s %s %s %s %s %s' % (mode, uniqueID, runDir, outputDir, cfg['shapeit'], cfg['plink'], cfg['gtool'], cfg['impute2'], cfg['sample_ref'], cfg['imputeDataDir'], cfg['impute_me'])
	subprocess.run(imputeCMD, shell = True)
	print('Finished! Check outputs in %s folder.' % (outputDir))

if __name__ == '__main__':
	if len(sys.argv[1:]) == 0:
		print('Missing input!')
		print('runImputation.py -i input.txt -t 23andme -m full')
		sys.exit(2)
	else:
		main(sys.argv[1:])
