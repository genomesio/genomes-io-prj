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

def main(argv):
	input = ''
	type = '23andme'
	output = 'output.txt'
	try:
		opts, args = getopt.getopt(argv,"i:t:o:h",["in", "type", "out", "help"])
	except getopt.GetoptError:
		print('runAdmix.py -i input.vcf -t vcf -o output.txt')
		sys.exit(2)

	for opt,arg in opts:
		if opt in ('-h','--help'):
			print('runAdmix.py -i <input file> -t <input type> -o <output file>')
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

	if type == "vcf":
		print('Convert VCF to 23andme...')
		if not os.path.exists('tmp'):
			os.makedirs('tmp')
		plinkCMD = '%s --vcf %s --snps-only --recode 23 --out tmp/tmp.23andme' % (cfg['plink'], input)
		subprocess.run(plinkCMD, shell = True)
		admixCMD = 'admix -f %s -v 23andme -o %s' % ("tmp/tmp.23andme.txt", output)
	else:
		admixCMD = 'admix -f %s -v 23andme -o %s' % (input, output)

	subprocess.call(admixCMD, shell = True)
	print('Finished! Check output %s.txt' % (output))

if __name__ == "__main__":
	if len(sys.argv[1:]) == 0:
		print('Missing input!')
		print('runAdmix.py -i input.vcf -t vcf -o output.txt')
		sys.exit(2)
	else:
		main(sys.argv[1:])
