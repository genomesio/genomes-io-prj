#!/usr/bin/env python
import os
import sys
import argparse
import yaml
import subprocess
import random
import string
from pathlib import Path

def randomStringDigits(stringLength = 6):
    lettersAndDigits = string.ascii_letters + string.digits
    return ''.join(random.choice(lettersAndDigits) for i in range(stringLength))

def checkFileExist(file):
	try:
		my_abs_path = Path(file).resolve(strict=True)
	except FileNotFoundError:
		sys.exit("%s not found" % file)

def load_config(config_file):
    with open(config_file, 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

def subprocess_cmd(commands):
	for cmd in commands:
		subprocess.call(cmd, shell = True)

def process_input(args):
    (infile, jobID, type, plink) = args
    rawInput = jobID + '_raw_data'
    if type == 'vcf':
        print('Convert VCF to 23andme...')
        removeContig = 'grep "^[#,1:22,X,Y,(MT)]" %s > imputation/tmp.vcf' % (infile)
        plinkCMD = '%s --vcf imputation/tmp.vcf --snps-only --recode 23 --out imputation/%s' % (plink, rawInput)
        removeVCF = 'rm imputation/tmp.vcf'
        subprocess_cmd((removeContig, plinkCMD, removeVCF))
    elif type == '23andme':
        cpCMD = 'cp %s imputation/%s\.txt' % (infile, rawInput)
        subprocess.run(cpCMD, shell = True)

def prepare_dir(args):
    (infile, mode, jobID, type, plink) = args
    if mode == 'full' or mode == 'impute':
        try:
            my_abs_path = Path('imputation').resolve(strict=True)
        except FileNotFoundError:
            Path('imputation').mkdir(parents = True, exist_ok = True)
        # [convert format and] copy infile to imputation folder
        process_input([infile, jobID, type, plink])
        # prepare folders
        runDir = os.getcwd() + '/imputation'
        outputDir = os.getcwd() + '/output/' + jobID
        Path(outputDir).mkdir(parents = True, exist_ok = True)
    elif mode == 'test':
        # get existing job ID
        outputDirTmp = os.getcwd() + '/output'
        checkFileExist(outputDirTmp)
        jobIDtmp = ''
        subFolder = os.listdir(outputDirTmp)
        if len(subFolder) < 1:
            print('No output folder found! Please check again or run with full mode!')
            sys.exit()
        elif len(subFolder) > 1:
            print('More than one output folders found:')
            for i in range(len(subFolder)):
                print (i + 1, end = " ")
                print (subFolder[i])
            jobIDtmp = input('Select one folder: ')
            checkFileExist(outputDirTmp + "/" + jobIDtmp)
        else:
            jobIDtmp = subFolder[0]
        runDir = os.getcwd() + '/imputation'
        outputDir = os.getcwd() + '/output/' + jobIDtmp
    return(runDir, outputDir)

def main():
    version = "1.0.1"
    parser = argparse.ArgumentParser(description="You are running runImputation.py version " + str(version) + ".")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('additional arguments')
    required.add_argument('-i', '--infile', help='Input file in vcf or 23andme format', action='store', default='', required=True)
    required.add_argument('-t', '--type', help='Type of input (vcf|23andme)', action='store', choices=['vcf', '23andme'], default='', required=True)
    optional.add_argument('-m', '--mode', help='Run mode (full|impute|test). Full will run both imputation and perform genetic testing.'
                                        'Test mode will run only the tests based on pre-computed imputation data. Default: full',
                                        choices=['full', 'impute', 'test'], action='store', default='full')
    optional.add_argument('-n', '--id', help='Job ID. If not given, a random string will be generated.', action='store', default='')
    optional.add_argument('--trail', help='Trail names for report. If not set, all available trails will be reported.', action='store', default='all')
    args = parser.parse_args()

    # get arguments and config paths
    cfg = load_config('config.yml')
    type = args.type
    mode = args.mode
    infile = args.infile
    checkFileExist(infile)
    jobID = args.id
    if (jobID == ''):
        jobID = 'id_' + randomStringDigits(9)
    trail = args.trail
    if not trail == "all":
        checkFileExist(cfg['imputeTrails'] + "/" + trail)

	# do imputation
    (runDir, outputDir) = prepare_dir([infile, mode, jobID, type, cfg['plink']])
    imputeCMD = 'Rscript imputation_a2z.R %s %s %s %s %s %s %s %s %s %s %s %s' % (trail, mode, jobID, runDir, outputDir, cfg['shapeit'], cfg['plink'], cfg['gtool'], cfg['impute2'], cfg['sample_ref'], cfg['imputeDataDir'], cfg['imputeTrails'])
    # print(imputeCMD)
    subprocess.run(imputeCMD, shell = True)
    print('Finished! Check outputs in %s folder.' % (outputDir))

if __name__ == '__main__':
    main()
