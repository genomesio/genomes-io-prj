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
    # check complete annotated chromosomes
    cmd = 'sed "s/^chr//g" %s | awk \'$3 != \".\"\' | grep "^[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X]" | cut -f1 | sort -u | wc -l' % vcfInput
    proc = subprocess.Popen([cmd], stdout = subprocess.PIPE, shell = True)
    (out, err) = proc.communicate()
    if not len(out) < 23:
        flag = 3
        if fileExt == 'gz':
            subprocess.call(['rm', 'temp.vcf'])
        sys.exit('Not all chromosomes are annotated!')

def process_input(args):
    (infile, jobID, type, outPath, plink) = args
    rawInput = jobID + '_raw_data'
    if type == 'vcf':
        print('Check VCF input...')
        check_vcf(infile)
        print('Convert VCF to 23andme...')
        fileExt = infile.split(".")[-1]
        if fileExt == 'gz':
            extractCmd = 'bgzip -dc %s | sed "s/^chr//g" | awk \'$3 != \".\"\' | grep "^[#,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,(MT)]" > %s/imputation/%s/tmp.vcf' % (infile, outPath, jobID)
            subprocess.call([extractCmd], shell = True)
        else:
            copyInput = 'sed "s/^chr//g" | awk \'$3 != \".\"\' | grep "^[#,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,(MT)]" > %s/imputation/%s/tmp.vcf' % (infile, outPath, jobID)
            subprocess.call([copyInput], shell = True)
        plinkCMD = '%s --vcf %s/imputation/%s/tmp.vcf --snps-only --recode 23 --out %s/imputation/%s/%s' % (plink, outPath, jobID, outPath, jobID, rawInput)
        removeTmpVCF = 'rm %s/imputation/%s/tmp.vcf' % (outPath, jobID)
        subprocess_cmd((plinkCMD, removeTmpVCF))
    elif type == '23andme':
        cpCMD = 'cp %s %s/imputation/%s/%s\.txt' % (infile, outPath, jobID, rawInput)
        subprocess.run(cpCMD, shell = True)

def prepare_dir(args):
    (infile, mode, jobID, type, outPath, plink) = args
    if mode == 'full' or mode == 'impute':
        try:
            my_abs_path = Path(outPath + '/imputation/' + jobID).resolve(strict=True)
        except FileNotFoundError:
            Path(outPath + '/imputation/' + jobID).mkdir(parents = True, exist_ok = True)
        # [convert format and] copy infile to imputation folder
        process_input([infile, jobID, type, outPath, plink])
        # prepare folders
        runDir = outPath + '/imputation/' + jobID
        outputDir = outPath + '/output/' + jobID
        Path(outputDir).mkdir(parents = True, exist_ok = True)
    elif mode == 'test':
        # get existing job ID
        outputDirTmp = outPath + '/output/' + jobID
        checkFileExist(outputDirTmp)
        runDir = outPath + '/imputation/' + jobID
        outputDir = outPath + '/output/' + jobID
    return(runDir, outputDir)

def main():
    version = "1.0.2"
    parser = argparse.ArgumentParser(description="You are running runImputation.py version " + str(version) + ".")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('additional arguments')
    required.add_argument('-i', '--infile', help='Input file in vcf or 23andme format', action='store', default='', required=True)
    required.add_argument('-t', '--type', help='Type of input (vcf|23andme)', action='store', choices=['vcf', '23andme'], default='', required=True)
    required.add_argument('-o', '--outPath', help='Output folder', action='store', default='', required=True)
    optional.add_argument('-m', '--mode', help='Run mode (full|impute|test). Full will run both imputation and perform genetic testing.'
                                        'Test mode will run only the tests based on pre-computed imputation data. Default: full',
                                        choices=['full', 'impute', 'test'], action='store', default='full')
    optional.add_argument('-n', '--id', help='Job ID. If not given, a random string will be generated.', action='store', default='')
    optional.add_argument('-c', '--cleanup', help='Clean up imputation files.', action='store_true')
    optional.add_argument('--trait', help='Trait names for report. If not set, all available traits will be reported.', action='store', default='all')
    args = parser.parse_args()

    # get arguments and config paths
    cfg = load_config('config.yml')
    type = args.type
    mode = args.mode
    infile = args.infile
    checkFileExist(infile)
    outPath = args.outPath
    jobID = args.id
    if (jobID == ''):
        jobID = 'id_' + randomStringDigits(9)
    cleanup = args.cleanup
    trait = args.trait
    if not trait == "all":
        checkFileExist(cfg['imputeTraits'] + "/" + trait)

    scriptPath = os.path.abspath(os.path.dirname(sys.argv[0]))

	# do imputation
    (runDir, outputDir) = prepare_dir([infile, mode, jobID, type, outPath, cfg['plink']])
    imputeCMD = 'Rscript %s/imputation_a2z.R %s %s %s %s %s %s %s %s %s %s %s %s' % (scriptPath, trait, mode, jobID, runDir, outputDir, cfg['shapeit'], cfg['plink'], cfg['gtool'], cfg['impute2'], cfg['sample_ref'], cfg['imputeDataDir'], cfg['imputeTraits'])
    # print(imputeCMD)
    subprocess.run(imputeCMD, shell = True)
    if cleanup:
        rmCMd = 'rm -rf %s' % (runDir)
        subprocess.run(rmCMd, shell = True)
    print('Finished! Check outputs in %s folder.' % (outputDir))

if __name__ == '__main__':
    main()
