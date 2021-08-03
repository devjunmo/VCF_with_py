import os
from glob import glob
import subprocess as sp


input_dir = r'/data_240/utuc/annotaion/'

input_format = r'*.vcf'
# input_format = r'hardFiltered_SNP_*'

output_dir_name = r'maf/'
tmp_dir = input_dir + 'vep_vcf/'
fasta_path = r'/data_240/refGenome/b37/human_g1k_v37.fasta'

SRC_DIR = r"/root/mskcc-vcf2maf-2235eed/"
SRC_PATH = SRC_DIR + "vcf2maf.pl"


output_dir = input_dir + output_dir_name

qsub_config = r'/data_240/qsub.5'

if os.path.isdir(output_dir) is False:
    os.mkdir(output_dir)

if os.path.isdir(tmp_dir) is False:
    os.mkdir(tmp_dir)

# if os.path.isdir(pbs_o) is False:
#     os.mkdir(pbs_o)


input_lst = glob(input_dir + input_format)

# os.chdir(r'/root/mskcc-vcf2maf-2235eed')

for i in range(len(input_lst)):

    sample_name = input_lst[i].split(r'/')[-1].split(r'.')[0] # 20S-31099-A4-5

    caller_name = input_lst[i].split(r'/')[-1].split(r'.')[1] # mutect1
    
    input_vcf_path = input_lst[i]
    output_maf_path = output_dir + sample_name + '_' + caller_name + '.maf' 

    # sp.call(rf"perl vcf2maf.pl --input-vcf {input_vcf_path} --output-maf {output_maf_path} --ref-fasta {fasta_path} --tmp-dir {tmp_dir}", shell=True)

    sp.call(f'perl {SRC_PATH} --input-vcf {input_vcf_path} --output-maf {output_maf_path} --ref-fasta {fasta_path} --tmp-dir {tmp_dir}', shell=True)
