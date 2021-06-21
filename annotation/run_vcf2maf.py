import os
from glob import glob
import subprocess as sp


# input_dir = r'/data_244/VCF/gatherd_WGS_WES/WES_specific/'
# input_dir = r'/data_244/VCF/gatherd_WGS_WES/common/'
input_dir = r'/data_244/VCF/gatherd_WGS_WES/WGS_specific/'


input_format = r'SNP_*.vcf'
output_dir_name = r'maf/'
tmp_dir = input_dir + 'vep_vcf/'
fasta_path = r'/data_244/refGenome/b37/human_g1k_v37.fasta'

SRC_DIR = r"/home/pbsuser/mskcc-vcf2maf-754d68a/"
SRC_PATH = SRC_DIR + "vcf2maf.pl"


## pbs config
pbs_N = "mk_maf.WES"
pbs_o = input_dir + "qsub_log/"
pbs_j = "oe"
pbs_l_core = 3


output_dir = input_dir + output_dir_name


if os.path.isdir(output_dir) is False:
    os.mkdir(output_dir)

if os.path.isdir(tmp_dir) is False:
    os.mkdir(tmp_dir)

if os.path.isdir(pbs_o) is False:
    os.mkdir(pbs_o)


input_lst = glob(input_dir + input_format)

# os.chdir(r'/root/mskcc-vcf2maf-2235eed')

for i in range(len(input_lst)):

    f_name = input_lst[i].split(r'/')[-1].split(r'.')[0].split(r'_')[-1] # teratoma-4
    f_type = input_lst[i].split(r'/')[-1].split(r'.')[0].split(r'_')[-2] # snp/indel
    
    input_vcf_path = input_lst[i]
    output_maf_path = output_dir + f_type + '_' + f_name + '.maf' 

    # sp.call(rf"perl vcf2maf.pl --input-vcf {input_vcf_path} --output-maf {output_maf_path} --ref-fasta {fasta_path} --tmp-dir {tmp_dir}", shell=True)

    sp.call(f'echo "perl {SRC_PATH} --input-vcf {input_vcf_path} --output-maf {output_maf_path} --ref-fasta {fasta_path} --tmp-dir {tmp_dir}" \
              | qsub -N {pbs_N} -o {pbs_o} -j {pbs_j} -l select={pbs_l_core} &', shell=True)
