import os
from glob import glob
import subprocess as sp


# input_dir = r'/data_244/VCF/gatherd_WGS_WES/WES_specific/'
# input_dir = r'/data_244/VCF/gatherd_WGS_WES/common/'
# input_dir = r'/data_244/VCF/gatherd_WGS_WES/WGS_specific/'
# input_dir = r'/data_244/WGS/HN00146173/gs/hardfiltered/'
# input_dir = r'/data_244/VCF/gatherd_WGS_WES_interval_apply/WES_specific/'
# input_dir = r'/data_244/VCF/gatherd_WGS_WES_interval_apply/WGS_specific/'
# input_dir = r'/data_244/VCF/gatherd_noDP_WGS_WES_interval_apply/WGS_specific/'
# input_dir = r'/data_244/WGS/HN00146173/gs/hardfiltered/interval_apply_vcf/'
# input_dir = r'/data_244/stemcell/WES/ips_recal_bam/tumor_only/filtered_vcf/pass_only/'
# input_dir = r'/data_244/stemcell/WES/ips_recal_bam/tumor_only/filtered_vcf/inc_germ/'
# input_dir = r'/data_244/stemcell/WES/ips_recal_bam/vardict_tumor_only/'
# input_dir = r'/data_244/stemcell/WES/ips_recal_bam/vardict_tumor_only/passonly/'
# input_dir = r'/data_244/stemcell/WES/hg38_pp/gs/merged/'
input_dir = r'/data_244/stemcell/WES/hg38_pp/gs/merged/'


input_format = r'*.vcf'
# input_format = r'hardFiltered_SNP_*'
output_suffix = r'_germline'

output_dir_name = r'maf/'
tmp_dir = input_dir + 'vep_vcf/'
# fasta_path = r'/data_244/refGenome/b37/human_g1k_v37.fasta'
fasta_path = r'/data_244/refGenome/hg38/v0/Homo_sapiens_assembly38.fasta'

SRC_DIR = r"/home/pbsuser/mskcc-vcf2maf-754d68a/"
SRC_PATH = SRC_DIR + "vcf2maf.pl"


## pbs config
pbs_N = "ips_ant_germline"
pbs_o = input_dir + "qsub_log/"
pbs_j = "oe"
pbs_l_core = 2


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

    # f_name = input_lst[i].split(r'/')[-1].split(r'.')[0].split(r'_')[1] # teratoma-4
    # f_name = input_lst[i].split(r'/')[-1].split(r'.')[0].split(r'_')[1]
    f_name = input_lst[i].split(r'/')[-1].split(r'.')[0].split(r'_')[0]
    # f_type = input_lst[i].split(r'/')[-1].split(r'.')[0].split(r'_')[-2] # snp/indel
    
    input_vcf_path = input_lst[i]
    output_maf_path = output_dir + f_name + output_suffix + '.maf'

    # sp.call(rf"perl vcf2maf.pl --input-vcf {input_vcf_path} --output-maf {output_maf_path} --ref-fasta {fasta_path} --tmp-dir {tmp_dir}", shell=True)

    sp.call(f'echo "perl {SRC_PATH} --input-vcf {input_vcf_path} --output-maf {output_maf_path} --ref-fasta {fasta_path} --tmp-dir {tmp_dir} \
                --tumor-id {f_name} --vcf-tumor-id {f_name}" \
              | qsub -N {pbs_N} -o {pbs_o} -j {pbs_j} -l select={pbs_l_core} &', shell=True)
