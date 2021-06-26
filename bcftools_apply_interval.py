import os
import subprocess as sp
from glob import glob

input_dir = r'/data_244/WGS/HN00146173/gs/hardfiltered/'
input_format = r'hardFiltered_SNP_*.vcf.gz'
output_dir_name = r'interval_apply_vcf/'
interval_path = r'/data_244/refGenome/b37/SureSelect_v6_processed.bed'
app_type = r'WGS'

output_dir = input_dir + output_dir_name

if os.path.isdir(output_dir) is False:
    os.mkdir(output_dir)

input_path_lst = glob(input_dir + input_format)

for i in range(len(input_path_lst)):
    input_name = input_path_lst[i].split(r'/')[-1].split(r'.')[0].split('_')[-1] # Teratoma-10
    input_variant_type = input_path_lst[i].split(r'/')[-1].split(r'.')[0].split('_')[-2] # SNP 
    
    output_path = output_dir + app_type + '_' + 'aplyInterval_' + input_variant_type + '_' + input_name + '.vcf.gz'

    sp.call(rf'bcftools view -O z -o {output_path} {input_path_lst[i]} -R {interval_path}', shell = True)

    # break