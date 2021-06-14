
import os
from glob import glob
import subprocess as sp


# input_dir = r'/myData/WES/data/vcf/hard/WES1_210420/'
# input_format = r'hardFiltered_*_T*.vcf.gz'

input_dir = r'/data_244/WES/hardfilterd_vcf/'
# input_format = r'hardFiltered_*_T*.vcf.gz'
input_format = r'hardFiltered_*_h*.vcf.gz'


output_dir = input_dir + r'dp_apply/'

# sample_type = "SNP" # "INDEL"

seq_type = "WES"

interval_path = r'/data_244/refGenome/b37/SureSelect_v6_processed.bed'

################################################################

input_path_lst = glob(input_dir + input_format)
print(input_path_lst)

print(os.getcwd())
# exit(0)
for i in range(len(input_path_lst)):
    
    sample_type = input_path_lst[i].split(r'/')[-1].split('.')[0].split('_')[1]
    sample_name = input_path_lst[i].split(r'/')[-1].split('.')[0].split('_')[2]
    
    output_path = output_dir + input_path_lst[i].split(r'/')[-1]

    sp.call(rf"sh ./variant_filteration.sh {input_path_lst[i]} {output_path} \
                                        {sample_type} {seq_type} {interval_path}", shell = True)



