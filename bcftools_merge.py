import glob
import os
import subprocess as sp
import natsort



############################################################

INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/'
LIST_FORMAT = r'*.list'

THREAD = 4

OUTPUT_DIR = INPUT_DIR

is_make_header_removed_file = True

############################################################

input_path_list = glob.glob(INPUT_DIR + LIST_FORMAT)
input_path_list = natsort.natsorted(input_path_list)

for i in input_path_list:
    group_name = i.split('.')[-2].split(r'/')[-1]
    type_name = group_name.split(r'_')[0]
    
    if type_name == 'SNP':
        type_option = 'snps'
    elif type_name == 'INDEL':
        type_option = 'indels'
    
    output_vcf = OUTPUT_DIR + f'{group_name}.vcf'

    sp.call(f"bcftools merge -l {i} -m {type_option} -o {output_vcf} --threads {THREAD}", shell=True)
    
    if is_make_header_removed_file == True:
        rm_header_vcf = OUTPUT_DIR + 'rmHeader_' + f'{group_name}.vcf'
        sp.call(f"sed '/^##/d' {output_vcf} > {rm_header_vcf}", shell=True)



