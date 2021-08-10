import subprocess as sp
from glob import glob
import os
import natsort


# SNP와 INDEL로 나뉜 VCF를 합쳐주는 코드

input_dir = r'/data_244/stemcell/WES/germline_vcf'

snp_format = r'hardFiltered_SNP*.vcf.gz'
indel_format = r'hardFiltered_INDEL*.vcf.gz'

output_dir_name = r'concat_vcf'

thread = 10

output_dir = os.path.join(input_dir, output_dir_name)

if os.path.isdir(output_dir) is False:
    os.mkdir(output_dir)


input_snp_lst = glob(os.path.join(input_dir, snp_format))
input_indel_lst = glob(os.path.join(input_dir, indel_format))

input_snp_lst = natsort.natsorted(input_snp_lst)
input_indel_lst = natsort.natsorted(input_indel_lst)

# print(input_snp_lst)
# print(input_indel_lst)




if len(input_snp_lst) == len(input_indel_lst):

    for i in range(len(input_snp_lst)):
        input_snp = input_snp_lst[i]
        input_indel = input_indel_lst[i]

        sample_name = input_snp.split(r'/')[-1].split(r'.')[0].split(r'_')[-1]

        out_name = 'concated-germline_' + sample_name + '.vcf'

        output_path = os.path.join(output_dir, out_name)

        sp.call(rf'bcftools concat --allow-overlaps --threads {thread} -O v -o {output_path} \
                    {input_snp} {input_indel}', shell=True)

else:
    print('len(input_snp_lst) != len(input_indel_lst)')