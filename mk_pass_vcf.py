import os
from glob import glob
import subprocess as sp

# 주석 이후부분이 안됨.. 다음에 다시 보자
# 이거 목적은 maf파일 만들때 PASS가 아닌게 껴있어서 짜증난다인데,
# 생각해보니 어짜피 pandas로 가져와서 p/np 솎아낼 수 있음


input_dir = r'/myData/WES/data/vcf/hard/WES1_210420/Teratoma_sample_vcf/'
input_format = r'hardFiltered_*_Teratoma*.vcf'
output_dir_name = r'teratoma_pass_vcf/'

output_dir = input_dir + output_dir_name

if os.path.isdir(output_dir) is False:
    os.mkdir(output_dir)


input_lst = glob(input_dir + input_format)

for i in range(len(input_lst)):

    f_name = input_lst[i].split(r'/')[-1]
    
    input_vcf = input_lst[i]
    output_vcf = output_dir + 'pass' + f_name
    
    sp.call(r"awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == 'PASS') print}' %s > %s" % (input_vcf, output_vcf)  \
    , shell=True)