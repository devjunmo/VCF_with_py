import vcfpy
import os
import glob
import pandas as pd

### READ
# from_path()로 파일 잡아서 레코드 형태로 가져옴
# fetch()로 원하는 포지션 만큼 잘라서 갖고올 수 있다.
# 



# test

print(os.getcwd())

# exit(0)

hard_path = r'/myData/WES/data/vcf/hard/WES1_210420/'
cnn_path = r'/myData/WES/data/vcf/cnn/WES1_210420/cle'

input_test_indel_vcf_name = r'hardFiltered_INDEL_Teratoma-10.vcf.gz'
reader_indels = vcfpy.Reader.from_path(hard_path + input_test_indel_vcf_name) # vcf를 레코드 단위로 가져옴
reader_indels # <vcfpy.reader.Reader object at 0x7f5b8ea669a0>

input_test_snp_vcf_name = r'hardFiltered_SNP_Teratoma-10.vcf.gz'
reader_snp = vcfpy.Reader.from_path(hard_path + input_test_snp_vcf_name)

header = ['#CHROM', 'POS', 'REF', 'ALT'] + reader_indels.header.samples.names # #CHROM  POS     REF     ALT     Teratoma-10
header
print('\t'.join(header))

# exit(0)

for record in reader_snp:
    # if record.is_snv():
    #     continue
    print(record)
    line = [record.CHROM, record.POS, record.REF]
    print(line)
    # print("record['CHROM'] =", record['CHROM']) 블가능
    
    print(record.ALT) 
    print(type(record.ALT)) # <class 'list'>
    print(record.ALT[0]) # Substitution(type_='INS', value='TCCCTGGAGGACC')
    # print(record.ALT[1]) # IndexError: list index out of range
    print(type(record.ALT[0])) # <class 'vcfpy.record.Substitution'>
    print(record.ALT[0].value) # TCCCTGGAGGACC
    print(record.ALT[0].type) # INS
    
    line += [alt.value for alt in record.ALT]
    print(line)

    print(record.calls) # [Call('Teratoma-10', {'GT': '0/1', 'AD': [18, 13], 'DP': 31, 'GQ': 99, 'PL': [491, 0, 707]})]
    line += [call.data.get('GT') or './.' for call in record.calls] # GT의 값을 받아오거나 ./.을 가져오거나
    line += [call.data.get('AD') or './.' for call in record.calls]
    print(type(record.calls)) # <class 'list'>
    print(record.calls[0]) # Call('Teratoma-10', {'GT': '0/1', 'AD': [18, 13], 'DP': 31, 'GQ': 99, 'PL': [491, 0, 707]})
                           # 붙어있는 샘플단위로 나뉘어진 리스트인듯
    print(record.calls[0].sample) # Teratoma-10
    print(record.calls[0].data) # {'GT': '0/1', 'AD': [18, 13], 'DP': 31, 'GQ': 99, 'PL': [491, 0, 707]}
    print(record.calls[0].called)  # True

    print(line) # ['1', 874950, 'T', 'TCCCTGGAGGACC', '0/1']
    
    mystr = '\t'.join(map(str, line))
    print(mystr) # 1       874950  T       TCCCTGGAGGACC   0/1
    print(type(mystr)) # <class 'str'>
    break

# for record in reader_indels:
#     if record.is_snv():
#         continue
#     line = [record.CHROM, record.POS, record.REF]
#     line += [alt.value for alt in record.ALT]
#     line += [call.data.get('GT') or './.' for call in record.calls] # GT의 값을 받아오거나 ./.을 가져오거나
#     line += [call.data.get('AD') or './.' for call in record.calls]
#     mystr = '\t'.join(map(str, line))

# mystr

# for record in reader_snp:
#     if not record.is_snv():
#         continue
#     line = [record.CHROM, record.POS, record.REF]
#     line += [alt.value for alt in record.ALT]
#     line += [call.data.get('GT') or './.' for call in record.calls]
#     print('\t'.join(map(str, line)))


# exit(0)



# 결론 = vcfpy + pandas 혼합 사용
# pandas -> 모든 레코드들의 pos정보 리스트에 담고, 중복제거 및 소팅 후 1, 2, 3, 4 컬럼에 넣기 
# 확보한 레코드들의 포지션정보에 해당하는 곳에 레코드의 info부분 넣어주기 


