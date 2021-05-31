import vcfpy
import os
import glob
import pandas as pd



print(os.getcwd())

# exit(0)

# hard_path = r'/myData/WES/data/vcf/hard/WES1_210420/'
# cnn_path = r'/myData/WES/data/vcf/cnn/WES1_210420/'
test_path = r'/myData/WES/data/vcf/test'

input_test_indel_vcf_name = r'hardFiltered_INDEL_Teratoma-10.vcf'
reader_indels = vcfpy.Reader.from_path(test_path + input_test_indel_vcf_name)

input_test_snp_vcf_name = r'hardFiltered_SNP_Teratoma-10.vcf'
reader_snp = vcfpy.Reader.from_path(test_path + input_test_snp_vcf_name)

header = ['#CHROM', 'POS', 'REF', 'ALT'] + reader_indels.header.samples.names # #CHROM  POS     REF     ALT     Teratoma-10
print('\t'.join(header))

output_file_name = 'output_test.vcf'


# Open input, add FILTER header, and open output file
print("Open input, add FILTER header, and open output file")
reader = vcfpy.Reader.from_path(test_path + input_test_indel_vcf_name)

reader.header.add_filter_line(vcfpy.OrderedDict([('ID', 'DP10'), ('Description', 'total DP < 10')]))

writer = vcfpy.Writer.from_path(test_path + output_file_name, reader.header)


# Add "DP10" filter to records having less than 10 reads
for record in reader:
    ad = sum(c.data.get('DP', 0) for c in record.calls)
    if ad < 10:
        record.add_filter('DP10')
    writer.write_record(record)
