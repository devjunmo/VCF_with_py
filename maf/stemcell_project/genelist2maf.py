
# venn_geneList 가져와서 caseNum pop시키고 maf파일 형태로 정렬해서 만들기

import pandas as pd 
import os


# root_dir = r'E:/stemcell_ips/HAP_MUT2_compare/unfiltered/tech2/nonpass_maf'
root_dir = r'E:/stemcell_ips/HAP_MUT2_compare/unfiltered/just_ac_apply_comp/hips29A-p49-tech2'

gene_lst_name = r'ips_A_p49_tech2_justAC.xlsx'

output_name = r'ips_A_p49_tech2_justAC.maf'

maf_col_name = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', \
                'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Classification', \
                'FILTER', 't_depth', 't_ref_count', 't_alt_count']

sheet_name = 'Gene data'


gene_lst_path = os.path.join(root_dir, gene_lst_name)


output_path = os.path.join(root_dir, output_name)


gene_lst_df = pd.read_excel(gene_lst_path, sheet_name = sheet_name)


gene_lst_df.pop('Case_number')


gene_lst_df = gene_lst_df.sort_values(by=['Chr', 'Start'])
gene_lst_df.reset_index(inplace=True, drop=True)

gene_lst_df.columns = maf_col_name



print(gene_lst_df)


gene_lst_df.to_csv(output_path, sep='\t', index=False)