
# venn_geneList 가져와서 caseNum pop시키고 maf파일 형태로 정렬해서 만들기

import pandas as pd 
import os


# root_dir = r'E:/stemcell_ips/HAP_MUT2_compare/unfiltered/tech2/nonpass_maf'
# root_dir = r'D:/stemcell/hg38/tech_comp/ips29-A-p49/unfiltered/DP_filtered/tech2/filtered_variant'
# root_dir = r'D:/stemcell/hg38/tech_comp/ips29-A-p49/unfiltered/DP_filtered/tech2/only_dpfilter_varinat'
root_dir = r'D:/stemcell/hg38/tech_comp/ips29-A-p49/unfiltered/for_opponent_maf/tech2'


gene_lst_name = r'ips_A_p49_tech2_comp_unfiltered_DP-tag.xlsx'
# gene_lst_name = r'ips_A_p49_tech1_muthap_comp_unfilter.xlsx'

# output_name = r'ips_A_p49_tech2_varinat_filtered.maf'
output_name = r'ips-A-p49-tech2-varinat-unfiltered-DPTag.maf'

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