
# venn_geneList 가져와서 caseNum pop시키고 maf파일 형태로 정렬해서 만들기

import pandas as pd 
import os



# root_dir = r'D:/stemcell/hg38/clone_comp/hiPS35/A/som_germ_merge'
# root_dir = r'D:/stemcell/hg38/clone_comp/hiPS29/som_germ_merge/unfiltered'
# root_dir = r'E:/stemcell_ips/gdc/tech/29A/filtered/tech2/filter_som_germ_merge'
# root_dir = r'E:/stemcell_ips/gdc/tech/29B/29B_tech1_muthap/filter_mut_hap_merge'
# root_dir = r'E:/stemcell_ips/gdc/tech/29B/29B_tech2_muthap/filter_mut_hap_merge'
# root_dir = r'E:/stemcell_ips/gdc/clone/hips29/29E_muthap/filter_mut_hap_merge'
# root_dir = r'E:/stemcell_ips/gdc/clone/hips66/66C_muthap/filter_mut_hap_merge'
root_dir = r'E:/stemcell_ips/gdc/clone/hips66/66C_muthap/unfilter_mut_hap_merge'


# gene_lst_name = r'hiPS66-C_varinat_filtered.xlsx'
gene_lst_name = r'hiPS66-C_varinat_unfiltered_DPTag.xlsx'


# output_name = r'hiPS66-C_varinat_filtered.maf'
output_name = r'hiPS66-C_varinat_unfiltered_DPTag.maf'


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