import pandas as pd
import os
import glob
import subprocess as sp

# _prefix_path_vcf: sample name 직전까지의 prefix 설정
# ex) hardFiltered_INDEL_hiPS68-A.vcf.gz -> _prefix_path_vcf = hardFiltered_INDEL_

class Mk_vcf_intersection():

    def __init__(self, _pair_names_df, _input_dir, \
                _SNP_path_prefix, _INDEL_path_prefix, \
                _root_output_dir_name_snp, _root_output_dir_name_indel, \
                _is_only_PASS=True, _filter_comp="PASS"):

        self.pair_names_df = _pair_names_df
        self.INPUT_DIR = _input_dir
        self._SNP_path_prefix = _SNP_path_prefix
        self._INDEL_path_prefix = _INDEL_path_prefix
        self.root_output_dir_name_snp = _root_output_dir_name_snp
        self.root_output_dir_name_indel = _root_output_dir_name_indel
        self._is_only_PASS = _is_only_PASS
        self.filter_comp = _filter_comp
        

    def run_isec(self):
        for rows in self.pair_names_df.itertuples():
            snp_tumor_data_path = self.INPUT_DIR + self._SNP_path_prefix + rows[1] + '.vcf.gz'
            snp_origin_data = self.INPUT_DIR + self._SNP_path_prefix + rows[2] + '.vcf.gz'

            indel_teratoma_data_path = self.INPUT_DIR + self._INDEL_path_prefix + rows[1] + '.vcf.gz'
            indel_origin_data = self.INPUT_DIR + self._INDEL_path_prefix + rows[2] + '.vcf.gz'

            output_dir_name = rows[1] + '_' + rows[2]

            snp_output_dir = self.INPUT_DIR + self.root_output_dir_name_snp + output_dir_name
            indel_output_dir = self.INPUT_DIR + self.root_output_dir_name_indel + output_dir_name

            if os.path.isdir(snp_output_dir) is False:
                if os.path.isdir(self.INPUT_DIR + self.root_output_dir_name_snp) is False:
                    os.mkdir(self.INPUT_DIR + self.root_output_dir_name_snp)
                os.mkdir(snp_output_dir)
            if os.path.isdir(indel_output_dir) is False:
                if os.path.isdir(self.INPUT_DIR + self.root_output_dir_name_indel) is False:
                    os.mkdir(self.INPUT_DIR + self.root_output_dir_name_indel)
                os.mkdir(indel_output_dir)
            

            if self.is_only_PASS:
                sp.call(f"bcftools isec -f {self.filter_comp} -p {snp_output_dir}/ {snp_tumor_data_path} {snp_origin_data}", shell=True)
                sp.call(f"bcftools isec -f {self.filter_comp} -p {indel_output_dir}/ {indel_teratoma_data_path} {indel_origin_data}", shell=True)
            
            else:
                sp.call(f"bcftools isec -p {snp_output_dir}/ {snp_tumor_data_path} {snp_origin_data}", shell=True)
                sp.call(f"bcftools isec -p {indel_output_dir}/ {indel_teratoma_data_path} {indel_origin_data}", shell=True)



