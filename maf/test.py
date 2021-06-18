# pbs test

import subprocess as sp


## config
pbs_N = "maf.test"
pbs_o = "/home/pbsuser/pbs_backup/maf_backup/"
pbs_j = "oe"
pbs_l_core = 2

SRC_DIR = r"/data_244/py_src/VCF_with_py/"

# print('hi')

for i in range(5):
    sp.call(f'echo "python {SRC_DIR}maf/qsub_test_worker.py" | qsub -N {pbs_N} -o {pbs_o} -j {pbs_j} -l select={pbs_l_core} &', shell=True)



sp.call(f'qstat', shell=True)