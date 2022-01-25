import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


"""
  root/ --- sample_group1/ --- maf1, maf2, ..
        |
        --- sample_group2/ --- maf1, maf2, ..
        |
        .               .
        .               .

# 기능

1. 그룹별 mutation의 평균치 plotting
2. 전체 mutation count plotting

# 만드는법:

1. root 디렉토리 받는다.
2. 그룹 디렉토리의 이름, path를 따서 리스트화 시킨다.
3. 하위 디렉토리의 maf파일들을 가져와서 df화 시키고 딕셔너리에 {그룹 tag : [평균 count, {sample_name: mut_count}, {}, ..]} 형태로 구조화 시킨다
4. 그룹 tag ~ 평균 count df 만들고 플로팅
5. idx = sample name / col1: mut count / col2: grp df 만들고, 그룹별 플로팅

"""

class GrpMutationProfile:
    
    def __init__(self) -> None:
        pass
