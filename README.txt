
SNP, INDEL로 나뉜 vcf파일을

bcftools_isec.py로 공통, 차이 변이로 각각 나눠주고,


counting_isec_data_and_mk_df.py로 각 디렉토리 내 파일들의 줄수(변이 포지션 수)를 
텍스트 파일로 옮겨줌


bcftools_mk_subset.py로 vcf파일을 서브셋화 시킬 수 있다

bcftools_vcf_filter_freq.py로 어떤 요소가 필터에 빈번하게 걸렸는지 확인할 수 있다.


gather_specific_isec_data.py로 차집합 or 교집합에 해당하는 
vcf파일들을 모아서(0000.vcf같은것) 따로 Teratoma_specifics같은 디렉토리에 모아 
저장할 수 있다. 인덱싱도 동시진행


extract_T_only.py로 1) Tonly와 2) T-specific이긴 한데 Origin이 탈락한 부분을 
얻고 / 비율구하고 / 걸린 필터 비율 구하기

dp30filter.py로 DP 30미만인 포지션에 DP30 찍어줌. 

mk_pass_vcf.py로 pass만 존재하는 vcf파일을 생성할 수 있다. 리눅스 커맨드로 코딩

maf_handling_test.py