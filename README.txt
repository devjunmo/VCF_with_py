## vcf / bcftools 활용

SNP, INDEL로 나뉜 vcf파일을

bcftools_isec.py로 공통, 차이 변이로 각각 나눠주고,


counting_isec_data_and_mk_df.py로 각 디렉토리 내 파일들의 줄수(변이 포지션 수)를 
텍스트 파일로 옮겨줌


bcftools_mk_subset.py로 vcf파일을 서브셋화 시킬 수 있다

bcftools_vcf_filter_freq.py로 어떤 요소가 필터에 빈번하게 걸렸는지 확인할 수 있다.


gather_specific_isec_data.py로 차집합 or 교집합에 해당하는 
vcf파일들을 모아서(0000.vcf같은것) 따로 Teratoma_specifics같은 디렉토리에 모아 
저장할 수 있다. 인덱싱도 동시진행


extract_T_only.py -> !! 개념이 틀림 !!

dp30filter.py로 DP 30미만인 포지션에 DP30 찍어줌. 

mk_pass_vcf.py로 pass만 존재하는 vcf파일을 생성할 수 있다. 리눅스 커맨드로 코딩

csv2bed로 bed파일생성


sort_bedfile.py로 chr, start, 기준으로 소팅한다.

bcftools_isec_WGS_WES.py  // gather_specific_isec_data_WGS_WES.py 
WES와 WGS를 비교하는 용도로 작성.



## maf 파일 다루기



## T-only 핸들링