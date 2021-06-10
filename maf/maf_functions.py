
def test():
    print('test')


def var_ref_df_to_dict(_refDf, _value_col_name='var_class'):
    _var_class_dict = _refDf.to_dict() # 딕셔너리 안의 딕셔너리
    # print(_var_class_dict)
    # print(_var_class_dict[_value_col_name])
    return _var_class_dict[_value_col_name]