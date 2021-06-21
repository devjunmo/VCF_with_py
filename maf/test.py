
tmp_lst = ['my_1', 'my_2']

# tmp_lst = list(map(split, ))


tmp_lst = [c.split('_')[0] for c in tmp_lst]

print(tmp_lst)