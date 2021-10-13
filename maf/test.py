
tmp_lst = ['my_1', 'my_2']

# tmp_lst = list(map(split, ))


tmp_lst = [c.split('_')[0] for c in tmp_lst]

print(tmp_lst)



in_test_lst = ['apple', 'banana', 'cherry']

if 'apple' not in in_test_lst:
    print('사과는 있다')
else:
    print('사과는 없다')