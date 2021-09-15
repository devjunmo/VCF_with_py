
test_null = None

if test_null is None:
    print('null!!')



tmp_lst = [1, 2, 3]

for i in tmp_lst:
    if i == (len(tmp_lst)-1):
        print(i)



sub_lst = [[('a', 'b'), ('c', 'd')], [('e', 'f'), ('g', 'h')], [('i', 'j'), ('k', 'l')]]

print(sub_lst[0][0][1])

for i in range(len(sub_lst)):

    isec_data = eval(",".join([f"'sub_lst[{i}][{j}][1]'" for j in range(len(sub_lst[i]))]))

    print(isec_data)