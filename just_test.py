from enum import Enum

tmp = 1

print(type(str(tmp)))

print(type(tmp))


 a = {'name':'pey', 'phone':'0119993323', 'birth': '1118'}
 print(a)


 for i in a.items():
    print(i[0])
    print(i[1])

# exit(0)

class TeratomaOrigin(Enum): # csv 파일의 컬럼명 순서로 구성
    Teratoma = 0
    orgin = 1

e = TeratomaOrigin
e.Teratoma.value
e.Teratoma.name

# exit(0)
a = 1
a