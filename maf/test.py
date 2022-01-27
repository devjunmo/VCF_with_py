

import pandas as pd

# fruit_dict = {
#     3: ['apple', 12],
#     2: ['banana', 11],
#     6: ['mango', 10]}

fruit_dict = {
    3: 'apple:10',
    2: 'banana, 11',
    6: 'mango, 10'}

print(pd.DataFrame(list(fruit_dict.items()),
                   columns=['Quantity', 'FruitName']))