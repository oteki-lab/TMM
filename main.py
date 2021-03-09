import csv
import shutil
import pandas as pd
from itertools import product
import TMM_wrapper as tmm
import bayesian as bo

import_flag = True

init_steps_num = 0

count_limit = 100

height_list = [300, 600, 900, 1200, 1500, 1800]
stack_list = [1, 2, 4, 6, 8, 10]

if import_flag:
    print("load initial steps")
    #init_steps = pd.read_csv('init_steps.csv', encoding='SHIFT-JIS')
    #print(init_steps)
else:
    init_steps = list(product(height_list, stack_list))
    for step in init_steps:
        score = tmm.calc(step[0],step[1])
        row = [step[0],step[1],score[1192]]
        with open('init_steps.csv', 'a', newline="") as f:
            writer = csv.writer(f)
            writer.writerow(row)
shutil.copyfile("init_steps.csv", "hp_steps.csv")

count = 0
while(count_limit>count):
    count = count+1
    print(count)
    next_step = bo.search('.', 'hp_table.csv', 'hp_steps.csv')
    next_step['score'] = tmm.calc(next_step['height'], next_step['stack'])
    with open('hp_steps.csv', 'a', newline="") as f:
        writer = csv.writer(f)
        writer.writerow([next_step['height'], next_step['stack'],next_step['score']])