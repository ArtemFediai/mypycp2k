import numpy as np
import csv
import os

folder_to_save_to = 'csv'

if not os.path.exists(folder_to_save_to):
    os.mkdir(folder_to_save_to)

last_number = 133885
first_number = 1

step = 1000

num_steps = last_number // 1000 + 1
last_iter = last_number - (num_steps - 1) * step
print(last_iter)

for i in range(num_steps):
    i1 = 1000*i +1
    i2 = 1000*(i+1)
    interval = np.linspace(i1, i2, i2-i1+1, dtype=int)
    #print(interval)
    print(len(interval))
    if i == num_steps-1:
        interval = np.linspace(i1, last_number, last_number - i1, dtype=int)
    with open(folder_to_save_to + '/' + str(i) + '.csv', mode='w') as stream:
        writer = csv.writer(stream)
        writer.writerow(['{:0>6}'.format(x) for x in interval])
print(interval)