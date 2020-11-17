import matplotlib.pyplot as plt
import sys
import numpy as np


with open(sys.argv[1], 'r') as pair_file:
    lines = pair_file.readlines()
    x = []
    y = []
    for line in lines:
        tmp = line.replace("\n", "").split(" ")
        x.append(int(tmp[0]))
        y.append(float(tmp[1]))
    plt.plot(x, y)
    plt.show()

     
#print(m)
