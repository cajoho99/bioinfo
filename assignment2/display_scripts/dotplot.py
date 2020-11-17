import matplotlib.pyplot as plt
import sys
import numpy as np

maxX = 0
maxY = 0
m = []
with open(sys.argv[1], 'r') as pair_file:
     lines = pair_file.readlines()
     coords = []
     for line in lines:
          tmp = line.replace("\n", "").split(" ")
          coords.append([int(tmp[0]), int(tmp[1])])
          if int(tmp[0]) > maxX: maxX = int(tmp[0])
          if int(tmp[1]) > maxY: maxY = int(tmp[1])
     print(maxX, maxY)
     for y in range(maxY):
          row = []
          for x in range(maxX):
               row.append(0)
          m.append(row)
     for coord in coords:
          print(coord)
          m[coord[0]-1][coord[1]-1] = 1
          #m[coord[0], coord[1]] = 1


     
#print(m)
plt.imshow(m)
plt.show()