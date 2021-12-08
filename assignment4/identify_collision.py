#!/usr/bin/env python
import matplotlib.pyplot as plt
import sys
import numpy as np


#
# Author: Carl Holmberg
# File: identify_collision.py
# Purpose: identifies the two collidiing atoms between two pdb files.
#

first = []
second = []

comparisons = 0

counter = 0;

class Point(object):
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
    
    def distance_to(self, other):
        global comparisons
        comparisons += 1
        return np.sqrt(np.power(self.x - other.x, 2) + np.power(self.y - other.y, 2) + np.power(self.z - other.z, 2))

class Atom(object):
    def __init__(self, atom_nr, point):
        self.atom_nr = atom_nr
        self.point = point

def build_array(lines):
    min_x = 1000000000
    max_x = -1000000000
    min_y = 1000000000
    max_y = -1000000000
    min_z = 1000000000
    max_z = -1000000000

    atoms = []
    global first

    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            split = [a for a in line.split(" ") if a != '']
            #print(split)
            x = float(split[6])
            y = float(split[7])
            z = float(split[8])
            point = Point(x,y,z)
            first.append(point)
            int_x = int(np.floor(x))
            int_y = int(np.floor(y))
            int_z = int(np.floor(z))
            atoms.append((int_x, int_y, int_z, point))
            if max_x < int_x:
                max_x = int_x
            
            if max_y < int_y:
                max_y = int_y
            
            if max_z < int_z:
                max_z = int_z
            
            if min_x > int_x:
                min_x = int_x

            if min_y > int_y:
                min_y = int_y
            
            if min_z > int_z:
                min_z = int_z

            #print(str(int_x), str(int_y), str(int_z))

    #print("min_x:" + str(min_x) + " | max_x: " + str(max_x))
    #print("min_y:" + str(min_y) + " | max_y: " + str(max_y))
    #print("min_z:" + str(min_z) + " | max_z: " + str(max_z))

    x_offset = 0 - min_x
    y_offset = 0 - min_y
    z_offset = 0 - min_z
    
    #print("x_offset: ", str(x_offset), " | y_offset: ", str(y_offset), " | z_offset: ", str(z_offset))
    #print("number of atoms in first file:", str(len(atoms)))

    ds = [[[[] for _ in range(min_z+z_offset, max_z+z_offset+1)] for _ in range(min_y+y_offset, max_y+y_offset+1)] for _ in range(min_x + x_offset, max_x + x_offset+1)]
    for (x, y, z, p) in atoms:
        #print(x, y, z)
        ds[x+x_offset][y+y_offset][z+z_offset].append(p)

    return (ds, x_offset, y_offset, z_offset)


def search_area(ds, b_x, b_y, b_z, atom):
    check_distance = []

    # Get dimensions
    d_x = len(ds)
    d_y = len(ds[0])
    d_z = len(ds[0][0])

    #print("d_x: " + str(d_x) + " | d_y: " + str(d_y) + " | d_z: " + str(d_z))


    
    for x in [0,1,-1,2,-2,3,-3,4,-4]:
        for y in [0,1,-1,2,-2,3,-3,4,-4]:
            for z in [0,1,-1,2,-2,3,-3,4,-4]:
                if not ((x == y and y == z) and (x == 4 or x == -4)):
                    if b_x+x>= 0 and b_x+x < d_x and b_y+y >= 0 and b_y+y < d_y and b_y+y >= 0 and b_z+z < d_z:
                        #print(str(b_x+x), str(b_y+y), str(b_z+z))
                        current = ds[b_x+x][b_y+y][b_z+z]
                        if len(current) > 0:
                            check_distance.extend(current)
                            break
    
    #print(check_distance)
    global counter

    for p in check_distance:
        #print(p)
        d = atom.point.distance_to(p)
        #print("Distance for atom " + str(atom.atom_nr) + " was " + str(d))
        if d <= 4:
            return atom.atom_nr
        else:
            counter += 1

    return None




def find_collisions(lines, ds, o_x, o_y, o_z):
    colliding = []
    global second
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            split = [a for a in line.split(" ") if a != '']
            #print(split)
            a_nr = int(split[5])
            x = float(split[6])
            y = float(split[7])
            z = float(split[8])
            point = Point(x,y,z)
            second.append(point)
            atom = Atom(a_nr, point)
            int_x = int(np.floor(x))
            int_y = int(np.floor(y))
            int_z = int(np.floor(z))

            c = search_area(ds, int_x + o_x, int_y + o_y, int_z + o_z, atom)
            if c != None:
                colliding.append(c)
    return colliding


def plot_debug():
    global first
    global second
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter([p.x for p in first], [p.y for p in first], [p.z for p in first], c="black")
    ax.scatter([p.x for p in second], [p.y for p in second], [p.z for p in second], c="red")
    plt.show()



with open(sys.argv[1], 'r') as input_file:
    processed = build_array(input_file.readlines())
    ds = processed[0]
    o_x = processed[1]
    o_y = processed[2]
    o_z = processed[3]
    with open(sys.argv[2], 'r') as comp_file:
        result = find_collisions(comp_file.readlines(), ds, o_x, o_y, o_z)
        for res in result:
            print(res)
        print("Collisions: " + str(len(result)))
        print("Comparisons: " + str(comparisons))
        #print("Bad comparisons: " + str(counter))
        plot_debug()
