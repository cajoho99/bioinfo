#!/usr/bin/env python
import matplotlib.pyplot as plt
import sys
import numpy as np


class Atom(object):
    def __init__(self, serial=None, x=0, y=0, z=0):
        self.serial = serial
        self.x = x
        self.y = y
        self.z = z

    def distance(self, other):
        return np.sqrt(np.power(self.x - other.x, 2) + np.power(self.y - other.y, 2) + np.power(self.z - other.z, 2))

    def __str__(self):
        return "Atom " + str(self.serial) + " has the coordinates (" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) + ")"

    def __repr__(self):
        return str(self) + "\n"


def process_input(lines):
    atoms = []
    for line in lines:
        data = line.strip().split("\t")
        seq = int(data[0])
        coords = data[1].split(" ")
        atom = Atom(seq, float(coords[0]), float(coords[1]), float(coords[2]))
        atoms.append(atom)
    result = identify_order(atoms)
    print(result)
    plot_result(atoms, result)


def within_distance_of(dist, size, diff):
    return dist < size+diff and dist > size-diff


def find_start_index(atoms):
    size = 3.8
    diff = 0.1
    combinations = [0 for _ in range(len((atoms)))]
    for i in range(len(atoms)):
        for j in range(i, len(atoms)):
            if i != j and within_distance_of(atoms[i].distance(atoms[j]), size, diff):
                combinations[i] += 1
                combinations[j] += 1
    return [i for i in range(len(combinations)) if combinations[i] == 1]


def identify_order(atoms):
    start_index_result = find_start_index(atoms)
    root_index = start_index_result[0]
    end_index = start_index_result[1]

    visited = [atoms[root_index].serial]
    current = atoms[root_index]
    while len(visited) != len(atoms):
        for atom in atoms:
            if atom.serial not in visited and within_distance_of(current.distance(atom), 3.8, 0.1):
                visited.append(atom.serial)
                current = atom

    return visited


def plot_result(atoms, result):
    x = []
    y = []
    z = []
    for atom in atoms:
        x.append(atom.x)
        y.append(atom.y)
        z.append(atom.z)

    resX = []
    resY = []
    resZ = []
    for r in result:
        resX.append(atoms[r-1].x)
        resY.append(atoms[r-1].y)
        resZ.append(atoms[r-1].z)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, c="black")
    ax.scatter([atoms[1].x, atoms[7].x], [atoms[1].y, atoms[7].y],
               [atoms[1].z, atoms[7].z], c="yellow", s=200)
    ax.plot(resX, resY, resZ, c="red")

    # ax.plot(x, y, z)
    plt.show()


with open(sys.argv[1], 'r') as input_file:
    process_input(input_file.readlines())
