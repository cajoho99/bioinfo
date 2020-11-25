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
    alpha_carbon = from_seq_to_atoms(filter_atoms(atoms), atoms)
    filtered = any_legal(alpha_carbon, atoms)
    print("any legal: ")
    print(filtered)
    result = identify_order(filtered)
    plot_result(atoms, result)


def within_distance_of(dist, size, diff):
    return dist > 0.0001 and dist < size+diff and dist > size-diff


def any_within_distance(atom, atoms):
    for a in atoms:
        if within_distance_of(atom.distance(a), 3.8, 0.1):
            return True
    return False


def any_legal(candidates, atoms):
    out = []
    for i in candidates:
        for j in candidates:
            if i != j:
                print("i.serial: " + str(i.serial) +
                      " j.serial: " + str(j.serial))
                center = (i.x + (i.x - j.x) / 2, i.y +
                          (i.y - j.y) / 2, i.z + (i.z - j.z) / 2)
                a_from = atoms_from(center, 3.8, atoms)
                if(len(a_from) > 4):
                    if i.serial not in out:
                        out.append(i.serial)
                    if j.serial not in out:
                        out.append(j.serial)


def atoms_from(xyz, r, atoms):
    out = []
    for a in atoms:
        if distance_to(xyz[0], xyz[1], xyz[2], a) < r:
            out.append(a)
    return out


def distance_to(x, y, z, atom):
    return np.sqrt(np.power(x - atom.x, 2) +
                   np.power(y - atom.y, 2) +
                   np.power(z - atom.z, 2))


def filter_atoms(atoms):
    return [a.serial for a in atoms if any_within_distance(a, atoms)]


def from_seq_to_atoms(seq, atoms):
    return [atoms[s-1] for s in seq]


def find_start_index(atoms):
    size = 3.8
    diff = 0.1
    combinations = [0 for _ in range(len(atoms))]
    for i in range(len(atoms)):
        for j in range(i, len(atoms)):
            if i != j and within_distance_of(atoms[i].distance(atoms[j]), size, diff):
                combinations[i] += 1
                combinations[j] += 1
    return [i for i in range(len(combinations)) if combinations[i] == 1]


def search_atom(atoms, current, excluded):
    for atom in atoms:
        print("c")
        if atom.serial not in excluded:
            print("a")
            if within_distance_of(current.distance(atom), 3.8, 0.1):
                print("b")
                return (True, atom)
    return (False, current)


def identify_order(atoms):
    start_index_result = find_start_index(atoms)
    root_index = start_index_result[0]
    end_index = start_index_result[1]

    visited = [atoms[root_index].serial]
    current = atoms[root_index]
    while (len(visited) != len(atoms)):
        res = search_atom(atoms, current, visited)
        if res[0]:
            visited.append(res[1].serial)
            current = res[1]

    return visited


def plot_result(atoms, result):
    print("number of residues: " + str(len(result)))

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
    plt.show()


with open(sys.argv[1], 'r') as input_file:
    process_input(input_file.readlines())
