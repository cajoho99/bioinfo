#
# Instructions
# ------------
#
# To compile the programs 'pdb_io', 'atom_array' and 'residue_array', type:
#
#	make
# or
#	make all
# or
#	make pdb_io
#	make atom_array
#	make residue_array
#
# To remove those files that can be recompiled from the source code, type:
#
#	make clean
#
#
#
# Description
# -----------
#
# This 'Makefile' contains instructions for creating the executable programs
# 'pdb_io', 'atom_array' and 'residue_array' from the C source file 'pdb_io.c',
# 'atom_array.c' and 'residue_array.c' using the GNU C compiler 'gcc'.
#
# When the user types 'make' at the UNIX prompt, the 'make' program looks
# for the first target in the 'Makefile', and builds that target.
# Targets are given at the start of a line, and are followed by a colon (':')
# then a tab character(s), then a list of files on which the target depends
# (the "dependencies").  The target will be built from the dependencies using
# the command(s) on the following line(s).
# If the target already exists, and is newer than all of the files on which
# it depends, then the target will not be rebuilt.
# For more information, see the UNIX manual page for 'make' (type 'man make').
#
# In this case, target 'all' depends on the executable files 'pdb_io',
# 'atom_array' and 'residue_array'. These in turn depend on the C source files.
#
# Compiler options used:
#
#	-Wall	report possible programming errors in the C source
#	-g	produce debugging information
#	-o	put the output (the executable file) in named file
#

CFLAGS= -Wall -g

all:	 make_distance_map domain_partition

make_distance_map: src/make_distance_map.c
	gcc $(CFLAGS) src/make_distance_map.c -o make_distance_map -lm 

domain_partition: src/domain_partition.c
	gcc $(CFLAGS) -D RESULT src/domain_partition.c -o domain_partition  -lm

domain_partition_debug: src/domain_partition.c
	gcc $(CFLAGS) -D DEBUG src/domain_partition.c -o domain_partition  -lm


domain_partition_graph: src/domain_partition.c
	gcc $(CFLAGS) -D GRAPH src/domain_partition.c -o domain_partition  -lm
	./domain_partition ./data/1cdh.pdb > ./out/1cdh.map
	./domain_partition ./data/2csn.pdb > ./out/2csn.map
	python3 ./display_scripts/graph_split.py ./out/1cdh.map &
	python3 ./display_scripts/graph_split.py ./out/2csn.map &

strudl: src/strudl_domain_partition.c
	gcc $(CFLAGS) src/strudl_domain_partition.c -o strudl -lm


clean:
	rm -f domain_partition
	rm -f make_distance_map
	rm -f strudl
