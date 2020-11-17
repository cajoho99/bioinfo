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

all:	 out/make_distance_map out/domain_partition

out/make_distance_map: src/make_distance_map.c
	gcc $(CFLAGS) src/make_distance_map.c -o make_distance_map -lm 

out/domain_partition: src/domain_partition.c
	gcc $(CFLAGS) src/domain_partition.c -o domain_partition  -lm

clean:
	rm -f domain_partition
	rm -f make_distance_map
