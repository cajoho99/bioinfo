/*
 * File:	domain_partition.c
 * Author	Carl Holmberg
 * Purpose:	Partitions a protein chain into two parts with parts of the DOMAK algorithm 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_ATOMS 10000
#define LINE_LENGTH 81
#define THRESHOLD 5

//#define DEBUG 1
//#define GRAPH 1
//#define RESULT 1

typedef struct
{
	double x, y, z;
} Point;

typedef struct
{
	int serial;
	char atomName[5];
	char altLoc[2];
	char resName[4];
	char chainID[2];
	int resSeq;
	char iCode[2];
	Point centre;
} Atom;

/*
 * Declare an array to hold data read from the ATOM records of a PDB file.
 */
Atom atom[MAX_ATOMS + 1];

int resA[MAX_ATOMS + 1];
int resB[MAX_ATOMS + 1];

int read_data(char *filename)
{
	FILE *stream;
	char line[LINE_LENGTH];
	char s_serial[6];
	char s_name[5];	  /* Alpha-carbon is " CA "; calcium is "CA  " */
	char s_altLoc[2]; /* Usually " " */
	char s_resName[4];
	char s_chainID[2];
	char s_resSeq[5];
	char s_iCode[2]; /* Usually " " */
	char s_x[9];
	char s_y[9];
	char s_z[9];

	int serial;
	int resSeq;
	double x;
	double y;
	double z;

	int i = 0;

	if ((stream = fopen(filename, "r")) == NULL)
	{
		(void)fprintf(stderr, "Unable to open %s\n", filename);
		exit(0);
	}

	while (fgets(line, LINE_LENGTH, stream))
	{
		if (strncmp(line, "ATOM  ", 6) == 0)
		{

			/*
                         * Split the line into its constituent fields.
                         * We are only interested in columns 1-54.
                         */

			strncpy(s_serial, &line[6], 5);
			s_serial[5] = '\0';
			strncpy(s_name, &line[12], 4);
			s_name[4] = '\0';
			strncpy(s_altLoc, &line[16], 1);
			s_altLoc[1] = '\0';
			strncpy(s_resName, &line[17], 3);
			s_resName[3] = '\0';
			strncpy(s_chainID, &line[21], 1);
			s_chainID[1] = '\0';
			strncpy(s_resSeq, &line[22], 4);
			s_resSeq[4] = '\0';
			strncpy(s_iCode, &line[26], 1);
			s_iCode[1] = '\0';
			strncpy(s_x, &line[30], 8);
			s_x[8] = '\0';
			strncpy(s_y, &line[38], 8);
			s_y[8] = '\0';
			strncpy(s_z, &line[46], 8);
			s_z[8] = '\0';

			/*
                         * Convert the numeric fields to integers or doubles.
                         * The library functions atoi() and atof() are
                         * described in the UNIX manual pages ('man atoi' and
                         * 'man atof').
                         */

			serial = atoi(s_serial);
			resSeq = atoi(s_resSeq);
			x = atof(s_x);
			y = atof(s_y);
			z = atof(s_z);

			/*
			 * Copy values to the next element in the atom array.
			 */

			//printf("NAME: '%s'\n", s_name);
			if (strncmp(s_name, " CA ", 4) == 0)
			{
				if (++i > MAX_ATOMS)
				{
					(void)fprintf(stderr, "Too many atoms read\n");
					exit(0);
				}
				atom[i].serial = serial;
				strcpy(atom[i].atomName, s_name);
				strcpy(atom[i].altLoc, s_altLoc);
				strcpy(atom[i].resName, s_resName);
				strcpy(atom[i].chainID, s_chainID);
				atom[i].resSeq = resSeq;
				strcpy(atom[i].iCode, s_iCode);
				atom[i].centre.x = x;
				atom[i].centre.y = y;
				atom[i].centre.z = z;
			}
		}
	}
	return i;
}

void write_pdb_atom(int serial, char *s_name, char *s_altLoc, char *s_resName, char *s_chainID,
					int resSeq, char *s_iCode, Point centre)
{
	printf("ATOM  %5d %s%s%s %s%4d%s   %8.3f%8.3f%8.3f\n",
		   serial,
		   s_name,
		   s_altLoc,
		   s_resName,
		   s_chainID,
		   resSeq,
		   s_iCode,
		   centre.x,
		   centre.y,
		   centre.z);
}

char residues_collide(Atom a, Atom b)
{

	double distance = sqrt(
		pow(a.centre.x - b.centre.x, 2) +
		pow(a.centre.y - b.centre.y, 2) +
		pow(a.centre.z - b.centre.z, 2));
	// printf("D: %f, a.x: %f, b.x: %f, a.y: %f, b.y: %f, a.z: %f, b.z: %f",
	// 	   distance,
	// 	   a.centre.x, b.centre.x,
	// 	   a.centre.y, b.centre.y,
	// 	   a.centre.z, b.centre.z);

	if (distance <= THRESHOLD)
	{
		return 1;
	}

	return 0;
}

int internal_collisions(int domainStart, int domainEnd, int *res)
{
	int internal = 0;
	for (int i = domainStart; i <= domainEnd; i++)
	{
		for (int j = i + 1; j <= domainEnd; j++)
		{

			if (residues_collide(atom[i], atom[j]))
			{
				//printf("index: %d\n", internal);
				res[i] = 1;
				res[j] = 1;
				internal++;
			}
		}
	}
	return internal;
}

int external_collisions(int numRes)
{
	int external = 0;
	for (int i = 0; i <= numRes; i++)
	{
		for (int j = 0; j <= numRes; j++)
		{
			if (resA[i] && resB[j] && residues_collide(atom[i], atom[j]) == 1)
			{
				external++;
			}
		}
	}
	return external;
}

double calculate_split(int numResidues, int splitIndex)
{
	memset(resA, 0, sizeof(resA));
	memset(resB, 0, sizeof(resB));

	int int_a = internal_collisions(0, splitIndex, resA);
	int int_b = internal_collisions(splitIndex + 1, numResidues, resB);
	int ext_ab = external_collisions(numResidues);

#if DEBUG
	printf("SPLIT: %d\n", splitIndex);
	printf("A: %d\n", int_a);
	printf("B: %d\n", int_b);
	printf("EXT: %d\n", ext_ab);
#endif
	if (int_a != 0 && int_b != 0)
	{
		double val = ((double)int_a / (double)ext_ab) * ((double)int_b / (double)ext_ab);
#if DEBUG
		printf("Val: %f\n\n", val);
#endif
		return val;
	}
	return 0;
}

void maximize_split(int numResidues)
{
	int maxIndex = 0;
	double maxValue = 0;
	for (int i = 0; i <= numResidues; i++)
	{
		double val = calculate_split(numResidues, i);
		if (val > maxValue)
		{
			maxValue = val;
			maxIndex = i;
		}
#if GRAPH
		printf("%d %f\n", i, val);
#endif
	}
#if DEBUG || RESULT
	printf("Split Residue: %d\n", maxIndex);
#endif
}

int main(int argc, char **argv)
{
	int numResidues;
	//int i, j;

	if (argc < 2)
	{
		(void)fprintf(stderr, "usage: residue_array file.pdb\n");
		exit(0);
	}

	numResidues = read_data(argv[1]);
	maximize_split(numResidues);

	// for (int i = 1; i <= numResidues; ++i)
	// {
	// 	for (int j = 1; j <= residue[i].numAtoms; ++j)
	// 	{
	// 		write_pdb_atom(
	// 			residue[i].atom[j].serial,
	// 			residue[i].atom[j].atomName,
	// 			residue[i].atom[j].altLoc,
	// 			residue[i].resName,
	// 			residue[i].chainID,
	// 			residue[i].resSeq,
	// 			residue[i].iCode,
	// 			residue[i].atom[j].centre);
	// 	}
	// }

	return 0;
}
