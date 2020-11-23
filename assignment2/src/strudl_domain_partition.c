/*
 * File:	strudl_domain_partition.c
 * Author: 	Carl Holmberg
 * Purpose:	A failed attempt at exercise 3
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_ATOMS 10000
#define LINE_LENGTH 81
#define THRESHOLD 7

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
	int flag;
	Point centre;
} Atom;

/*
 * Declare an array to hold data read from the ATOM records of a PDB file.
 */
Atom atom[MAX_ATOMS + 1];

int resA[MAX_ATOMS + 1];
int resB[MAX_ATOMS + 1];

int results[MAX_ATOMS + 1];
int saves[MAX_ATOMS + 1][MAX_ATOMS + 1];

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
				atom[i].flag = 0;
			}
		}
	}
	return i;
}

void print_result(int numResidues)
{
	int i = 0, j = 0;
	printf("   A    B\n");
	while (i <= numResidues || j <= numResidues)
	{
		while (resA[i] == 0 && i <= numResidues)
		{
			i++;
		}
		while (resB[j] == 0 && j <= numResidues)
		{
			j++;
		}
		int a = atom[i].resSeq;
		int b = atom[j].resSeq;
		//printf("B: %d %d\n", j, b);

		printf("%4d %4d\n", a, b);
		i++;
		j++;
	}
}

int residues_collide(Atom a, Atom b)
{

	double distance = sqrt(
		pow(a.centre.x - b.centre.x, 2) +
		pow(a.centre.y - b.centre.y, 2) +
		pow(a.centre.z - b.centre.z, 2));

	if (distance <= THRESHOLD)
	{
		return 1;
	}

	return 0;
}

void move(int index, int initial)
{
	if (resA[index] == 1)
	{
		resA[index] = 0;
		resB[index] = 1;
		//printf("index %d moved to B\n", index);
	}
	else if (resB[index] == 1)
	{
		resB[index] = 0;
		resA[index] = 1;
		//  Flag each residue leaving B
		if (initial == 0)
		{
			atom[index].flag = 1;
		}
		//printf("index %d moved to A\n", index);
	}
}

int find_min_a(int numResidues)
{
	int minIndex = -1;
	int minValue = sizeof(int);
	for (int i = 0; i <= numResidues; i++)
	{
		if (resA[i] == 1 && atom[i].flag == 0)
		{
			int collisions = 0;
			for (int j = 0; j <= numResidues; j++)
			{
				//printf("I: %d [%d] | J: %d [%d]\n", i, resA[i], j, resB[j]);
				if (resA[j] == 1 && atom[i].flag == 0)
				{
					if (residues_collide(atom[i], atom[j]))
					{
						collisions++;
					}
				}
			}
			if (collisions < minValue)
			{
				minIndex = i;
				minValue = collisions;
			}
		}
	}
	return minIndex;
}

int find_min_b(int numResidues)
{
	int minIndex = -1;
	int minValue = 100000000;
	for (int i = 0; i <= numResidues; i++)
	{
		if (resB[i] == 1)
		{
			int collisions = 0;
			for (int j = 0; j <= numResidues; j++)
			{
				//printf("I: %d [%d] | J: %d [%d]\n", i, resA[i], j, resB[j]);
				if (resB[j] == 1)
				{
					if (residues_collide(atom[i], atom[j]))
					{
						collisions++;
					}
				}
			}
			if (collisions < minValue)
			{
				minIndex = i;
				minValue = collisions;
			}
		}
	}
	return minIndex;
}

int check_flagged(int numResidues)
{
	for (int i = 0; i <= numResidues; i++)
	{

		if (resA[i] && atom[i].flag == 0)
		{
			return 1;
		}
	}
	return 0;
}

int contact_area(int numResidues)
{
	//printf("start\n");
	int coll = 0;
	for (int i = 0; i <= numResidues; i++)
	{
		if (resA[i] == 1)
		{
			for (int j = 0; j <= numResidues; j++)
			{
				if (resB[j] == 1)
				{
					if (residues_collide(atom[i], atom[j]))
					{
						//printf("coll: %d | i: %d | j: %d\n", coll, i, j);
						coll++;
					}
				}
			}
		}
	}
	//print_result(numResidues);
	//printf("end: %d", coll);
	return coll;
}

double calculate_partitions(int k, int numResidues)
{
	memset(resA, 0, sizeof(resA));
	memset(resB, 0, sizeof(resB));

	for (int i = 0; i < numResidues; i++)
	{
		resB[i] = 1;
	}

	int a = 0, b = 0;
	// for (int i = 0; i < numResidues; i++)
	// {
	// 	printf("%d %d\n", i, resB[i]);
	// }
	// Randomly select a as seed for A
	a = rand() % numResidues + 1;
	//printf("RAND: %d\n", a);
	move(a, 1);

	// Increase size of A by adding residues til |A| = k
	int i = 0;
	while (i++ <= k)
	{
		b = find_min_b(numResidues);
		move(b, 1);
	}
	int counter = 0;
	results[counter] = contact_area(numResidues);
	for (int i = 0; i <= numResidues; i++)
	{
		saves[counter][i] = resA[i];
	}
	counter++;
	//print_result(numResidues);
	// While not all residues in A are flagged
	while (check_flagged(numResidues) == 1)
	{
		// 	Choose best a and b for exchange (keeping contact area as small as possible)
		a = find_min_a(numResidues);
		b = find_min_b(numResidues);
		if (a == -1 || b == -1)
		{
			break;
		}
		//printf("min_a: %d | min_b: %d\n", a, b);
		// 	Make exchange
		move(a, 0);
		move(b, 0);
		// 	Record contact area between A and B at each step
		// TODO: implement this
		results[counter] = contact_area(numResidues);
		for (int i = 0; i <= numResidues; i++)
		{
			saves[counter][i] = resA[i];
		}
		counter++;
	}
	//print_result(numResidues);

	int minValue = 1000000000;
	int minIndex = 0;
	for (int i = 0; i < counter; i++)
	{
		if (results[i] < minValue)
		{
			minValue = results[i];
			minIndex = i;
		}
	}

	for (int i = 0; i <= numResidues; i++)
	{
		if (saves[minIndex][i] == 1)
		{
			resA[i] = 1;
			resB[i] = 0;
		}
		else
		{
			resB[i] = 1;
			resA[i] = 0;
		}
	}
	//print_result(numResidues);
	// Return contact area.
	return contact_area(numResidues);
}

int main(int argc, char **argv)
{
	int numResidues;
	//int i, j;

	time_t t;

	/* Intializes random number generator */
	srand((unsigned)time(&t));

	if (argc < 2)
	{
		(void)fprintf(stderr, "usage: residue_array file.pdb\n");
		exit(0);
	}
	numResidues = read_data(argv[1]);

	int minIndex = 0;
	int minValue = 100000000;
	for (int i = 1; i <= numResidues / 2; i++)
	{
		int val = calculate_partitions(i, numResidues);
		printf("k [%d] | val [%d]\n", i, val);
		//printf("%d %d\n", i, val);
		if (val < minValue)
		{
			minValue = val;
			minIndex = i;
		}
	}
	printf("\n\nIndex: %d | Value: %d\n", minIndex, minValue);

	calculate_partitions(minIndex, numResidues);
	print_result(numResidues);
	//int val = calculate_partitions(50, numResidues);
	//printf("VAL: %d", val);
	return 0;
}
