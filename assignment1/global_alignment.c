/*
 * To compile this C program, placing the executable file in 'global', type:
 *
 *      gcc -o global global_alignment.c
 *
 * To run the program, type:
 *
 *      ./global
 */

#include <stdio.h>

#define MAX_LENGTH 100

#define MATCH_SCORE 2
#define MISMATCH_SCORE -1
#define GAP_PENALTY 2

#define STOP 0
#define UP 1
#define LEFT 2
#define DIAG 3

#define percent_identity(mChars, pNom) (mChars * 100 / pNom)

int main() {
  int i, j;
  int m, n;
  int barIndex, matchingChars = 0, percentNom, hamming = 0;
  int alignmentLength, score, tmp;
  int maxLength, minLength;
  int branches = 0;
  char X[MAX_LENGTH + 1] = "ATCGAT";
  char Y[MAX_LENGTH + 1] = "ATACGT";
  char bars[MAX_LENGTH + 1];

  int F[MAX_LENGTH + 1][MAX_LENGTH + 1];     /* score matrix */
  int trace[MAX_LENGTH + 1][MAX_LENGTH + 1]; /* trace matrix */
  char alignX[MAX_LENGTH * 2];               /* aligned X sequence */
  char alignY[MAX_LENGTH * 2];               /* aligned Y sequence */

  /*
   * Find lengths of (null-terminated) strings X and Y
   */
  m = 0;
  n = 0;
  while (X[m] != 0) {
    m++;
  }
  while (Y[n] != 0) {
    n++;
  }

  if(m > n) {
    maxLength = m;
    minLength = n;
  } else {
    minLength = m;
    maxLength = n;
  }
  barIndex = maxLength;
  // I chose to use the method of calculating the percent where you chose
  // the length of the shortest string as the denominator
  percentNom = minLength;

  /*
   * Initialise matrices
   */

  F[0][0] = 0;
  trace[0][0] = STOP;
  for (i = 1; i <= m; i++) {
    F[i][0] = F[i - 1][0] - GAP_PENALTY;
    trace[i][0] = STOP;
  }
  for (j = 1; j <= n; j++) {
    F[0][j] = F[0][j - 1] - GAP_PENALTY;
    trace[0][j] = STOP;
  }

  /*
   * Fill matrices
   */

  for (i = 1; i <= m; i++) {

    for (j = 1; j <= n; j++) {

      if (X[i - 1] == Y[j - 1]) {
        score = F[i - 1][j - 1] + MATCH_SCORE;
      } else {
        score = F[i - 1][j - 1] + MISMATCH_SCORE;
      }
      trace[i][j] = DIAG;

      tmp = F[i - 1][j] - GAP_PENALTY;
      if (tmp > score) {
        score = tmp;
        trace[i][j] = UP;
      } else if(tmp == score) {
        branches++;
      }


      tmp = F[i][j - 1] - GAP_PENALTY;
      if (tmp > score) {
        score = tmp;
        trace[i][j] = LEFT;
      } else if(tmp == score) {
        branches++;
      }

      F[i][j] = score;
    }
  }

  /*
   * Print score matrix
   */
  printf("------ global_alignment.c -----\n");
  printf("Score matrix:\n      ");
  for (j = 0; j < n; ++j) {
    printf("%5c", Y[j]);
  }
  printf("\n");
  for (i = 0; i <= m; i++) {
    if (i == 0) {
      printf(" ");
    } else {
      printf("%c", X[i - 1]);
    }
    for (j = 0; j <= n; j++) {
      printf("%5d", trace[i][j]);
    }
    printf("\n");
  }
  printf("\n");

  /*
   * Trace back from the lower-right corner of the matrix
   */

  i = m;
  j = n;
  alignmentLength = 0;

  while (trace[i][j] != STOP) {

    switch (trace[i][j]) {

    case DIAG:
      alignX[alignmentLength] = X[i - 1];
      alignY[alignmentLength] = Y[j - 1];
      i--;
      j--;
      bars[barIndex] = '|';
      matchingChars++;
      alignmentLength++;
      barIndex--;
      break;

    case LEFT:
      alignX[alignmentLength] = '-';
      alignY[alignmentLength] = Y[j - 1];
      j--;
      bars[barIndex] = ' ';
      alignmentLength++;
      barIndex--;
      break;

    case UP:
      alignX[alignmentLength] = X[i - 1];
      alignY[alignmentLength] = '-';
      i--;
      bars[barIndex] = ' ';
      alignmentLength++;
      barIndex--;
      break;
    }
  }

  /*
   * Unaligned beginning
   */

  while (i > 0) {
    alignX[alignmentLength] = X[i - 1];
    alignY[alignmentLength] = '-';
    if (barIndex >= 0) {
      bars[barIndex] = ' ';
      barIndex--;
    }
    i--;
    alignmentLength++;
  }

  while (j > 0) {
    alignX[alignmentLength] = '-';
    alignY[alignmentLength] = Y[j - 1];
    if (barIndex >= 0) {
      bars[barIndex] = ' ';
      barIndex--;
    }
    j--;
    alignmentLength++;
  }

  int p_identity = percent_identity(matchingChars, percentNom);

  for ( int index = 0; index < minLength; index++) {
    if(X[index] == Y[index]) {
      hamming++;
    }  
  }
  /*
   * Print alignment
   */

  for (i = alignmentLength - 1; i >= 0; i--) {
    printf("%c", alignX[i]);
  }
  printf("\n");
  for (i = 0; i < alignmentLength; i++) {
   printf("%c", bars[i]);
  }
  printf("\n");
  for (i = alignmentLength - 1; i >= 0; i--) {
    printf("%c", alignY[i]);
  }
  printf("\n");
  printf("%d%s\n", matchingChars, " matching characters");
  printf("%d%s\n", p_identity, "% percent identity");
  printf("%d%s\n", hamming, " hamming distance");
  //printf("%d%s", branches / 2, " different optimal paths");
  return (1);
}
