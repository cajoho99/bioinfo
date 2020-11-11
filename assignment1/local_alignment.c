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
  int barIndex;
  int alignmentLength, score, tmp;
  int maxLength, minLength;
  int maxX = 0, maxY = 0, maxValue = 0;
  char X[MAX_LENGTH + 1] = "PAWHEAE";
  char Y[MAX_LENGTH + 1] = "HDAGAWGHEQ";
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

  /*
   * Initialise matrices
   */

  F[0][0] = 0;
  trace[0][0] = STOP;
  for (i = 1; i <= m; i++) {
    F[i][0] = 0;
    trace[i][0] = STOP;
  }
  for (j = 1; j <= n; j++) {
    F[0][j] = 0;
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
      }

      tmp = F[i][j - 1] - GAP_PENALTY;
      if (tmp > score) {
        score = tmp;
        trace[i][j] = LEFT;
      }
      if(score < 0) {
        score = 0;
        trace[i][j] = STOP;
      }
      if(score > maxValue) {
        maxValue = score;
        maxX = i;
        maxY = j;
      }

      F[i][j] = score;
    }
  }

  /*
   * Print score matrix
   */

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
      printf("%5d", F[i][j]);
    }
    printf("\n");
  }
  printf("\n");

  /*
   * Trace back from the lower-right corner of the matrix
   */

  i = maxX;
  j = maxY;
  alignmentLength = 0;
  if(i > j) {
    barIndex = i;
  } else {
    barIndex = j;
  }
  while (trace[i][j] != STOP) {

    switch (trace[i][j]) {

    case DIAG:
      alignX[alignmentLength] = X[i - 1];
      alignY[alignmentLength] = Y[j - 1];
      i--;
      j--;
      bars[barIndex] = '|';
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
    }
  }

  /*
   * Print alignment
   */

  for (i = alignmentLength - 1; i >= 0; i--) {
    printf("%c", alignX[i]);
  }
  printf("\n");
  for (i = barIndex + 1; i < barIndex + alignmentLength + 1; i++) {
    printf("%c", bars[i]);
  }
  printf("\n");
  for (i = alignmentLength - 1; i >= 0; i--) {
    printf("%c", alignY[i]);
  }
  printf("\n");
  return (1);
}
