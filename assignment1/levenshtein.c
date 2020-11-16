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
  int score, tmp;
  char X[MAX_LENGTH + 1] = "ATCGAT";
  char Y[MAX_LENGTH + 1] = "ATACGT";

  int F[MAX_LENGTH + 1][MAX_LENGTH + 1];     /* score matrix */

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

  /*
   * Initialise matrices
   */

  F[0][0] = 0;
  for (i = 1; i <= m; i++) {
    F[i][0] = F[i - 1][0] + 1;
  }
  for (j = 1; j <= n; j++) {
    F[0][j] = F[0][j - 1] + 1;
  }

  /*
   * Fill matrices
   */

  for (i = 1; i <= m; i++) {

    for (j = 1; j <= n; j++) {

      if (X[i - 1] == Y[j - 1]) {
        score = F[i - 1][j - 1];
      } else {
        score = F[i - 1][j - 1] + 1;
      }

      tmp = F[i - 1][j] + 1;
      if (tmp < score) {
        score = tmp;
      }

      tmp = F[i][j - 1] + 1;
      if (tmp < score) {
        score = tmp;
      }

      F[i][j] = score;
    }
  }

  /*
   * Print score matrix
   */
  printf("----- levenshtein.c ----- \n");
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
   * Print the levenshtein distance
   */
  printf("%s%d", "Levenshtein distance: ", F[m][n]);
  return (1);
}
