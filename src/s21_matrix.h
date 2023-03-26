#include <stdio.h>
#include <stdlib.h>
// #include "s21_matrix_support.h"

#define SUCCESS 1
#define FAILURE 0
#define EPS 1e-7

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);
// helpful
void output(matrix_t A);
void input(matrix_t *A);
void get_minor(int rows, int cols, matrix_t A, matrix_t *minor);
void fill_matrix_by_list(matrix_t *matrix, double list[]);
void fill_matrix(matrix_t *matrix, double number);
int is_different(matrix_t A, matrix_t B);
int is_null(matrix_t A);