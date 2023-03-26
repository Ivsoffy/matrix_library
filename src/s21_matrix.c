#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int error_code = 0;
  double **matrix = NULL;
  if (rows < 1 || columns < 1)
    error_code = 1;
  else if ((matrix = (double **)malloc(rows * sizeof(double *) +
                                       rows * columns * sizeof(double))) ==
           NULL)
    error_code = 1;
  else {
    matrix[0] = (double *)(matrix + rows);
    for (int i = 1; i < rows; ++i) matrix[i] = matrix[0] + i * columns;
    for (int r = 0; r < rows; ++r)
      for (int c = 0; c < columns; ++c) matrix[r][c] = 0;
    result->rows = rows;
    result->columns = columns;
    result->matrix = matrix;
  }
  return error_code;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int error_code = FAILURE, flag = 0;
  if (A && B && A->columns == B->columns && A->rows == B->rows) {
    error_code = SUCCESS;
    for (int i = 0; i < A->rows; ++i)
      for (int j = 0; j < A->columns; ++j)
        if (A->matrix[i][j] != B->matrix[i][j]) flag = 1;
    if (flag) error_code = FAILURE;
  }
  return error_code;
}

int is_null(matrix_t A) { return !(A.rows && A.columns); }

int is_different(matrix_t A, matrix_t B) {
  return (A.rows != B.rows) || (A.columns != B.columns);
}

void s21_remove_matrix(matrix_t *A) {
  if (A->matrix) {
    free(A->matrix);
    A->matrix = NULL;
  }
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error_code = 0;
  if ((A->rows == B->rows) && (A->columns == B->columns)) {
    if (s21_create_matrix(A->rows, B->columns, result))
      error_code = 1;
    else {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = (A->matrix[i][j] + B->matrix[i][j]);
        }
      }
    }
  } else {
    error_code = 2;
  }
  return error_code;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error_code = 0;
  if ((A->rows == B->rows) && (A->columns == B->columns)) {
    if (s21_create_matrix(A->rows, B->columns, result))
      error_code = 1;
    else {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = (A->matrix[i][j] - B->matrix[i][j]);
        }
      }
    }
  } else {
    error_code = 2;
  }
  return error_code;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int error_code = 0;
  if (s21_create_matrix(A->rows, A->columns, result))
    error_code = 1;
  else {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  }
  return error_code;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error_code = 0;
  if ((A == NULL) || (B == NULL) || (A->rows <= 0) || A->columns <= 0 ||
      (B->rows <= 0) || B->columns <= 0)
    error_code = 1;
  else if (A->columns != B->rows)
    error_code = 2;
  else {
    if (s21_create_matrix(A->rows, B->columns, result))
      error_code = 1;
    else {
      for (int i = 0; i < A->columns; i++)
        for (int j = 0; j < A->rows; j++)
          for (int k = 0; k < B->columns; k++)
            result->matrix[j][k] += A->matrix[j][i] * B->matrix[i][k];
    }
  }
  return error_code;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int error_code = 0;
  if (s21_create_matrix(A->columns, A->rows, result))
    error_code = 1;
  else {
    for (int i = 0; i < A->columns; ++i)
      for (int j = 0; j < A->rows; ++j) result->matrix[i][j] = A->matrix[j][i];
  }
  return error_code;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int error_code = 0;
  matrix_t minor;
  double addition;
  if (A->columns != A->rows)
    error_code = 2;
  else {
    if (s21_create_matrix(A->rows, A->columns, result) ||
        s21_create_matrix(A->rows - 1, A->columns - 1, &minor))
      error_code = 1;
    else {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          get_minor(i, j, *A, &minor);
          s21_determinant(&minor, &addition);
          addition = addition * (((i + j) % 2) ? -1 : 1);
          result->matrix[i][j] = addition;
        }
      }
      s21_remove_matrix(&minor);
    }
  }
  // output(minor);
  return error_code;
}

int s21_determinant(matrix_t *A, double *result) {
  int error_code = 0;
  double addition = 0, addition_r = 0;
  matrix_t minor;
  if ((A == NULL) || (A->matrix == NULL))
    // if (A == NULL)
    error_code = 1;
  else if (A->columns != A->rows)
    error_code = 2;
  else {
    if (A->rows == 2)
      *result =
          A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
    else if (A->rows == 1)
      *result = A->matrix[0][0];
    else {
      for (int j = 0; j < A->columns; j++) {
        if (s21_create_matrix(A->rows - 1, A->columns - 1, &minor))
          error_code = 1;
        else {
          get_minor(0, j, *A, &minor);
          s21_determinant(&minor, &addition);
          addition_r += A->matrix[0][j] * addition * ((j % 2) ? -1 : 1);
          s21_remove_matrix(&minor);
        }
      }
      *result = addition_r;
    }
  }
  return error_code;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int error_code = 0;
  if (A->rows != A->columns)
    error_code = 2;
  else {
    double det = 0;
    s21_determinant(A, &det);
    if (!det)
      error_code = 2;
    else {
      matrix_t minor, tran;
      s21_calc_complements(A, &minor);
      s21_transpose(&minor, &tran);
      s21_mult_number(&tran, 1 / det, result);
      s21_remove_matrix(&minor);
      s21_remove_matrix(&tran);
    }
  }
  return error_code;
}

void get_minor(int rows, int cols, matrix_t A, matrix_t *minor) {
  int m = 0, n = 0, fl = 0;
  for (int i = 0; i < A.rows; i++) {
    fl = 0;  // ok
    for (int j = 0; j < A.columns; j++) {
      if (i != rows && j != cols) {
        minor->matrix[m][n] = A.matrix[i][j];
        n++;
        fl++;
      }
    }
    if (fl) m++;
    n = 0;
  }
}