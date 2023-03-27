#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  resulting_code res = OK;
  double **matrix = NULL;
  if (rows < 1 || columns < 1) {
    res = INCORRECT_MATRIX;
    result->matrix = NULL;
  } else if ((matrix = (double **)malloc(rows * sizeof(double *) +
                                         rows * columns * sizeof(double))) ==
             NULL)
    res = INCORRECT_MATRIX;
  else {
    matrix[0] = (double *)(matrix + rows);
    for (int i = 1; i < rows; ++i) matrix[i] = matrix[0] + i * columns;
    for (int r = 0; r < rows; ++r)
      for (int c = 0; c < columns; ++c) matrix[r][c] = 0;
    result->rows = rows;
    result->columns = columns;
    result->matrix = matrix;
  }
  return res;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int res = SUCCESS;
  if (is_valid(A) && is_valid(B) && A->rows == B->rows &&
      A->columns == B->columns) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        if ((fabs(A->matrix[i][j] - B->matrix[i][j])) >= 1e-7) {
          res = FAILURE;
          break;
        }
      }
    }
  } else
    res = FAILURE;
  return res;
}

int is_valid(matrix_t *matrix) {
  int res = 1;
  if (!matrix || !matrix->matrix || matrix->rows <= 0 || matrix->columns <= 0)
    res = 0;
  return res;
}

void s21_remove_matrix(matrix_t *A) {
  if (A->matrix) {
    free(A->matrix);
    A->matrix = NULL;
  }
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  resulting_code res = OK;
  if ((A->rows == B->rows) && (A->columns == B->columns)) {
    if (s21_create_matrix(A->rows, B->columns, result))
      res = INCORRECT_MATRIX;
    else {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = (A->matrix[i][j] + B->matrix[i][j]);
        }
      }
    }
  } else {
    res = CALCULATION_ERROR;
  }
  return res;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  resulting_code res = OK;
  if ((A->rows == B->rows) && (A->columns == B->columns)) {
    if (s21_create_matrix(A->rows, B->columns, result))
      res = INCORRECT_MATRIX;
    else {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = (A->matrix[i][j] - B->matrix[i][j]);
        }
      }
    }
  } else {
    res = CALCULATION_ERROR;
  }
  return res;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  resulting_code res = OK;
  if (s21_create_matrix(A->rows, A->columns, result))
    res = INCORRECT_MATRIX;
  else {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  }
  return res;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  resulting_code res = OK;
  if ((A == NULL) || (B == NULL) || (A->rows <= 0) || A->columns <= 0 ||
      (B->rows <= 0) || B->columns <= 0)
    res = INCORRECT_MATRIX;
  else if (A->columns != B->rows)
    res = CALCULATION_ERROR;
  else {
    if (s21_create_matrix(A->rows, B->columns, result))
      res = INCORRECT_MATRIX;
    else {
      for (int i = 0; i < A->columns; i++)
        for (int j = 0; j < A->rows; j++)
          for (int k = 0; k < B->columns; k++)
            result->matrix[j][k] += A->matrix[j][i] * B->matrix[i][k];
    }
  }
  return res;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  resulting_code res = OK;
  if (s21_create_matrix(A->columns, A->rows, result))
    res = INCORRECT_MATRIX;
  else {
    for (int i = 0; i < A->columns; ++i)
      for (int j = 0; j < A->rows; ++j) result->matrix[i][j] = A->matrix[j][i];
  }
  return res;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  resulting_code res = OK;
  matrix_t minor;
  double addition;
  if (A->columns != A->rows)
    res = CALCULATION_ERROR;
  else {
    if (s21_create_matrix(A->rows, A->columns, result) ||
        s21_create_matrix(A->rows - 1, A->columns - 1, &minor))
      res = INCORRECT_MATRIX;
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
  return res;
}

int s21_determinant(matrix_t *A, double *result) {
  resulting_code res = OK;
  double addition = 0, addition_r = 0;
  matrix_t minor;
  if ((A == NULL) || (A->matrix == NULL))
    res = INCORRECT_MATRIX;
  else if (A->columns != A->rows)
    res = CALCULATION_ERROR;
  else {
    if (A->rows == 2)
      *result =
          A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
    else if (A->rows == 1)
      *result = A->matrix[0][0];
    else {
      for (int j = 0; j < A->columns; j++) {
        if (s21_create_matrix(A->rows - 1, A->columns - 1, &minor))
          res = INCORRECT_MATRIX;
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
  return res;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  resulting_code res = OK;
  double det = 0.0;
  res = s21_determinant(A, &det);
  if (!res) {
    if (fabs(det) > 1e-7) {
      matrix_t minor = {0};
      res = s21_calc_complements(A, &minor);
      if (!res) {
        matrix_t transp_minor = {0};
        res = s21_transpose(&minor, &transp_minor);
        if (!res) {
          if (!res) res = s21_mult_number(&transp_minor, 1 / det, result);
          s21_remove_matrix(&transp_minor);
        }
        s21_remove_matrix(&minor);
      }
    } else
      res = CALCULATION_ERROR;
  }
  return res;
}

void get_minor(int rows, int cols, matrix_t A, matrix_t *minor) {
  int m = 0, n = 0, fl = 0;
  for (int i = 0; i < A.rows; i++) {
    fl = 0;
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
