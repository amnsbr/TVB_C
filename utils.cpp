#include <iostream>
#include <cmath>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics_double.h>


// Helpers
void gsl_matrix_print(FILE * f, gsl_matrix * m) {
    int i, j;
    for (i = 0; i < (m->size1); i++) {
        for (j = 0; j < (m->size2); j++) {
            fprintf(f, "%f", gsl_matrix_get(m, i, j));
            if (j != (m->size2 - 1)) {
                fprintf(f, " ");
            }
        }
        fprintf(f, "\n");
    }
}

void gsl_matrix_print(FILE * f, gsl_matrix_complex * m) {
    int i, j;
    for (i = 0; i < (m->size1); i++) {
        for (j = 0; j < (m->size2); j++) {
            fprintf(f, "%f+%fj", 
                GSL_REAL(gsl_matrix_complex_get(m, i, j)),
                GSL_IMAG(gsl_matrix_complex_get(m, i, j)));
            if (j != (m->size2 - 1)) {
                fprintf(f, " ");
            }
        }
        fprintf(f, "\n");
    }
}

void gsl_matrix_print(gsl_matrix * m) {
    gsl_matrix_print(stdout, m);
}

void gsl_matrix_print(gsl_matrix_complex * m) {
    gsl_matrix_print(stdout, m);
}

gsl_vector * repeat(double a, int size) {
    gsl_vector * out = gsl_vector_alloc(size);
    gsl_vector_set_all(out, a);
    return out;
}

gsl_vector * vector_scale(gsl_vector * v, double a) {
    gsl_vector * out = gsl_vector_alloc(v->size);
    gsl_vector_memcpy(out, v);
    gsl_vector_scale(out, a);
    return out;
}

gsl_matrix * mul_eye(double a, int size) {
    gsl_matrix * out = gsl_matrix_alloc(size, size);
    gsl_matrix_set_identity(out);
    gsl_matrix_scale(out, a);
    return out;
}

gsl_matrix * zeros(int size1, int size2) {
    gsl_matrix * out = gsl_matrix_alloc(size1, size2);
    gsl_matrix_set_zero(out);
    return out;
}

gsl_matrix * zeros(int size) {
    return zeros(size, size);
}

gsl_matrix * make_diag(gsl_vector * v) {
    int size = v->size;
    gsl_vector * ones = gsl_vector_alloc(size);
    gsl_vector_set_all(ones, 1.0);
    gsl_matrix * out = gsl_matrix_alloc(size, size);
    gsl_blas_dger(1.0, ones, v, out); // outer product of 1 * v
    gsl_matrix * eye = gsl_matrix_alloc(size, size);
    gsl_matrix_set_identity(eye);
    gsl_matrix_mul_elements(out, eye);
    return out;
}

gsl_matrix * dot(gsl_matrix * A, gsl_matrix * B) {
    int size1 = A->size1;
    int size2 = A->size2;
    gsl_matrix * out = gsl_matrix_alloc(size1, size2);
    gsl_blas_dgemm(
        CblasNoTrans, CblasNoTrans,
        1.0, A, B, 0.0, out);
    return out;
}

gsl_vector * dot(gsl_matrix * m, gsl_vector * v) {
    gsl_vector * out = gsl_vector_alloc(m->size1);
    gsl_blas_dgemv(
        CblasNoTrans,
        1.0, m, v, 0.0, out);
    return out;
}

gsl_matrix_complex * dot(gsl_matrix_complex * A, gsl_matrix_complex * B) {
    int size1 = A->size1;
    int size2 = B->size2;
    gsl_complex alpha;
    GSL_SET_COMPLEX(&alpha, 1.0, 0);
    gsl_complex beta;
    GSL_SET_COMPLEX(&beta, 0.0, 0);
    gsl_matrix_complex * out = gsl_matrix_complex_alloc(size1, size2);
    gsl_blas_zgemm(
        CblasNoTrans, CblasNoTrans,
        alpha, A, B, beta, out);
    return out;
}

gsl_matrix * concat(std::vector<gsl_matrix*> submatrices, int axis) {
    /* Concatenates nested list of gsl_matrix objects that have an equal shape */
    // TODO: remove the constraint of equal shape, only the concatenation axis
    // must be equal
    int size1 = 0; int size2 = 0;
    int n = submatrices.size();
    if (axis == 0) {
        size2 = submatrices[0]->size2;
        for (gsl_matrix * submatrix: submatrices) {
            size1 += submatrix->size1;
            if (submatrix->size2 != size2) {
                std::cout << "\nInput matrices must have the same col size\n";
            }
        }
        gsl_matrix * out = gsl_matrix_alloc(size1, size2);
        int corner_i = 0;
        for (int i=0; i<n; i++) {
            gsl_matrix_view curr_mat_view = gsl_matrix_submatrix(out, corner_i, 0, submatrices[i]->size1, size2);
            gsl_matrix_memcpy(&curr_mat_view.matrix, submatrices[i]);
            corner_i += submatrices[i]->size1;
        }
        return out;
    } else {
        size1 = submatrices[0]->size1;
        for (gsl_matrix * submatrix: submatrices) {
            size2 += submatrix->size2;
            if (submatrix->size1 != size1) {
                std::cout << "\nInput matrices must have the same row size\n";
            }
        }
        gsl_matrix * out = gsl_matrix_alloc(size1, size2);
        int corner_j = 0;
        for (int j=0; j<n; j++) {
            gsl_matrix_view curr_mat_view = gsl_matrix_submatrix(out, 0, corner_j, size1, submatrices[j]->size2);
            gsl_matrix_memcpy(&curr_mat_view.matrix, submatrices[j]);
            corner_j += submatrices[j]->size2;
        }
        return out;
    }
}