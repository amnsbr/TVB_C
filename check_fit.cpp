#include <iostream>
#include <fstream>
#include <string>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
/*
wget gsl-<...>
tar -xvzf gsl-<...>.tar.gz
cd gsl-<...>
./configure --prefix=<dir>/gsl_build
make && make install
g++ -o check_fit check_fit.cpp ../gsl_build/lib/libgsl.a ../gsl_build/lib/libgslcblas.a -I/data/project/ei_development/tools/gsl_build/include -lm
*/
extern "C" {
#include "ks.h"
};
/*
git clone https://github.com/agentlans/libks.git
cd libks
make
cd ../TVB_C
g++ -o check_fit check_fit.cpp ../libks/libks.so ../gsl_build/lib/libgsl.a ../gsl_build/lib/libgslcblas.a -I/data/project/ei_development/tools/gsl_build/include -I /data/project/ei_development/tools/libks/include -lm && ./check_fit
*/


gsl_vector * fc_tril(gsl_matrix * bold, const int n_pairs) {
    /*
     Given empirical/simulated bold (n_vols x n_regions) returns
     the lower triangle of the FC
    */
    int i, j;
    double corr;
    gsl_vector * FC_tril = gsl_vector_alloc(n_pairs);
    int curr_idx = 0;
    for (i = 0; i<(bold->size2); i++) {
        for (j = 0; j<(bold->size2); j++) {
            if (i > j) {
                gsl_vector_view ts_i = gsl_matrix_column(bold, i);
                gsl_vector_view ts_j = gsl_matrix_column(bold, j);
                corr = gsl_stats_correlation(
                    ts_i.vector.data, ts_i.vector.stride,
                    ts_j.vector.data, ts_j.vector.stride,
                    (bold->size1)
                );
                gsl_vector_set(FC_tril, curr_idx, corr);
                curr_idx ++;
            }
        }
    }
    return FC_tril;
}

gsl_vector * fcd_tril(gsl_matrix * bold, const int step, const int window_size) {
    /*
     Calculates the functional connectivity dynamics matrix (lower triangle)
     given BOLD, step and window size. Note that the actual window size is +1 higher.
     The FCD matrix shows similarity of FC patterns between the windows.
    */
    const int n_vols = bold->size1;
    const int n_regions = bold->size2;
    const int n_pairs = ((n_regions) * (n_regions - 1)) / 2;

    int curr_center = 0;
    // too lazy to calculate the exact number of windows, but this should be the maximum bound (not sure though)
    gsl_matrix * window_FC_trils = gsl_matrix_alloc(n_pairs, ((n_vols + window_size) / step));
    // gsl_matrix_set_all(window_FC_trils, -1.5);
    // calculate the FC of each window
    int n_windows = 0;
    while (curr_center <= n_vols) {
        int start = curr_center - (window_size/2);
        if (start < 0)
            start = 0;
        int end = curr_center + (window_size/2);
        if (end >= n_vols)
            end = n_vols-1;
        gsl_matrix_view bold_window =  gsl_matrix_submatrix(bold, start, 0, end-start+1, n_regions);
        gsl_vector * window_FC_tril = fc_tril(&bold_window.matrix, n_pairs);
        gsl_matrix_set_col(window_FC_trils, n_windows, window_FC_tril);
        curr_center += step;
        n_windows ++;
    }
    // calculate the FCD matrix (lower triangle)
    int window_i, window_j;
    double corr;
    int n_window_pairs = ((n_windows) * (n_windows - 1)) / 2;
    gsl_vector * FCD_tril = gsl_vector_alloc(n_window_pairs);
    int curr_idx = 0;
    for (window_i=0; window_i<n_windows; window_i++) {
        for (window_j=0; window_j<n_windows; window_j++) {
            if (window_i > window_j) {
                gsl_vector_view FC_i = gsl_matrix_column(window_FC_trils, window_i);
                gsl_vector_view FC_j = gsl_matrix_column(window_FC_trils, window_j);
                // gsl_vector_fprintf(stdout, &FC_i.vector, "%f");
                corr = gsl_stats_correlation(
                    FC_i.vector.data, FC_i.vector.stride,
                    FC_j.vector.data, FC_j.vector.stride,
                    n_pairs
                );
                gsl_vector_set(FCD_tril, curr_idx, corr);
                curr_idx ++;
            }
        }
    }
    return FCD_tril;
}

int main() {
    const int n_regions = 84;
    const int n_pairs = ((n_regions) * (n_regions - 1)) / 2;
    // TODO: allow different duration, TR and n_vols for empirical and simualted data
    const int duration = 100000;
    const int TR = 720;
    const int n_vols = duration / TR;
    int FCD_step = 5;
    int FCD_window_size = 10; // not including the center
    const char* emp_FC_path = "output/emp_FC_tril.txt";
    const bool write_sim_FC = true;
    // Read or create the empirical FC (lower triangle)
    gsl_vector * emp_FC_tril = gsl_vector_alloc(n_pairs);
    // if (FILE * old_emp_FC_file = fopen(emp_FC_path, "r")) {
    //     gsl_vector_fscanf(old_emp_FC_file, emp_FC_tril);
    // } else {
        gsl_matrix * emp_bold = gsl_matrix_alloc(n_vols, n_regions);
        FILE * emp_bold_file = fopen("output/pseudo_empBOLD.txt", "r");
        gsl_matrix_fscanf(emp_bold_file, emp_bold);
        emp_FC_tril = fc_tril(emp_bold, n_pairs);
        fclose(emp_bold_file);
        FILE * new_emp_FC_file = fopen(emp_FC_path, "w");
        gsl_vector_fprintf(new_emp_FC_file, emp_FC_tril, "%f");
        fclose(new_emp_FC_file);
    // }
    // Read simulated bold and create simulated FC (lower triangle)
    gsl_matrix * sim_bold = gsl_matrix_alloc(n_vols, n_regions);
    FILE * sim_bold_file = fopen("output/simBOLD_param_set_1.txt", "r");
    gsl_matrix_fscanf(sim_bold_file, sim_bold);
    fclose(sim_bold_file);
    gsl_vector * sim_FC_tril = fc_tril(sim_bold, n_pairs);
    if (write_sim_FC) {
        FILE * sim_FC_file = fopen("output/sim_FC_tril.txt", "w");
        gsl_vector_fprintf(sim_FC_file, sim_FC_tril, "%f");
        fclose(sim_FC_file);
    }
    double emp_sim_FC_corr = gsl_stats_correlation(
        emp_FC_tril->data, emp_FC_tril->stride,
        sim_FC_tril->data, sim_FC_tril->stride,
        emp_FC_tril->size
    );
    std::cout << "FC correlation: " << emp_sim_FC_corr << "\n";

    // Calculate FCD
    gsl_vector * sim_FCD_tril = fcd_tril(sim_bold, FCD_step, FCD_window_size);
    FILE * sim_FCD_file = fopen("output/sim_FCD_tril_param_set_1.txt", "w");
    gsl_vector_fprintf(sim_FCD_file, sim_FCD_tril, "%f");
    fclose(sim_FCD_file);
    gsl_vector * emp_FCD_tril = fcd_tril(emp_bold, FCD_step, FCD_window_size);
    double ks = ks_stat(emp_FCD_tril->data, sim_FCD_tril->data, emp_FCD_tril->size, sim_FCD_tril->size, NULL);
    std::cout << "FCD KS: " << ks << "\n";
    return 0;
}