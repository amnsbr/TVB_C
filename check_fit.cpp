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

gsl_vector * fcd_tril(gsl_matrix * bold, const int step, const int window_size, bool drop_edges) {
    /*
     Calculates the functional connectivity dynamics matrix (lower triangle)
     given BOLD, step and window size. Note that the actual window size is +1 higher.
     The FCD matrix shows similarity of FC patterns between the windows.
    */
    const int n_vols = bold->size1;
    const int n_regions = bold->size2;
    const int n_pairs = ((n_regions) * (n_regions - 1)) / 2;
    int first_center, last_center, window_center, window_start, window_end;
    if (drop_edges) {
        first_center = window_size / 2;
        last_center = n_vols - 1 - (window_size / 2);
    } else {
        first_center = 0;
        last_center = n_vols - 1;
    }

    // too lazy to calculate the exact number of windows, but this should be the maximum bound (not sure though)
    gsl_matrix * window_FC_trils = gsl_matrix_alloc(n_pairs, ((n_vols + window_size) / step));
    // gsl_matrix_set_all(window_FC_trils, -1.5);
    // calculate the FC of each window
    int n_windows = 0;
    window_center = first_center;
    while (window_center <= last_center) {
        window_start = window_center - (window_size/2);
        if (window_start < 0)
            window_start = 0;
        window_end = window_center + (window_size/2);
        if (window_end >= n_vols)
            window_end = n_vols-1;
        gsl_matrix_view bold_window =  gsl_matrix_submatrix(bold, window_start, 0, window_end-window_start+1, n_regions);
        gsl_vector * window_FC_tril = fc_tril(&bold_window.matrix, n_pairs);
        gsl_matrix_set_col(window_FC_trils, n_windows, window_FC_tril);
        window_center += step;
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

double check_fit(int n_regions, int duration, int TR, int FCD_step, int FCD_window_size,
        bool FCD_drop_edges, bool FC_only, char * emp_FC_tril_path, char * emp_FCD_tril_path, char * sim_bold_path) {
    const int n_pairs = ((n_regions) * (n_regions - 1)) / 2;
    const int n_vols = duration / TR;
    const bool write_sim_FC = true;
    // Read simulated bold and create simulated FC (lower triangle)
    gsl_matrix * sim_bold = gsl_matrix_alloc(n_vols, n_regions);
    FILE * sim_bold_file = fopen(sim_bold_path, "r");
    gsl_matrix_fscanf(sim_bold_file, sim_bold);
    fclose(sim_bold_file);
    gsl_vector * sim_FC_tril = fc_tril(sim_bold, n_pairs);
    if (write_sim_FC) {
        FILE * sim_FC_file = fopen("output/sim_FCtril.txt", "w");
        gsl_vector_fprintf(sim_FC_file, sim_FC_tril, "%f");
        fclose(sim_FC_file);
    }
    // Read  empirical FC (lower triangle)
    FILE * emp_FC_tril_file = fopen(emp_FC_tril_path, "r");
    if (emp_FC_tril_file==NULL) {
        printf("Empirical FC file not found. Please include lower triangle of FC (excluding diagonal) in `<input_folder>/emp_FCtril.txt`");
        exit(0);
    }
    gsl_vector * emp_FC_tril = gsl_vector_alloc(n_pairs);
    gsl_vector_fscanf(emp_FC_tril_file, emp_FC_tril);
    fclose(emp_FC_tril_file);
    // Check fit to FC
    double emp_sim_FC_corr = gsl_stats_correlation(
        emp_FC_tril->data, emp_FC_tril->stride,
        sim_FC_tril->data, sim_FC_tril->stride,
        emp_FC_tril->size
    );
    std::cout << "FC correlation: " << emp_sim_FC_corr << "\n";
    if (FC_only) {
        return -emp_sim_FC_corr;
    }
    // Calculate simulated FCD (lower triangle)
    gsl_vector * sim_FCD_tril = fcd_tril(sim_bold, FCD_step, FCD_window_size, FCD_drop_edges);
    FILE * sim_FCD_file = fopen("output/sim_FCDtril.txt", "w");
    gsl_vector_fprintf(sim_FCD_file, sim_FCD_tril, "%f");
    fclose(sim_FCD_file);

    // Read empirical FCD (lower triangle) (size is not known)
    FILE * emp_FCD_tril_file = fopen(emp_FCD_tril_path, "r");
    if (emp_FCD_tril_file==NULL) {
        printf("Empirical FCD file not found. Please include lower triangle of FCD matrix (excluding diagonal) in `<input_folder>/emp_FCDtril.txt`");
        exit(0);
    }
    double tmp;
    double emp_FCD_tril_arr[100000];
    // current compiling settings does not allow creating a dynamic
    // array (I get segmentation error). So I set the maximum size
    // to 100000 elements
    int emp_fcd_tril_length = 0;
    while (fscanf(emp_FCD_tril_file, "%lf", &tmp) != EOF) {
        emp_FCD_tril_arr[emp_fcd_tril_length] = tmp;
        emp_fcd_tril_length++;
    }
    gsl_vector_view emp_FCD_tril_view = gsl_vector_view_array(emp_FCD_tril_arr, emp_fcd_tril_length);
    fclose(emp_FCD_tril_file);

    double ks = ks_stat(
        (&emp_FCD_tril_view.vector)->data, 
        sim_FCD_tril->data, 
        (&emp_FCD_tril_view.vector)->size, 
        sim_FCD_tril->size, 
        NULL);
    std::cout << "FCD KS: " << ks << "\n";
    return ks - emp_sim_FC_corr;
}