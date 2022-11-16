#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
// #include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <random>
#include "fsolve.cpp"
#include "utils.cpp"

/* Model data global variable
  It must be global so that fsolve and the _inh_curr_fixed_pts 
  function can properly use the model parameters. If this function
  was implemented inside the Model class definition, fsolve was
  not able to execute it correclty */
// TODO: convert some of these variables to private Model variables
struct ModelData {
    gsl_matrix *sc, *_Q, *_K_EE, *_K_EI, *_K_IE, *_K_II, *_jacobian,
               *_jacobian_bold, *_Q_bold;
    double w_EE, w_EI, G, w_IE_scale;
    int nc, curr_node_FIC, curr_t_step;
    bool fic_only, _unstable;
    double _W_I, _W_E, _tau_I, _tau_E, _d_I, _d_E, _b_I, _b_E, _a_I, _a_E,
            _gamma_I, _gamma, _sigma, _tau_E_reset, _tau_I_reset, _gamma_reset,
            _gamma_I_reset;
    gsl_vector *_w_II, *_w_IE, *_w_EI, *_w_EE, *_I0, *_J_NMDA, *_I_ext,
                *_I0_E, *_I0_I, *_I_E_ss, *_I_I_ss, *_S_E_ss, *_S_I_ss,
                *_r_E_ss, *_r_I_ss, *_I_E, *_I_I, *_S_E, *_S_I, *_r_E, *_r_I,
                *_dSEdt;
    gsl_vector_complex *_evals, *evals_bold;
    gsl_matrix_complex *_evects, *evects_bold;

    // Excitatory gating variables
    double a_E = 310.0;    // [nC^-1]
    double b_E = 125.0;    // [Hz]
    double d_E = 0.16;    // [s]
    double tau_E = 0.1;   // [s] (NMDA)
    double W_E = 1.0;     // excitatory external input weight

    // Inhibitory gating variables
    double a_I = 615.0;    // [nC^-1]
    double b_I = 177.0;    // [Hz]
    double d_I = 0.087;   // [s]
    double tau_I = 0.01;  // [s] (GABA)
    double W_I = 0.7;     // inhibitory external input weight

    // Other variables from text
    double I_task = 0.02;        // [nA] for (non-RSN case, not currently implemented)
    double gamma = 0.641;        // kinetic conversion factor (typo in text)
    double J_NMDA = 0.15;        // [nA] (excitatory synaptic coupling)
    double I0 = 0.382;           // [nA] (overall effective external input)
    double sigma = 0.01;         // [nA] (noise amplitude)

    // Additional long-range projecting current term (used in single-node case)
    double I_ext = 0.0;  // [nA]

    // Synaptic weight connections
    double w_II = 1.0;     // I->I self-coupling
    double w_IE = 1.0;     // E->I coupling

    // Steady-state solutions in isolated case
    double r_I_ss = 3.9218448633;  // Hz
    double r_E_ss = 3.0773270642;  // Hz
    double I_I_ss = 0.2528951325;  // nA
    double I_E_ss = 0.3773805650;  // nA
    double S_I_ss = 0.0392184486;  // dimensionless
    double S_E_ss = 0.1647572075;  // dimensionless
};
ModelData md;

// transfer function and derivatives for Excitatory and Inhibitory populations

double _phi_E(double IE) {
    return ((md._a_E * IE) - md._b_E) / (1 - exp(-md._d_E * ((md._a_E * IE) - md._b_E)));
}

double _dphi_E(double IE) {
    return (
        (md._a_E * (1 - exp(-1 * md._d_E * ((md._a_E * IE) - md._b_E))))
        - (md._a_E * md._d_E * exp(-1 * md._d_E * ((md._a_E * IE) - md._b_E)) * ((md._a_E * IE) - md._b_E))
        ) / pow((1 - exp(-1 * md._d_E * ((md._a_E * IE) - md._b_E))), 2);
}

double _phi_I(double II) {
    return ((md._a_I * II) - md._b_I) / (1 - exp(-1 * md._d_I * ((md._a_I * II) - md._b_I)));
}

double _dphi_I(double II) {
    return (
        md._a_I * (1 - exp(-1 * md._d_I * ((md._a_I * II) - md._b_I)))
        - md._a_I * md._d_I * exp(-1 * md._d_I * ((md._a_I * II) - md._b_I)) * ((md._a_I * II) - md._b_I)
        ) / pow((1 - exp(-1 * md._d_I * ((md._a_I * II) - md._b_I))), 2);
}

/* Eq.10 in Demirtas which would be used in `fsolve`
 to find the steady-state inhibitory synaptic gating variable
 and the suitable w_IE weight according to the FIC algorithm */

void _inh_curr_fixed_pts(int n, double x[], double fx[]) {
    fx[0] = md._I0_I->data[md.curr_node_FIC] + md._w_EI->data[md.curr_node_FIC] * md._S_E_ss->data[md.curr_node_FIC] -
            md._w_II->data[md.curr_node_FIC] * md._gamma_I * md._tau_I * _phi_I(x[0]) - x[0];
    return;
}

/* Just to make sure fsolve works correctly */

void test_fsolve(int n, double x[], double fx[]) {
    fx[0] = (x[0] * x[0]) - 9;
}

class Model {
    public:
        Model(gsl_matrix * sc, double G, double wee, double wei, double wie_scale) {
            md.sc = sc;
            md.nc = sc->size1;
            md.G = G;
            md.w_EE = wee;
            md.w_EI = wei;
            md.w_IE_scale = wie_scale;

            // // # Unstable if Jacobian has eval > 0
            md._unstable = false;

            // // # Initialize model outputs to None
            md._jacobian = gsl_matrix_alloc(md.nc*2, md.nc*2);
            // // self._cov = None
            // // self._corr = None
            // // self._cov_bold = None
            // // self._corr_bold = None
            // // self._full_cov = None

            // // # Initialize state members to None
            md._I_E = gsl_vector_alloc(md.nc);
            md._I_I = gsl_vector_alloc(md.nc);
            md._S_E = gsl_vector_alloc(md.nc);
            md._S_I = gsl_vector_alloc(md.nc);
            md._r_E = gsl_vector_alloc(md.nc);
            md._r_I = gsl_vector_alloc(md.nc);

            // Various model parameters
            md._w_II = repeat(md.w_II, md.nc);
            md._w_IE = repeat(md.w_IE, md.nc);
            md._w_EI = repeat(md.w_EI, md.nc);
            md._w_EE = repeat(md.w_EE, md.nc);

            md._I0 = repeat(md.I0, md.nc);
            md._J_NMDA = repeat(md.J_NMDA, md.nc);
            md._sigma = md.sigma;
            md._gamma = md.gamma;
            md._W_I = md.W_I;
            md._W_E = md.W_E;
            md._tau_I = md.tau_I;
            md._tau_E = md.tau_E;
            md._d_I = md.d_I;
            md._d_E = md.d_E;
            md._b_I = md.b_I;
            md._b_E = md.b_E;
            md._a_I = md.a_I;
            md._a_E = md.a_E;
            md._I_ext = repeat(md.I_ext, md.nc);

            md._gamma_I = 1.0;

            md._tau_E_reset = md._tau_E;
            md._tau_I_reset = md._tau_I;
            md._gamma_reset = md._gamma;
            md._gamma_I_reset = md._gamma_I;

            // Baseline input currents
            md._I0_E = vector_scale(md._I0, md._W_E);
            md._I0_I = vector_scale(md._I0, md._W_I);

            // Steady state values for isolated node
            md._I_E_ss = repeat(md.I_E_ss, md.nc);
            md._I_I_ss = repeat(md.I_I_ss, md.nc);
            md._S_E_ss = repeat(md.S_E_ss, md.nc);
            md._S_I_ss = repeat(md.S_I_ss, md.nc);
            md._r_E_ss = repeat(md.r_E_ss, md.nc);
            md._r_I_ss = repeat(md.r_I_ss, md.nc);

            // Noise covariance matrix
            md._Q = mul_eye(md._sigma * md._sigma, 2 * md.nc);
        }

        gsl_vector * _analytic_FIC() {
            /*
                Implements feedback inhibition control analytically.
                Sets w_IE in each node so that its excitatory activity
                is ~ 3 Hz
                Credit: This is hbnm.model.dmf.Model._analytic_FIC function
                ported to C++
            */
            gsl_vector * out = repeat(md.w_IE, md.nc);
            double x[1], fx[1]; 
            for (md.curr_node_FIC=0; md.curr_node_FIC<md.nc; md.curr_node_FIC++) {
                // fsolve usage based on 
                // https://people.sc.fsu.edu/~jburkardt/cpp_src/fsolve_test/fsolve_test.cpp
                x[0] = md.I_I_ss;
                double tol = 0.000000014912; // same as scipy default
                double *wa;
                int n = 1; // number of variables/equations
                int lwa = ( n * ( 3 * n + 13 ) ) / 2;
                wa = new double[lwa];
                // fsolve(test_fsolve, 1, x, fx, tol, wa, lwa);
                fsolve(_inh_curr_fixed_pts, 1, x, fx, tol, wa, lwa);
                gsl_vector_set(md._I_I, md.curr_node_FIC, x[0]);
                gsl_vector_set(md._I_I_ss, md.curr_node_FIC, x[0]);
                gsl_vector_set(md._r_I, md.curr_node_FIC, _phi_I(x[0]));
                gsl_vector_set(md._r_I_ss, md.curr_node_FIC, _phi_I(x[0]));
                gsl_vector_set(md._S_I_ss, md.curr_node_FIC, 
                               _phi_I(x[0]) * md._tau_I * md._gamma_I);
                gsl_vector * _K_EE_row = gsl_vector_alloc(md.nc);
                gsl_matrix_get_row(_K_EE_row, md._K_EE, md.curr_node_FIC);
                double _K_EE_dot_S_E_ss;
                gsl_blas_ddot(_K_EE_row, md._S_E_ss, &_K_EE_dot_S_E_ss);
                double J = (-1 / md._S_I_ss->data[md.curr_node_FIC]) *
                           (md._I_E_ss->data[md.curr_node_FIC] - 
                            md._I_ext->data[md.curr_node_FIC] - 
                            md._I0_E->data[md.curr_node_FIC] -
                            _K_EE_dot_S_E_ss);
                if (J < 0) {
                    md._unstable = true;
                    std::cout << "FIC calculation led to negative J values!\n";
                    exit(0);
                }
                gsl_vector_set(out, md.curr_node_FIC, J);
            }
            return out;
        }

        void set_jacobian() {
            // Excitatory connection weights
            md._K_EE = gsl_matrix_alloc(md.nc, md.nc);
            gsl_matrix_memcpy(md._K_EE, md.sc);
            gsl_matrix_scale(md._K_EE, md.G * md.J_NMDA);
            gsl_matrix * _w_EE_matrix = mul_eye(md.w_EE, md.nc); // for heterogeneous model this won't work
            gsl_matrix_add(md._K_EE, _w_EE_matrix);
            md._K_EI = mul_eye(md.w_EI, md.nc); // for heterogeneous model this won't work

            // FIC
            md._w_IE = _analytic_FIC();
            // multiply it by the parameter w_IE_scale (disrupt FIC)
            if (md.w_IE_scale != 1.0) {
                gsl_vector_scale(md._w_IE, md.w_IE_scale);
            }

            if (md.fic_only) {
                return;
            }
            // Inhibitory connection weights
            md._K_IE = make_diag(vector_scale(md._w_IE, -1));
            md._K_II = mul_eye(-1 * md.w_II, md.nc);

            // Derivatives of transfer function for each cell type
            // at steady state value of current
            gsl_matrix * dr_E = gsl_matrix_alloc(md.nc, md.nc);
            gsl_matrix * dr_I = gsl_matrix_alloc(md.nc, md.nc);
            for (int i=0; i<md.nc; i++) {
                gsl_matrix_set(dr_E, i, i, 
                            _dphi_E(md._I_E_ss->data[i]));
                gsl_matrix_set(dr_I, i, i, 
                            _dphi_I(md._I_I_ss->data[i]));
            }
            // A_EE = (-1. / self._tau_E - (self._gamma * self._r_E_ss)) * eye + \
            //    ((-self._gamma * (self._S_E_ss - 1.)) * eye).dot(dr_E.dot(self._K_EE))
            gsl_vector * _A_EE_a = gsl_vector_alloc(md.nc);
            gsl_vector_set_all(_A_EE_a, -1 / md._tau_E);
            gsl_vector_add(_A_EE_a, vector_scale(md._r_E_ss, -1 * md._gamma));
            gsl_matrix * _A_EE_a_diag = make_diag(_A_EE_a);
            gsl_vector * _A_EE_b = gsl_vector_alloc(md.nc);
            gsl_vector_memcpy(_A_EE_b, md._S_E_ss);
            gsl_vector_add_constant(_A_EE_b, -1.0);
            gsl_vector_scale(_A_EE_b, -1 * md._gamma);
            gsl_matrix * _A_EE_b_diag = make_diag(_A_EE_b);
            gsl_matrix * A_EE = dot(_A_EE_b_diag, dot(dr_E, md._K_EE));
            gsl_matrix_add(A_EE, _A_EE_a_diag);
            // A_IE = ((self._gamma * (1. - self._S_E_ss)) * eye).dot(dr_E.dot(self._K_IE))
            gsl_vector * _A_IE_a = vector_scale(md._S_E_ss, -1);
            gsl_vector_add_constant(_A_IE_a, 1);
            gsl_vector_scale(_A_IE_a, md._gamma);
            gsl_matrix * A_IE = dot(make_diag(_A_IE_a), dot(dr_E, md._K_IE));
            // A_EI = self._gamma_I * dr_I.dot(self._K_EI)
            gsl_matrix * A_EI = dot(dr_I, md._K_EI);
            gsl_matrix_scale(A_EI, md._gamma_I);
            // A_II = (-1. / self._tau_I) * eye + self._gamma_I * dr_I.dot(self._K_II)
            gsl_matrix * A_II = mul_eye((-1 / md._tau_I), md.nc);
            gsl_matrix * A_II_b = dot(dr_I, md._K_II);
            gsl_matrix_scale(A_II_b, md._gamma_I);
            gsl_matrix_add(A_II, A_II_b);
            // Stack blocks to form full Jacobian
            // col1 = np.vstack((A_EE, A_EI)); col2 = np.vstack((A_IE, A_II))
            // self._jacobian = np.hstack((col1, col2))
            md._jacobian = concat({
                concat({A_EE, A_IE}, 1),
                concat({A_EI, A_II}, 1)
            }, 0);

            // Eigenvalues of Jacobian matrix
            // self._evals, self._evects = eig(self._jacobian)
            // Following gsl docs example for nonsymmetric matrices :
            // https://www.gnu.org/software/gsl/doc/html/eigen.html
            // first create a copy because eigenvalue calculations alters the matrix (don't know why!)
            int jacobian_size = md.nc*2;
            gsl_matrix * jacobian_copy = gsl_matrix_alloc(jacobian_size, jacobian_size);
            gsl_matrix_memcpy(jacobian_copy, md._jacobian);

            md._evals = gsl_vector_complex_alloc(jacobian_size);
            md._evects = gsl_matrix_complex_alloc(jacobian_size, jacobian_size);
            gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc(jacobian_size);
            gsl_eigen_nonsymmv(jacobian_copy, md._evals, md._evects, w);
            gsl_eigen_nonsymmv_free(w);
            // Check stability using eigenvalues
            // self._max_eval = np.real(self._evals.max())
            // self._unstable = self._max_eval >= 0.0
            gsl_vector_view eval_real_view = gsl_vector_complex_real(md._evals);
            double max_eval = gsl_vector_max(&eval_real_view.vector);
            if (max_eval >= 0) {
                md._unstable = true;
                std::cout << "Model is unstable based on Jacobian eigenvalues";
                exit(0);
            }
        }
};

void do_fic(float * J_i, int n_regions, char * sc_path, float G, float w_EE_bias, float w_EI_bias, float J_i_scale) {
    md.fic_only = true;
    gsl_matrix * sc = gsl_matrix_alloc(n_regions, n_regions);
    FILE * sc_file = fopen(sc_path, "r");
    gsl_matrix_fscanf(sc_file, sc);
    // initialize model with parameters G, wee, wei, wie_scale
    Model m(sc, G, w_EE_bias, w_EI_bias, J_i_scale);
    m.set_jacobian();
    for (int i=0; i<n_regions; i++) {
        J_i[i] = (float)gsl_vector_get(md._w_IE, i);
    }
    // TODO: consider the case that FIC fails (negative Ji)
    return;
}