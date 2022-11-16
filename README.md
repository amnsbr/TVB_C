# TVB_C -- A fast implementation of The Virtual Brain's simulation core in C

The code implements a brain network model composed of connected ReducedWongWang neural mass models (Wong & Wang, 2006) with feedback inhibition control (FIC). For more information on the model and FIC please see Deco et al. (2014), Schirner et al. (2018) and Shen et al. (2019).

## Modifications
On top of the original code in https://github.com/amandakeasson/TVB_C I have made some modifications:
1. FIC is calculated analytically based on the approach suggested by [Demirtaş et al. 2019 Neuron](https://doi.org/10.1016/j.neuron.2019.01.017), and porting the Python code in https://github.com/murraylab/hbnm to C++.
2. Added an option for "disrupting" FIC from the "healthy" regime, i.e., reducing or increasing optimal `J_i` values (identified by the FIC algorithm) by a `J_i_scale` parameter, to be able to model different levels of E-I imbalance, following [Yang et al. 2016 PNAS](https://doi.org/10.1073/pnas.1508436113). The `J_i_scale` parameter by default is 1.0, which reflects no disruption in the FIC.
3. Added an option for introducing regional heterogeneity of J_NMDA in the model based on a heterogeneity map and a `heterogeneity_scale` parameter. This is inspired by (but not exactly implement) the approach used in [Deco et al. 2021 Sci Adv](https://doi.org/10.1126/sciadv.abf4752) and [Demirtaş et al. 2019 Neuron](https://doi.org/10.1016/j.neuron.2019.01.017). (This currently does not work)
4. Added a `check_fit` function (included in `check_fit.cpp`) which calculates the fit of simulated and empirical FC and FCD matrices.

# Usage

```
./tvbii <parameter_file> <path_to_input_folder> <n_all_nodes> <n_real_nodes>
```

Example
```
./tvbii param_set_1 ./test_input 68 68
```

- The first argument specifies a text file that contains parameters (see the file 'param_set_1' for an example and description of paramters)

- The second argument specifies the subject input folder including empirical SC and FC data. It should include the following files:
    - `SC_regionids.txt`: showing which regions are connected to each region in a sparse matrix format
    - `SC_strengths.txt`: showing strength of SC connections (same order as `SC_regionids`)
    - `SC_distances.txt`: showing length of SC connections (same order as `SC_regionids`)
    - `SC_strengths_matrix.txt`: strength of SC connection in the full non-sparse matrix
    - `emp_FCtril.txt`: lower triangle of empirical FC matrix
    - `emp_FCDtril.txt`: lower triangle of empirical FCD matrix (only needed if in `check_fit` `FC_only` is set to false)
    - `heterogeneity.txt`: regional heterogeneity map which is ideally Z-scored (optional)
 
- Results are written into folder _'output'_

- Due to optimization reasons, the number of nodes must be divisible by four. If the number of nodes is not divisible by four, "fake" regions must be added that contain zero coupling to other nodes (except in `SC_strengths_matrix.txt`). The number of total (fake+real) and real nodes must be specified as arguments to the program.


# Dependencies
- [GSL](https://www.gnu.org/software/gsl/) 2.7.1
- [libks](https://github.com/agentlans/libks)

  
# Compilation
  
```
g++  -o tvb -Wall -msse2 -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -m64 -lm \
    tvb.cpp <libks-path>/libks.so <gsl-path>/lib/libgsl.a <gsl-path>/lib/libgslcblas.a \
    -I <gsl-path>/include \
    -I <libks-path>/include
```