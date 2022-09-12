# TVB_C -- A fast implementation of The Virtual Brain's simulation core in C

The code implements a brain network model composed of connected ReducedWongWang neural mass models (Wong & Wang, 2006) with feedback inhibition control (FIC). For more information on the model and FIC please see Deco et al. (2014), Schirner et al. (2018) and Shen et al. (2019).

For more information on The Virtual Brain (TVB) please see 
www.thevirtualbrain.org

For questions and other enquiries please contact 

Michael Schirner (m.schirner@fu-berlin.de) or 
Petra Ritter (petra.ritter@charite.de).

## Modifications
On top of the original code in https://github.com/amandakeasson/TVB_C I have made some modifications:
1. Added an option for "disrupting" FIC from the "healthy" regime, i.e., reducing or increasing optimal `J_i` values (identified by the FIC algorithm) by a `J_i_scale` parameter, to be able to model different levels of E-I imbalance, following [Yang et al. 2016 PNAS](https://doi.org/10.1073/pnas.1508436113). The `J_i_scale` parameter by default is 1.0, which reflects no disruption in the FIC.
2. Added an option for introducing regional heterogeneity of J_NMDA in the model based on a heterogeneity map and a `heterogeneity_scale` parameter. This is inspired by (but not exactly implement) the approach used in [Deco et al. 2021 Sci Adv](https://doi.org/10.1126/sciadv.abf4752) and [Demirtaş et al. 2019 Neuron](https://doi.org/10.1016/j.neuron.2019.01.017).

# Usage

```
./tvbii <parameter_file> <subject_id>
```

Example
```
./tvbii param_set_1 UE_20120803
```

• The first argument specifies a text file that contains parameters (see the file 'param_set_1' for an example and description of paramters)

• The second argument specifies the subject-id for the input files contained in the folder _'input'_. Each of the three input files must have <subject_id> as prefix and as suffix either _"\_SC\_strengths.txt"_, _"\_SC\_distances.txt"_ or, _"\_SC\_regionids.txt"_. The three files specify structural connectivity of the brain network model in a sparse matrix format. The enclosed Matlab script _generate_input_SC.m_ generates these input files from Matlab's standard matrix format.

• Results are written into folder _'output'_. File-schema: BOLD_<parameter_file>.txt; the first n lines (n=number of regions) contain two columns each that contain the sum of all input strengths for that region and the J_i value found during FIC tuning, respectively. The following t lines contain n columns each and contain simulated fMRI BOLD activity for t time points.

• The relative folder structure, i.e., the location of the folders 'input' and 'output' relative to the program binary needs to remain stable, otherwise the program won't be able to read or write data.

• Due to optimization reasons, the number of nodes must be divisible by four. If the number of nodes is not divisible by four, "fake" regions must be added that contain zero coupling to other nodes (all zeros in strength matrix).

• To introduce heterogeneity of J_NMDA in the model, a txt file with <subject_id> as prefix and `_heterogeneity.txt` should be added to the input folder, with one value in each line for each region. A value of 0 in a region means no divergence from the J_NMDA_bias, but positive/negative values make the J_NMDA of a region higher/lower than the bias term: `J_NMDA[i] = J_NMDA_bias * (1 + (heterogeneity[i] * heterogeneity_bias))`, with `heterogeneity` being the heterogeneity map and `heterogeneity_bias` the free parameters that determines the degree to which there's regional heterogeneity of J_NMDA. If this parameter is set to 0 or the txt file does not exist a homogeneous model will be simulated.
  
# Compilation
  
  Compile the C code with a compiler with support for SSE instructions enabled. I found the following combination of flags yields good performance (in terms of simulation speed) with the GNU C compiler

```
gcc -Wall -std=c99 -msse2 -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -m64 -lm main.c -o tvb
```
  
# References
  
  Deco, G., Ponce-Alvarez, A., Hagmann, P., Romani, G. L., Mantini, D., & Corbetta, M. (2014). How local excitation–inhibition ratio impacts the whole brain dynamics. Journal of Neuroscience, 34(23), 7886-7898.
  
  Schirner, M., McIntosh, A. R., Jirsa, V., Deco, G., & Ritter, P. (2018). Inferring multi-scale neural mechanisms with brain network modelling. Elife, 7, e28927.
  
  Shen, K., Bezgin, G., Schirner, M., Ritter, P., Everling, S., McIntosh, A. R. (2019) A macaque connectome for large-scale network simulations in TheVirtualBrain. (under review)
  
  Wong, K. F., & Wang, X. J. (2006). A recurrent network mechanism of time integration in perceptual decisions. Journal of Neuroscience, 26(4), 1314-1328.
