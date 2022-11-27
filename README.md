# Calibration-of-Functional-Input


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
These codes provide Matlab implementation of the examples in Section 6, Appendix F and Appendix G of the paper.

Implementation of the example in Section 6 uses codes contained in the folders 'Example_parabolic' and 'Example_parabolic_sparse_obs'.
'Example_parabolic': All codes for Section 6 except those for comparison of design criteria based on the sparse observation grid (P^b*T).
'Example_parabolic_sparse_obs': Codes for comparison of design criteria based on the sparse observation grid (P^b*T) as in section 6.2. This gives results summarized in Figure 2(b) in the paper.

Implementation of the example in Appendix F uses codes contained in the folder 'Example_supplement_elliptic_F'.

Implementation of the example in Appendix G uses codes contained in the folder 'Example_supplement_toy_G'.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
To obtain results in Section 6.2, open the folder titled "Example_parabolic" and then run the .m file 'compare_design_size' and 'compare_design_criteria'. 
	                        Next, open the folder titled "Example_parabolic_sparse" and then run the .m file 'compare_design_criteria.m'.
To obtain results in Section 6.3, open the folder titled "Example_parabolic" and then run the .m file 'compare_prior_specification.m'.
To obtain results in Section 6.4, open the folder titled "Example_parabolic" and then run the .m file 'compare_MCMC.m'.

To obtain results in Appendix F.2, open the folder titled "Example_supplement_elliptic_F" and then run the .m file 'compare_design_size.m' and 'compare_design_criteria.m'.
To obtain results in Appendix F.3, open the folder titled "Example_supplement_elliptic_F" and then run the .m file 'compare_prior_specification.m'.
To obtain results in Appendix F.4, open the folder titled "Example_supplement_elliptic_F" and then run the .m file 'compare_different_M.m'.
To obtain results in Appendix F.5, open the folder titled "Example_supplement_elliptic_F" and then run the .m file 'compare_MCMC.m'.

To obtain results in Appendix G, open the folder titled 'Example_supplement_toy_G' and then run the .m file 'compare_with_exact_method.m' and 'compare_with_exact_method_MCMC.m';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
After Matlab has finished with the computations, the results are stored in the root folder in .mat files： 'RES_*.mat' and 'ESS_*.mat', where '*' can be 'prior', 'criteria', 'size', 'M', 'MCMC', 'exact'  which correspond to the different simulations mentioned above.
'ESS_*.mat'： gives the relative L2 errors.
'RES_*.mat'： includes other data arising from the simulation. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
When MATLAB finishes, the variables below are important. 

Important variables after running 'compare_prior_specification.m': 
ESS_post： stores the relative L2 errors of the posterior mode of the functional input. 
ESS_post(i,j): the relative L2 error of the posterior mode of the functional input for the i-th simulation run, and j-th prior model. j=1: Unif-GP, j=2: L-GP; j=3: M-GP; j=4: U-GP.

Important variables after running 'compare_design_criteria.m': 
ESS_post： stores the relative L2 errors of the posterior mode of the functional input. 
ESS_post(i,j): the relative L2 error of the posterior mode of the functional input for the i-th simulation run, and j-th design criterion. j=1: WPV; j=2: VL; j=3: AEI.

Important variables after running 'compare_design_size.m': 
ESS_post： stores the relative L2 errors of the posterior mode of the functional input. 
ESS_post(i,j): the relative L2 error of the posterior mode of the functional input for the i-th simulation run, and j-th follow-up design size. j=1: q_1=0 ; j=2: q_1=20 ; j=3: q_1=50.

Important variables after running 'compare_different_M.m': 
ESS_post： stores the relative L2 errors of the posterior mode of the functional input. 
ESS_post(i,j): the relative L2 error of the posterior mode of the functional input for the i-th simulation run, and j-th value of M. j=1: M=3; j=2: M=9 ; j=3: M=12 ; j=4: M=15.

Important variables after running 'compare_with_exact_method.m':
ESS_post： stores the relative L2 errors of the posterior mode of the functional input. 
ESS_post(i,j): the relative L2 error of the posterior mode of the functional input for the i-th simulation run, and j-th method for functional parameter calibration. j=1: proposed method; j=2: no-emulator method; j=3: exact method.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
If you want to plot the figures, do the following AFTER running all .m files mentioned above: 

1. Copy the .mat file 'ESS_criteria_sparse.mat' from folder 'Example_parabolic_sparse_obs' to folder 'Example_parabolic'. Then, 
run 'plot_parabolic.m' in folder 'Example_parabolic'. All figures in Section 6 will be plotted and stored in the root folder.

2. Run 'plot_elliptic.m' in folder 'Example_supplement_elliptic_F'. All figures in Appendix F will be plotted and stored in the root folder.

3. Run 'plot_toy.m' in folder 'Example_supplement_toy_G'. All figures in Appendix G will be plotted and stored in the root folder.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Note: For the codes in the folders 'Example_parabolic' and 'Example_parabolic_sparse_obs', i.e., the codes for reproducing results in Section 6,
the unit for the length dimension is m (instead of km). For example,
the unit for x is m;
the unit for T is m^2/s;
the unit for the PDE solution u is m;
the unit for the Dirichlet boundary condition is m;
the unit for h(x) is m/s.
...

For the codes in the folder "Example_supplement_elliptic_F", i.e., the codes for reproducing results in Appendix F,
the unit for the length dimension is m (instead of km). For example, 
the unit for x is m;
the unit for T is m^2/day;
the unit for the PDE solution u is m;
the unit for the Dirichlet boundary condition is m;
the unit for h(x) is m/day.
...
