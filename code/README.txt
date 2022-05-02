***********************************************************************
Readme for Nested Sampling of Radial Velocity Posterior code
***********************************************************************

Applied Math 207, Spring 2013

=== Team Members ===

   Marion Dierickx
   Xinyi Guo
   Philip Mocz

=== Description ===

Final project. Fit multiple planet RV data in Bayesian approach using nested sampling

=== File Listing ===

 code/                       Matlab .m files
   compare_data_ext.m                  plot extended RV data
   diver_extra_figures.m               generate additional figures from the results
   driver_mock_data.m                  * main script file for fitting synthetic data
   driver_planet_mass_analysis.m       calculate planet mass and semi-major axis after fitting orbital parameters
   driver_real_data.m                  * main script file for fitting real data
   generate_synthetic_rv_data.m        generates synthetic radial velocity observations
   load_result_plot.m                  additional plots for real data
   nested_sample.m                     performs nested sampling of posterior
   plot_posterior.m                    plots posterior
   rv_model.m                          calculates radial velocity of planet given orbital parameters
   sample_prior.m                      samples prior distributions


 code/output                 all outputs from runs kept here           

      