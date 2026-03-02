MATLAB scripts used to produce figures relating to the correlation model

Cyclizability_Data folder contains experimental data of sequences and their respective cyclizabilities originally produced in https://doi.org/10.1038/s41586-020-03052-3

distributions.m produces the combined distributions of the Random and Tiling libraries' cyclizabilities

W_ls.m produces figures of mean values and standard deviation for the linear model matrix W over 100 random shuffles, both for the Random and the Tiling libraries, using least-squares and with an option to calculate it through the correlation approach with the Random Library

linear_prediction.m calculates the correlation between the predicted/actual cyclizability with the linear model, either with the least-squares or the correlation function approach

prediction_with_eigenvectors.m calculates the Pearson's correlation mean and standard deviation for the full correlation model and the first 2 modes

pairwise_model folder is code meant for figures in the main manuscript body

eigenvalues.m produces the spectrum of the actual and shuffled data over 10 random train/test splits 

sequence_features.m produces the quadratic sequence features with 90/10 train/test splits

comparison_of_random_tiling.m produces the scatterplot of matrix values achieved with the correlation model on the Random Library vs the least-squares on the Tiling Library, as well as eigenvectors (modes) 





