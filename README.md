BEST_urchin_model.R is an R file, open this first. Should be able to open this, make sure Rstan is installed. Otherwise just step through this file,
first part is setting up, second part is running stan model, third part saving that (or reading from already run), then the final (majority of the remaining code) is to 
generate all the figures in the paper.

best_stan_urchin_model.stan is the stan file for the main urchin model.

historical_catch.csv is historical catch data

survey_data.csv is the data to which the model is fit - fisheries-independent surveyed urchin densities in 9 regions down the east coast of Tasmania for 3 separate survey years

size_transition_matrix is the size transition matrix for growth of Longspined sea urchins through time
You can also create this matrix by running
sizetransition_code.R
