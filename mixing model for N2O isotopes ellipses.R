# Mixing Model for N2O Ellipse dual-isotope paper: Snider, Venkiteswaran, Schiff, and Spoelstra
# https://github.com/jjvenky/Global-N2O-Ellipses
#
# Mixing model MixSIAR
# Code modified from https://github.com/brianstock/MixSIAR/blob/master/R/mixsiar_script.R (Brian Stock, July 2014)



###
# Load packages
library(R2jags)
library(plyr)
library(ggplot2)
require(MASS)
require(RColorBrewer)
require(reshape)
library(devtools) # used to load the code direct from Github
###



# Generates one random number (else JAGS can complain)
runif(1)



# Source the model files for MixSIAR instead of downloading and unzipping or pulling from Github
source_url("https://raw.github.com/brianstock/MixSIAR/v2.1.2/R/load_mix_data.r")
source_url("https://raw.github.com/brianstock/MixSIAR/v2.1.2/R/load_source_data.r")
source_url("https://raw.github.com/brianstock/MixSIAR/v2.1.2/R/load_discr_data.r")
# source_url("https://raw.github.com/brianstock/MixSIAR/v2.1.2/R/plot_data.r")
source_url("https://raw.github.com/jjvenky/MixSIAR/master/R/plot_data.r")
source_url("https://raw.github.com/brianstock/MixSIAR/v2.1.2/R/write_JAGS_model.r")
source_url("https://raw.github.com/brianstock/MixSIAR/v2.1.2/R/run_model.r")
source_url("https://raw.github.com/brianstock/MixSIAR/v2.1.2/R/output_JAGS.r")
source_url("https://raw.github.com/brianstock/MixSIAR/v2.1.2/R/plot_continuous_var.r")



# Have to create four sets of models scenarios:
# Model scenario 1: 2 isotopes (δ15N, δ18O), 3 end-members (Freshwater, Marine, Soil)
# Model scenario 2: 3 isotopes (δ15N, δ18O, SP), 3 end-members (Freshwater, Marine, Soil)
# Model scenario 3: 2 isotopes (δ15N, δ18O), 2 end-members (Continental, Marine)
# Model scenario 4: 3 isotopes (δ15N, δ18O, SP), 2 end-members (Continental, Marine)
# Continental == Freshwater + Soils + Urban Wastewater




###
# BEGIN: Model scenario 1: 2 isotopes (δ15N, δ18O), 3 end-members (Freshwater, Marine, Soil)

# Save the necessary data as csv files since MixSIAR likes to load from csv files
# Consumer is the tropospheric value we want. It is the Modern antthro+natural value from Röckmann et al. 2003
# Sources are the subsetting Freshwater, Marine, and Soils data
# Discrimination is isotope fractionation between Sources and Consumer -- set at 0
write.csv(subset(globModSol["10",], select = c("d15N", "d18Ovsmow")), 
          file = "N2O_consumer_2-isotopes.csv", row.names = FALSE)
write.csv(rename(na.omit(subset(subFMS, select = c("Category", "d15N", "d18Ovsmow"))), 
                 c("Category" = "Source")), 
          file = "N2O_sources_2-isotopes.csv", row.names = FALSE)
write.csv(data.frame(Source = c("Freshwater", "Marine", "Soil"),
                     Meand15N = 0, SDd15N = 0, Meand18Ovsmow = 0, SDd18Ovsmow = 0),
          file = "N2O_discrimination_2-isotopes.csv", row.names = FALSE)


# Load the data from csv files using MixSIAR scripts
mix <- load_mix_data(filename="N2O_consumer_2-isotopes.csv", iso_names=c("d15N", "d18Ovsmow"), 
                     random_effects=NULL, cont_effects=NULL, fixed_effects=NULL)
source <- load_source_data(filename="N2O_sources_2-isotopes.csv", 
                           source_factors=NULL, conc_dep=FALSE, data_type="raw", mix)    
discr <- load_discr_data(filename="N2O_discrimination_2-isotopes.csv", mix)


# Plot data to check it was loaded correctly
plot_data(filename="N2O_isospace_plot_2-isotopes", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)



# Write JAGS model file
model_filename <- "MixSIAR_model_N2O_2-isotopes.txt"   # Name of the JAGS model file
indiv_effect <- FALSE                    # Include Individual as a random effect in the model?
nested <- FALSE                         # If there are 2 random effects, is the 2nd nested in the 1st (hierarchical)?
write_JAGS_model(model_filename, indiv_effect, nested, resid_err=TRUE, mix,source)



# Run model
# JAGS output will be saved as 'jags.1isotopes'
#
# MCMC run options:
# run <- "test" # list(chainLength=1000, burn=500, thin=1, chains=3, calcDIC=TRUE)
# run <- "very short" # list(chainLength=10000, burn=5000, thin=5, chains=3, calcDIC=TRUE)
# run <- "short" # list(chainLength=50000, burn=25000, thin=25, chains=3, calcDIC=TRUE)
# run <- "normal" # list(chainLength=100000, burn=50000, thin=50, chains=3, calcDIC=TRUE)
# run <- "long" # list(chainLength=300000, burn=200000, thin=100, chains=3, calcDIC=TRUE)
# run <- "very long" # list(chainLength=1000000, burn=700000, thin=300, chains=3, calcDIC=TRUE)
# run <- "extreme" # list(chainLength=3000000, burn=2700000, thin=300, chains=3, calcDIC=TRUE)

jags.1isotopes <- run_model(run="normal", indiv_effect,mix,source,discr,model_filename)



# Process JAGS output
output_options.1isotopes <- list(summary_save = TRUE,                   # Save the summary statistics as a txt file?
                                 summary_name = "N2O_summary_statistics_2-isotopes",    # If yes, specify the base file name (.txt will be appended later)
                                 sup_post = FALSE,                       # Suppress posterior density plot output in R?
                                 plot_post_save_pdf = TRUE,            # Save posterior density plots as pdfs?
                                 plot_post_name = "N2O_posterior_density_2-isotopes",   # If yes, specify the base file name(s) (.pdf/.png will be appended later)
                                 sup_pairs = FALSE,                      # Suppress pairs plot output in R?
                                 plot_pairs_save_pdf = TRUE,            # Save pairs plot as pdf?
                                 plot_pairs_name = "N2O_pairs_plot_2-isotopes",         # If yes, specify the base file name (.pdf/.png will be appended later)
                                 sup_xy = FALSE,                         # Suppress xy/trace plot output in R?
                                 plot_xy_save_pdf = TRUE,               # Save xy/trace plot as pdf?
                                 plot_xy_name = "N2O_xy_plot_2-isotopes",               # If yes, specify the base file name (.pdf/.png will be appended later)
                                 gelman = TRUE,                          # Calculate Gelman-Rubin diagnostic test?
                                 heidel = FALSE,                         # Calculate Heidelberg-Welch diagnostic test?
                                 geweke = TRUE,                          # Calculate Geweke diagnostic test?
                                 diag_save = TRUE,                       # Save the diagnostics as a txt file?
                                 diag_name = "N2O_diagnostics_2-isotopes",              # If yes, specify the base file name (.txt will be appended later)
                                 indiv_effect = indiv_effect,            # Is Individual a random effect in the model? (already specified)
                                 plot_post_save_png = FALSE,             # Save posterior density plots as pngs?
                                 plot_pairs_save_png = FALSE,            # Save pairs plot as png?
                                 plot_xy_save_png = FALSE)               # Save xy/trace plot as png?

output_JAGS(jags.1isotopes, mix, source, output_options.1isotopes)

# END: Model scenario 1: 2 isotopes (δ15N, δ18O), 3 end-members (Freshwater, Marine, Soil)
###



###
# BEGIN: Model scenario 2: 3 isotopes (δ15N, δ18O, SP), 3 end-members (Freshwater, Marine, Soil)

# Save the necessary data as csv files since MixSIAR likes to load from csv files
# Consumer is the tropospheric value we want. It is the Modern antthro+natural value from Röckmann et al. 2003
# Sources are the subsetting Freshwater, Marine, and Soils data
# Discrimination is isotope fractionation between Sources and Consumer -- set at 0
write.csv(cbind(subset(globModSol["10",], select = c("d15N", "d18Ovsmow")), "SP" = 12.7), 
          file = "N2O_consumer_3-isotopes.csv", row.names = FALSE)
write.csv(rename(na.omit(subset(subFMS, select = c("Category", "d15N", "d18Ovsmow", "SP"))), 
                 c("Category" = "Source")), 
          file = "N2O_sources_3-isotopes.csv", row.names = FALSE)
write.csv(data.frame(Source = c("Freshwater", "Marine", "Soil"),
                     Meand15N = 0, SDd15N = 0, Meand18Ovsmow = 0, SDd18Ovsmow = 0, MeanSP = 0, SDSP = 0),
          file = "N2O_discrimination_3-isotopes.csv", row.names = FALSE)


# Load the data from csv files using MixSIAR scripts
mix <- load_mix_data(filename="N2O_consumer_3-isotopes.csv", iso_names=c("d15N", "d18Ovsmow", "SP"), 
                     random_effects=NULL, cont_effects=NULL, fixed_effects=NULL)
source <- load_source_data(filename="N2O_sources_3-isotopes.csv", 
                           source_factors=NULL, conc_dep=FALSE, data_type="raw", mix)    
discr <- load_discr_data(filename="N2O_discrimination_3-isotopes.csv", mix)


# Plot data to check it was loaded correctly
plot_data(filename="N2O_isospace_plot_3-isotopes", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)



# Write JAGS model file
model_filename <- "MixSIAR_model_N2O_3-isotopes.txt"   # Name of the JAGS model file
indiv_effect <- FALSE                    # Include Individual as a random effect in the model?
nested <- FALSE                         # If there are 2 random effects, is the 2nd nested in the 1st (hierarchical)?
write_JAGS_model(model_filename, indiv_effect, nested, resid_err=TRUE, mix,source)



# Run model
# JAGS output will be saved as 'jags.2isotopes'
#
# MCMC run options:
# run <- "test" # list(chainLength=1000, burn=500, thin=1, chains=3, calcDIC=TRUE)
# run <- "very short" # list(chainLength=10000, burn=5000, thin=5, chains=3, calcDIC=TRUE)
# run <- "short" # list(chainLength=50000, burn=25000, thin=25, chains=3, calcDIC=TRUE)
# run <- "normal" # list(chainLength=100000, burn=50000, thin=50, chains=3, calcDIC=TRUE)
# run <- "long" # list(chainLength=300000, burn=200000, thin=100, chains=3, calcDIC=TRUE)
# run <- "very long" # list(chainLength=1000000, burn=700000, thin=300, chains=3, calcDIC=TRUE)
# run <- "extreme" # list(chainLength=3000000, burn=2700000, thin=300, chains=3, calcDIC=TRUE)

jags.2isotopes <- run_model(run="long", indiv_effect,mix,source,discr,model_filename)



# Process JAGS output
output_options.2isotopes <- list(summary_save = TRUE,                   # Save the summary statistics as a txt file?
                                 summary_name = "N2O_summary_statistics_3-isotopes",    # If yes, specify the base file name (.txt will be appended later)
                                 sup_post = FALSE,                       # Suppress posterior density plot output in R?
                                 plot_post_save_pdf = TRUE,            # Save posterior density plots as pdfs?
                                 plot_post_name = "N2O_posterior_density_3-isotopes",   # If yes, specify the base file name(s) (.pdf/.png will be appended later)
                                 sup_pairs = FALSE,                      # Suppress pairs plot output in R?
                                 plot_pairs_save_pdf = TRUE,            # Save pairs plot as pdf?
                                 plot_pairs_name = "N2O_pairs_plot_3-isotopes",         # If yes, specify the base file name (.pdf/.png will be appended later)
                                 sup_xy = FALSE,                         # Suppress xy/trace plot output in R?
                                 plot_xy_save_pdf = TRUE,               # Save xy/trace plot as pdf?
                                 plot_xy_name = "N2O_xy_plot_3-isotopes",               # If yes, specify the base file name (.pdf/.png will be appended later)
                                 gelman = TRUE,                          # Calculate Gelman-Rubin diagnostic test?
                                 heidel = FALSE,                         # Calculate Heidelberg-Welch diagnostic test?
                                 geweke = TRUE,                          # Calculate Geweke diagnostic test?
                                 diag_save = TRUE,                       # Save the diagnostics as a txt file?
                                 diag_name = "N2O_diagnostics_3-isotopes",              # If yes, specify the base file name (.txt will be appended later)
                                 indiv_effect = indiv_effect,            # Is Individual a random effect in the model? (already specified)
                                 plot_post_save_png = FALSE,             # Save posterior density plots as pngs?
                                 plot_pairs_save_png = FALSE,            # Save pairs plot as png?
                                 plot_xy_save_png = FALSE)               # Save xy/trace plot as png?

output_JAGS(jags.2isotopes, mix, source, output_options.2isotopes)

# END: Model scenario 2: 3 isotopes (δ15N, δ18O, SP), 3 end-members (Freshwater, Marine, Soil)
###



###
# BEGIN: Model scenario 3: 2 isotopes (δ15N, δ18O), 2 end-members (Continental, Marine)

# Save the necessary data as csv files since MixSIAR likes to load from csv files
# Consumer is the tropospheric value we want. It is the Modern antthro+natural value from Röckmann et al. 2003
# Sources are the subsetting Continental (Freshwater, Soils, Urban Wastewater) and Marine data
# Discrimination is isotope fractionation between Sources and Consumer -- set at 0
subContinental <- rbind(subset(subDatasetS1, Category == "Freshwater" | Category == "Soil"), subset(DatasetS1, Category == "Urban Wastewater"))
Continental$Category <- "Continental"
write.csv(subset(globModSol["10",], select = c("d15N", "d18Ovsmow")), 
          file = "N2O_consumer_2-isotopes_2-endmembers.csv", row.names = FALSE)
write.csv(rename(na.omit(subset(subContinental, select = c("Category", "d15N", "d18Ovsmow"))), 
                 c("Category" = "Source")), 
          file = "N2O_sources_2-isotopes_2-endmembers.csv", row.names = FALSE)
write.csv(data.frame(Source = c("Continental", "Marine"),
                     Meand15N = 0, SDd15N = 0, Meand18Ovsmow = 0, SDd18Ovsmow = 0),
          file = "N2O_discrimination_2-isotopes_2-endmembers.csv", row.names = FALSE)


# Load the data from csv files using MixSIAR scripts
mix <- load_mix_data(filename="N2O_consumer_2-isotopes_2-endmembers.csv", iso_names=c("d15N", "d18Ovsmow"), 
                     random_effects=NULL, cont_effects=NULL, fixed_effects=NULL)
source <- load_source_data(filename="N2O_sources_2-isotopes_2-endmembers.csv", 
                           source_factors=NULL, conc_dep=FALSE, data_type="raw", mix)    
discr <- load_discr_data(filename="N2O_discrimination_2-isotopes_2-endmembers.csv", mix)


# Plot data to check it was loaded correctly
plot_data(filename="N2O_isospace_plot_2-isotopes_2-endmembers", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)



# Write JAGS model file
model_filename <- "MixSIAR_model_N2O_2-isotopes_2-endmembers.txt"   # Name of the JAGS model file
indiv_effect <- FALSE                    # Include Individual as a random effect in the model?
nested <- FALSE                         # If there are 2 random effects, is the 2nd nested in the 1st (hierarchical)?
write_JAGS_model(model_filename, indiv_effect, nested, resid_err=TRUE, mix,source)



# Run model
# JAGS output will be saved as 'jags.3isotopes'
#
# MCMC run options:
# run <- "test" # list(chainLength=1000, burn=500, thin=1, chains=3, calcDIC=TRUE)
# run <- "very short" # list(chainLength=10000, burn=5000, thin=5, chains=3, calcDIC=TRUE)
# run <- "short" # list(chainLength=50000, burn=25000, thin=25, chains=3, calcDIC=TRUE)
# run <- "normal" # list(chainLength=100000, burn=50000, thin=50, chains=3, calcDIC=TRUE)
# run <- "long" # list(chainLength=300000, burn=200000, thin=100, chains=3, calcDIC=TRUE)
# run <- "very long" # list(chainLength=1000000, burn=700000, thin=300, chains=3, calcDIC=TRUE)
# run <- "extreme" # list(chainLength=3000000, burn=2700000, thin=300, chains=3, calcDIC=TRUE)

jags.3isotopes <- run_model(run="normal", indiv_effect,mix,source,discr,model_filename)



# Process JAGS output
output_options.3isotopes <- list(summary_save = TRUE,                   # Save the summary statistics as a txt file?
                                 summary_name = "N2O_summary_statistics_2-isotopes_2-endmembers",    # If yes, specify the base file name (.txt will be appended later)
                                 sup_post = FALSE,                       # Suppress posterior density plot output in R?
                                 plot_post_save_pdf = TRUE,            # Save posterior density plots as pdfs?
                                 plot_post_name = "N2O_posterior_density_2-isotopes_2-endmembers",   # If yes, specify the base file name(s) (.pdf/.png will be appended later)
                                 sup_pairs = FALSE,                      # Suppress pairs plot output in R?
                                 plot_pairs_save_pdf = TRUE,            # Save pairs plot as pdf?
                                 plot_pairs_name = "N2O_pairs_plot_2-isotopes_2-endmembers",         # If yes, specify the base file name (.pdf/.png will be appended later)
                                 sup_xy = FALSE,                         # Suppress xy/trace plot output in R?
                                 plot_xy_save_pdf = TRUE,               # Save xy/trace plot as pdf?
                                 plot_xy_name = "N2O_xy_plot_2-isotopes_2-endmembers",               # If yes, specify the base file name (.pdf/.png will be appended later)
                                 gelman = TRUE,                          # Calculate Gelman-Rubin diagnostic test?
                                 heidel = FALSE,                         # Calculate Heidelberg-Welch diagnostic test?
                                 geweke = TRUE,                          # Calculate Geweke diagnostic test?
                                 diag_save = TRUE,                       # Save the diagnostics as a txt file?
                                 diag_name = "N2O_diagnostics_2-isotopes_2-endmembers",              # If yes, specify the base file name (.txt will be appended later)
                                 indiv_effect = indiv_effect,            # Is Individual a random effect in the model? (already specified)
                                 plot_post_save_png = FALSE,             # Save posterior density plots as pngs?
                                 plot_pairs_save_png = FALSE,            # Save pairs plot as png?
                                 plot_xy_save_png = FALSE)               # Save xy/trace plot as png?

output_JAGS(jags.3isotopes, mix, source, output_options.3isotopes)

# END: Model scenario 3: 2 isotopes (δ15N, δ18O), 2 end-members (Continental, Marine)
###



###
# BEGIN: Model scenario 4: 3 isotopes (δ15N, δ18O, SP), 2 end-members (Continental, Marine)

# Save the necessary data as csv files since MixSIAR likes to load from csv files
# Consumer is the tropospheric value we want. It is the Modern antthro+natural value from Röckmann et al. 2003
# Sources are the subsetting Freshwater, Marine, and Soils data
# Discrimination is isotope fractionation between Sources and Consumer -- set at 0
subContinental <- rbind(subDatasetS1, subset(DatasetS1, Category == "Urban Wastewater"))
subContinental$Category <- revalue(subContinental$Category, c("Freshwater" = "Continental", "Soil" = "Continental", "Urban Wastewater"= "Continental"))
write.csv(cbind(subset(globModSol["10",], select = c("d15N", "d18Ovsmow")), "SP" = 12.7), 
          file = "N2O_consumer_3-isotopes_2-endmembers.csv", row.names = FALSE)
write.csv(rename(na.omit(subset(subContinental, select = c("Category", "d15N", "d18Ovsmow", "SP"))), 
                 c("Category" = "Source")), 
          file = "N2O_sources_3-isotopes_2-endmembers.csv", row.names = FALSE)
write.csv(data.frame(Source = c("Continental", "Marine"),
                     Meand15N = 0, SDd15N = 0, Meand18Ovsmow = 0, SDd18Ovsmow = 0, MeanSP = 0, SDSP = 0),
          file = "N2O_discrimination_3-isotopes_2-endmembers.csv", row.names = FALSE)


# Load the data from csv files using MixSIAR scripts
mix <- load_mix_data(filename="N2O_consumer_3-isotopes_2-endmembers.csv", iso_names=c("d15N", "d18Ovsmow", "SP"), 
                     random_effects=NULL, cont_effects=NULL, fixed_effects=NULL)
source <- load_source_data(filename="N2O_sources_3-isotopes_2-endmembers.csv", 
                           source_factors=NULL, conc_dep=FALSE, data_type="raw", mix)    
discr <- load_discr_data(filename="N2O_discrimination_3-isotopes_2-endmembers.csv", mix)


# Plot data to check it was loaded correctly
plot_data(filename="N2O_isospace_plot_3-isotopes_2-endmembers", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)



# Write JAGS model file
model_filename <- "MixSIAR_model_N2O_3-isotopes_2-endmembers.txt"   # Name of the JAGS model file
indiv_effect <- FALSE                    # Include Individual as a random effect in the model?
nested <- FALSE                         # If there are 2 random effects, is the 2nd nested in the 1st (hierarchical)?
write_JAGS_model(model_filename, indiv_effect, nested, resid_err=TRUE, mix,source)



# Run model
# JAGS output will be saved as 'jags.4isotopes'
#
# MCMC run options:
# run <- "test" # list(chainLength=1000, burn=500, thin=1, chains=3, calcDIC=TRUE)
# run <- "very short" # list(chainLength=10000, burn=5000, thin=5, chains=3, calcDIC=TRUE)
# run <- "short" # list(chainLength=50000, burn=25000, thin=25, chains=3, calcDIC=TRUE)
# run <- "normal" # list(chainLength=100000, burn=50000, thin=50, chains=3, calcDIC=TRUE)
# run <- "long" # list(chainLength=300000, burn=200000, thin=100, chains=3, calcDIC=TRUE)
# run <- "very long" # list(chainLength=1000000, burn=700000, thin=300, chains=3, calcDIC=TRUE)
# run <- "extreme" # list(chainLength=3000000, burn=2700000, thin=300, chains=3, calcDIC=TRUE)

jags.4isotopes <- run_model(run="long", indiv_effect,mix,source,discr,model_filename)



# Process JAGS output
output_options.4isotopes <- list(summary_save = TRUE,                   # Save the summary statistics as a txt file?
                                 summary_name = "N2O_summary_statistics_3-isotopes_2-endmembers",    # If yes, specify the base file name (.txt will be appended later)
                                 sup_post = FALSE,                       # Suppress posterior density plot output in R?
                                 plot_post_save_pdf = TRUE,            # Save posterior density plots as pdfs?
                                 plot_post_name = "N2O_posterior_density_3-isotopes_2-endmembers",   # If yes, specify the base file name(s) (.pdf/.png will be appended later)
                                 sup_pairs = FALSE,                      # Suppress pairs plot output in R?
                                 plot_pairs_save_pdf = TRUE,            # Save pairs plot as pdf?
                                 plot_pairs_name = "N2O_pairs_plot_3-isotopes_2-endmembers",         # If yes, specify the base file name (.pdf/.png will be appended later)
                                 sup_xy = FALSE,                         # Suppress xy/trace plot output in R?
                                 plot_xy_save_pdf = TRUE,               # Save xy/trace plot as pdf?
                                 plot_xy_name = "N2O_xy_plot_3-isotopes_2-endmembers",               # If yes, specify the base file name (.pdf/.png will be appended later)
                                 gelman = TRUE,                          # Calculate Gelman-Rubin diagnostic test?
                                 heidel = FALSE,                         # Calculate Heidelberg-Welch diagnostic test?
                                 geweke = TRUE,                          # Calculate Geweke diagnostic test?
                                 diag_save = TRUE,                       # Save the diagnostics as a txt file?
                                 diag_name = "N2O_diagnostics_3-isotopes_2-endmembers",              # If yes, specify the base file name (.txt will be appended later)
                                 indiv_effect = indiv_effect,            # Is Individual a random effect in the model? (already specified)
                                 plot_post_save_png = FALSE,             # Save posterior density plots as pngs?
                                 plot_pairs_save_png = FALSE,            # Save pairs plot as png?
                                 plot_xy_save_png = FALSE)               # Save xy/trace plot as png?

output_JAGS(jags.4isotopes, mix, source, output_options.4isotopes)

# END:  Model scenario 4: 3 isotopes (δ15N, δ18O, SP), 2 end-members (Continental, Marine)
###

# EOF