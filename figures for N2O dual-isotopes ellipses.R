# Figures for N2O Ellipse dual-isotope paper: Snider, Venkiteswaran, Schiff, and Spoelstra
# https://github.com/jjvenky/Global-N2O-Ellipses
#


###
# Load packages
library(ggplot2)
library(reshape)
library(RColorBrewer)
library(scales)
library(grid)
library(siar)
library(plyr)
###



###
# Load the data
DatasetS1 <- read.csv('Dataset S1.csv', head = TRUE, skip = 1)
names(DatasetS1) <- c("Reference", "Category", "d15N", "SP", "d18Ovsmow", "Notes", "Filter")
# Remove any blank lines
DatasetS1 <- subset.data.frame(DatasetS1, subset = DatasetS1$Category != "")
# Check that everything looks good
str(DatasetS1)
# Will need Freshwater & Soil for a "Continental" calculation
Continental <- subset.data.frame(DatasetS1, subset = (Category == "Freshwater" | Category == "Soil") & Filter == "yes")
Continental$Category <- "Continental"
###



###
# Check if data look good
# Basic figure -- have to reorganise to plot with ellipses
ggplot(data = DatasetS1, aes(x = d15N, y = d18Ovsmow, colour = Category)) + 
  geom_point(size = 2) + 
  labs(x = expression(paste("\u03b4" ^15, "N-N" [2], "O (\u2030 air N" [2], ")")), 
       y = expression(paste("\u03b4" ^18, "O-N" [2], "O (\u2030 VSMOW)")) ) + 
  scale_color_brewer(palette = "Dark2")
###



###
# Calculate ellipses metrics
#
### BEGIN siber code
# code modified from Jackson & Parnell
# http://www.tcd.ie/Zoology/research/research/theoretical/Rpodcasts.php#siber
# siber package is presented in Parnell, A.C., Inger R., & Bearhop, S. 2011.
# Comparing isotopic niche widths among and within communities: SIBER – Stable
# Isotope Bayesian Ellipses in R. Journal of Animal Ecology, 80, 595-602.
# doi: http://dx.doi.org/10.1111/j.1365-2656.2011.01806.x
#
# Loop through the data by Category and calculate the ellipses and other metrics
# then reassemble the results and fix the factors
attach(DatasetS1)
ngroups <- length(unique(Category))
# split the isotope data based on group
spx <- split(d15N, Category)
spy <- split(d18Ovsmow, Category)
# create some empty vectors for recording our metrics
SEA <- numeric(ngroups)
SEAc <- numeric(ngroups)
TA <- numeric(ngroups)
axisA <- numeric(ngroups)
axisB <- numeric(ngroups)
theta <- numeric(ngroups)
ellipses <- data.frame(d15N = NA, d18Ovsmow = NA, set = NA)
for (j in as.numeric(unique(Category))){
  # Fit a standard ellipse to the data
  SE <- standard.ellipse(spx[[j]], spy[[j]], steps = 1)
  # Extract the estimated SEA and SEAc from this object
  SEA[j] <- SE$SEA
  SEAc[j] <- SE$SEAc
  # Extract the axes (a and b) and angle (theta)
  axisA[j] <- SE$a
  axisB[j] <- SE$bc
  theta[j] <- SE$theta
  # The ellipses themselves
  ellipses <- rbind(ellipses, data.frame(d15N = SE$xSEAc, d18Ovsmow = SE$ySEAc, set = j))
  # Calculate the convex hull for the jth group's isotope values
  # held in the objects created using split() called spx and spy
  CH <- convexhull(spx[[j]], spy[[j]])
  # Extract the area of the convex hull from this object
  TA[j] <- CH$TA
}
detach(DatasetS1)
# Stich together the metrics
ellipseMetrics <- data.frame(Category = levels(DatasetS1$Category), SEA, SEAc, TA, axisA, axisB, theta, slope = tan(theta))
# Fix Category in the ellipses
ellipses$Category <- factor(ellipses$set, labels = c("Antarctic", "Freshwater", "Groundwater", "Marine", "Soil", "Stratosphere", "Troposphere", "Urban Wastewater"))
### END siber code



###
# TABLE 1: Descriptive metrics for the ellipses
table1 <- merge(ellipseMetrics, ddply(DatasetS1, .(Category), summarise,
                                      N = length(d15N),
                                      meand15N = mean(d15N),
                                      sd15N = sd(d15N),
                                      meand18O = mean(d18Ovsmow),
                                      sd18O = sd(d18Ovsmow),
                                      cor = cor(d15N, d18Ovsmow)),
                by = "Category"
)
print(table1, digits = 4)
#           Category       SEA     SEAc        TA   axisA   axisB  theta  slope         Category   N meand15N  sd15N meand18O   sd18O    cor
# 1        Antarctic 2994.8576 3085.611  7820.250 34.6752 27.9055 0.8609 1.1636        Antarctic  35  -40.843 30.749    29.03 31.8224 0.2256
# 2       Freshwater  202.6786  202.954  2058.840 12.0565  5.3547 0.7016 0.8450       Freshwater 738   -4.649  9.835    41.77  8.7897 0.6656
# 3      Groundwater  767.0175  768.470 10449.040 20.1799 12.1101 0.9326 1.3481      Groundwater 530  -13.967 15.459    45.34 17.7381 0.4552
# 4           Marine   91.6201   91.806  1463.705  9.7072  3.0074 1.3752 5.0474           Marine 495    6.629  3.499    47.35  9.5400 0.4866
# 5             Soil  296.0861  296.422  2694.265 14.0310  6.7209 0.6300 0.7292             Soil 884  -14.849 12.008    31.23  9.8891 0.6083
# 6     Stratosphere   43.3487   43.500   584.720 27.7854  0.4975 0.7256 0.8870     Stratosphere 288   20.315 20.789    56.39 18.4413 0.9994
# 7      Troposphere    0.4729    0.475     3.565  0.5009  0.3012 0.4318 0.4608      Troposphere 225    6.546  0.472    44.40  0.3441 0.3758
# 8 Urban Wastewater  539.4780  545.472  2878.350 15.3783 11.2283 0.9615 1.4330 Urban Wastewater  92  -11.564 12.701    31.51 14.1378 0.2922
# write.csv(table1, file = "table1.csv")
###


###
# Basic figure with ellipses to test if everything looks good and makes sense
# Add transparency to the points to make it easier to see the density of points
ggplot(mapping = aes(x = d15N, y = d18Ovsmow, colour = Category)) + 
  geom_point(data = DatasetS1, size = 2, alpha = 1/3) + 
  geom_path(data = ellipses, size = 1) + 
  labs(x = expression(paste("\u03b4" ^15, "N-N" [2], "O (\u2030 air N" [2], ")")), 
       y = expression(paste("\u03b4" ^18, "O-N" [2], "O (\u2030 VSMOW)")) ) + 
  scale_color_brewer(palette = "Dark2")
###


### ---------- Figures for paper below ----------


###
# Will make separate variables for the colours in each figure
# Need to subset some of the data to make plotting easier
allButST <- subset(DatasetS1, Category != "Stratosphere" & Category != "Tropposphere")
allButST.ell <- subset(ellipses, Category != "Stratosphere" & Category != "Troposphere")
Stratosphere <- subset(DatasetS1, Category == "Stratosphere")
Troposphere <- data.frame(d15N = 6.9, d18Ovsmow = 44.5, Category = "Troposphere")
StratoTropo <-  subset(DatasetS1, Category == "Stratosphere" | Category == "Troposphere")
###



###
# THIS WILL BE FIGURE 1
# Save at 5.00" x 6.83" and base_size = 12 and Arial for PLOS ONE
cols.Fig1 <- brewer_pal(pal = "Dark2")(3)
cols.Fig1 <- c(cols.Fig1[c(1, 3, 2, 1, 2)], "grey", "black", "black")
ggplot(mapping = aes(x = d15N, y = d18Ovsmow, linetype = Category, 
                     shape = Category, colour = Category)) +
  geom_point(data = allButST, size = 2, alpha = 1/3, aes(colour = Category, shape = Category)) +
  geom_point(data = Stratosphere, size = 2, aes(colour = Category, shape = Category)) +
  geom_point(data = Troposphere, size = 2, aes(colour = Category, shape = Category)) +
  geom_path(data = allButST.ell, size = 1, aes(linetype = Category, colour = Category)) +
  geom_path(data = allButST.ell, size = 0.2, aes(linetype = Category), colour = "black") +
  scale_x_continuous(limits = c(-100, 50)) +
  scale_y_continuous(limits = c(-5, 100)) +
  scale_colour_manual(values = cols.Fig1) +
  scale_shape_manual(values = c(16, 16, 16, 17, 17, 3, 8, 17)) +
  scale_linetype_manual(values = c(1, 1, 1, 2, 2, 0, 0, 2)) +
  labs(x = expression(paste("\u03b4" ^15, "N-N" [2], "O (\u2030 air N" [2], ")")), 
       y = expression(paste("\u03b4" ^18, "O-N" [2], "O (\u2030 VSMOW)")), 
       linetype = "", shape = "", colour = "", fill = "") +
  coord_fixed() +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(legend.position = c(.15, 0.80), legend.background = element_rect(fill = NA), 
        legend.margin = unit(0.5, "cm"), legend.key.size = unit(0.5, "cm"), 
        legend.key = element_rect(colour = NA),
        plot.margin = unit(c(0,0,0,0), "mm"))
###



###
# THIS WILL BE FIGIURE 2
# Figure generated in different software
# This is not a complete version
UW <- subset(DatasetS1, Category == "Urban Wastewater")
UW.ell <- subset(ellipses, Category == "Urban Wastewater")
UW.means <- data.frame(d15N = mean(UW$d15N), d18Ovsmow = mean(UW$d18Ovsmow))
UW.lm1 <- lm(d18Ovsmow~d15N, data = UW)
UW.lm1coefs <- data.frame(a = coef(UW.lm1)[1], b = coef(UW.lm1)[2])
UW.lm2 <- lm(d15N~d18Ovsmow, data = UW)
UW.lm2coefs <- data.frame(a = coef(UW.lm2)[1], b = coef(UW.lm2)[2])
UW.ell.line <- data.frame(a = (UW.means$d18Ovsmow - 1.66*UW.means$d15N), b = 1.66)

ggplot(mapping = aes(x = d15N, y = d18Ovsmow)) + geom_point(data = UW) +
  geom_path(data = UW.ell, size = 2) +
  geom_ribbon(data = UW.ell, aes(xmin=min(d15N), xmax=max(d15N), ymin=min(d18Ovsmow), ymax=max(d18Ovsmow)), 
              fill = NA, colour = "black", size = 2) +
  geom_abline(data = UW.lm1coefs, aes(intercept = a, slope = b)) +
  geom_abline(data = UW.lm2coefs, aes(intercept = -a/b, slope = 1/b)) +
  geom_abline(data = UW.ell.line, aes(intercept = a, slope = b)) +
  geom_point(data = UW.means, colour = "black", size = 6) +
  labs(x = expression(paste("\u03b4" ^15, "N-N" [2], "O (\u2030 air N" [2], ")")), y = expression(paste("\u03b4" ^18, "O-N" [2], "O (\u2030 VSMOW)"))) +
  coord_fixed() +
  theme_classic()
###



###
# Will make separate variables for the colours in each figure
# Need to subset some of the data to make plotting easier
FMSUW <- subset(DatasetS1, Category == "Freshwater" | Category == "Marine" | Category == "Soil" | Category == "Urban Wastewater")
FMSUW.ell <- subset(ellipses, Category == "Freshwater" | Category == "Marine" | Category == "Soil" | Category == "Urban Wastewater")
###



###
# THIS WILL BE FIGURE 3
# It is the same as Figure 1 but with narrow axes and without Antarctic and Groundwater
# Numbers by Category: ddply(subset(DatasetS1, d15N >= -40 & d15N <= 20 & d18Ovsmow >= 10 & d18Ovsmow <= 70), .(Category), summarise, N = length(d15N))
# Save at 6.83" x 6.83" and base_size = 12 and Arial for PLOS ONE
cols.Fig3 <- brewer_pal(pal = "Dark2")(3)
cols.Fig3 <- c(cols.Fig3[c(3, 1, 2)], "grey", "black", "black")
ggplot(mapping = aes(x = d15N, y = d18Ovsmow, linetype = Category, 
                     shape = Category, colour = Category)) +
  geom_point(data = FMSUW, size = 2, alpha = 1/3, aes(colour = Category, shape = Category)) +
  geom_point(data = Stratosphere, size = 2, aes(colour = Category, shape = Category)) +
  geom_point(data = Troposphere, size = 2, aes(colour = Category, shape = Category)) +
  geom_path(data = FMSUW.ell, size = 1, aes(linetype = Category, colour = Category)) +
  geom_path(data = FMSUW.ell, size = 0.2, aes(linetype = Category), colour = "black") +
  scale_x_continuous(limits = c(-40, 20)) +
  scale_y_continuous(limits = c(10, 70)) +
  scale_colour_manual(values = cols.Fig3) +
  scale_shape_manual(values = c(16, 17, 15, 3, 8, 17)) +
  scale_linetype_manual(values = c(1, 2, 2, 0, 0, 2)) +
  labs(x = expression(paste("\u03b4" ^15, "N-N" [2], "O (\u2030 air N" [2], ")")), 
       y = expression(paste("\u03b4" ^18, "O-N" [2], "O (\u2030 VSMOW)")), 
       linetype = "", shape = "", colour = "", fill = "") +
  coord_fixed() +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(legend.position = c(.20, 0.88), legend.background = element_rect(fill = NA), 
        legend.margin = unit(0.5, "cm"), legend.key.size = unit(0.5, "cm"), 
        legend.key = element_rect(colour = NA),
        plot.margin = unit(c(0,0,0,0), "mm"))
###



###
# THESE WILL BE FIGURE 4A AND 4B
# It is the same as Figure 1 but with narrow axes and without Antarctic and Groundwater
# Site Preference Figures
# Save at 3.27" x 3.27" and base_size = 12 and Arial for PLOS ONE

# New to make ellipses for each SP plot
DatasetS1SP <- subset.data.frame(DatasetS1, subset = is.na(SP) == FALSE)
### BEGIN siber code
attach(DatasetS1SP)
ngroups <- length(levels(Category))
# split the isotope data based on group
spx <- split(d15N, Category)
spy <- split(SP, Category)
# create some empty vectors for recording our metrics
SEA <- numeric(ngroups)
SEAc <- numeric(ngroups)
TA <- numeric(ngroups)
axisA <- numeric(ngroups)
axisB <- numeric(ngroups)
theta <- numeric(ngroups)
ellipsesSPN <- data.frame(d15N = NA, SP = NA, set = NA)
for (j in as.numeric(unique(Category))){
  SE <- standard.ellipse(spx[[j]],spy[[j]],steps = 1)
  SEA[j] <- SE$SEA
  SEAc[j] <- SE$SEAc
  axisA[j] <- SE$a
  axisB[j] <- SE$bc
  theta[j] <- SE$theta
  ellipsesSPN <- rbind(ellipsesSPN, data.frame(d15N = SE$xSEAc, SP = SE$ySEAc, set = j))
  CH <- convexhull(spx[[j]],spy[[j]])
  TA[j] <- CH$TA
}
detach(DatasetS1SP)
ellipseMetricsSPN <- data.frame(Category = levels(DatasetS1SP$Category), SEA, SEAc, TA, axisA, axisB, theta, slope = tan(theta))
# Fix Category in the ellipses
ellipsesSPN$Category <- factor(ellipsesSPN$set, labels = c("Antarctic", "Freshwater", "Groundwater", "Marine", "Soil", "Stratosphere", "Troposphere", "Urban Wastewater"))
### END siber code

### BEGIN siber code
attach(DatasetS1SP)
ngroups <- length(levels(Category))
# split the isotope data based on group
spx <- split(d18Ovsmow, Category)
spy <- split(SP, Category)
# create some empty vectors for recording our metrics
SEA <- numeric(ngroups)
SEAc <- numeric(ngroups)
TA <- numeric(ngroups)
axisA <- numeric(ngroups)
axisB <- numeric(ngroups)
theta <- numeric(ngroups)
ellipsesSPO <- data.frame(d18Ovsmow = NA, SP = NA, set = NA)
for (j in as.numeric(unique(Category))){
  SE <- standard.ellipse(spx[[j]],spy[[j]],steps = 1)
  SEA[j] <- SE$SEA
  SEAc[j] <- SE$SEAc
  axisA[j] <- SE$a
  axisB[j] <- SE$bc
  theta[j] <- SE$theta
  ellipsesSPO <- rbind(ellipsesSPO, data.frame(d18Ovsmow = SE$xSEAc, SP = SE$ySEAc, set = j))
  CH <- convexhull(spx[[j]],spy[[j]])
  TA[j] <- CH$TA
}
detach(DatasetS1SP)
ellipseMetricsSPO <- data.frame(Category = levels(DatasetS1SP$Category), SEA, SEAc, TA, axisA, axisB, theta, slope = tan(theta))
# Fix Category in the ellipses
ellipsesSPO$Category <- factor(ellipsesSPO$set, labels = c("Antarctic", "Freshwater", "Groundwater", "Marine", "Soil", "Stratosphere", "Troposphere", "Urban Wastewater"))
### END siber code



# THESE WILL BE FIGURE 4A AND 4B
# It is the same as Figure 1 but with narrow axes and without Antarctic and Groundwater
# Site Preference Figures
# Save at 3.27" x 3.27" and base_size = 8 and Arial for PLOS ONE
FMSUW.SPN.ell <- subset(ellipsesSPN, Category == "Freshwater" | Category == "Marine" | Category == "Soil" | Category == "Urban Wastewater")
FMSUW.SPO.ell <- subset(ellipsesSPO, Category == "Freshwater" | Category == "Marine" | Category == "Soil" | Category == "Urban Wastewater")
FMSUW.SP <- subset(DatasetS1SP, Category == "Freshwater" | Category == "Marine" | Category == "Soil" | Category == "Urban Wastewater")

cols.Fig4 <- brewer_pal(pal = "Dark2")(3)
cols.Fig4 <- c(cols.Fig4[c(3, 1, 2)], "grey", "black", "black")
ggplot(mapping = aes(x = d15N, y = SP, linetype = Category, 
                     shape = Category, colour = Category)) +
  geom_point(data = FMSUW.SP, size = 1, alpha = 1/2, aes(colour = Category, shape = Category)) +
  geom_point(data = StratoTropo, size = 1, aes(colour = Category, shape = Category)) +
  geom_path(data = FMSUW.SPN.ell, size = 1/2, aes(linetype = Category, colour = Category)) +
  geom_path(data = FMSUW.SPN.ell, size = 0.2, aes(linetype = Category), colour = "black") +
  scale_x_continuous(limits = c(-40, 20)) +
  scale_y_continuous(limits = c(-10, 50)) +
  scale_colour_manual(values = cols.Fig4) +
  scale_shape_manual(values = c(16, 17, 15, 3, 8, 17)) +
  scale_linetype_manual(values = c(1, 2, 2, 0, 0, 2)) +
  labs(x = expression(paste("\u03b4" ^15, "N-N" [2], "O (\u2030 air N" [2], ")")), 
       y = expression(paste("\u03b4" ^15, "N" ^{α}, "-N" [2], "O  \u2212  \u03b4" ^15, "N" ^{β}, "-N" [2], "O (\u2030 air N" [2], ")")), 
       linetype = "", shape = "", colour = "", fill = "") +
  coord_fixed() +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(legend.position = c(.20, 0.88), legend.background = element_rect(fill = NA), 
        legend.margin = unit(0.5, "cm"), legend.key.size = unit(0.25, "cm"), 
        legend.key = element_rect(colour = NA),
        plot.margin = unit(c(0,0,0,0), "mm"))

ggplot(mapping = aes(x = d18Ovsmow, y = SP, linetype = Category, 
                     shape = Category, colour = Category)) +
  geom_point(data = FMSUW.SP, size = 1, alpha = 1/2, aes(colour = Category, shape = Category)) +
  geom_point(data = StratoTropo, size = 1, aes(colour = Category, shape = Category)) +
  geom_path(data = FMSUW.SPO.ell, size = 1/2, aes(linetype = Category, colour = Category)) +
  geom_path(data = FMSUW.SPO.ell, size = 0.2, aes(linetype = Category), colour = "black") +
  scale_x_continuous(limits = c(10, 70)) +
  scale_y_continuous(limits = c(-10, 50)) +
  scale_colour_manual(values = cols.Fig4) +
  scale_shape_manual(values = c(16, 17, 15, 3, 8, 17)) +
  scale_linetype_manual(values = c(1, 2, 2, 0, 0, 2)) +
  labs(x = expression(paste("\u03b4" ^18, "O-N" [2], "O (\u2030 VSMOW)")), 
       y = expression(paste("\u03b4" ^15, "N" ^{α}, "-N" [2], "O  \u2212  \u03b4" ^15, "N" ^{β}, "-N" [2], "O (\u2030 air N" [2], ")")), 
       linetype = "", shape = "", colour = "", fill = "") +
  coord_fixed() +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(legend.position = c(.20, 0.88), legend.background = element_rect(fill = NA), 
        legend.margin = unit(0.5, "cm"), legend.key.size = unit(0.25, "cm"), 
        legend.key = element_rect(colour = NA),
        plot.margin = unit(c(0,0,0,0), "mm"))
###



###
# Before Figure 5, have to subset the data and run it through siber to get the ellipses
subDatasetS1 <- subset.data.frame(DatasetS1, subset = Filter == 'yes')
subFMS <- subset(subDatasetS1, Category == "Freshwater" | Category == "Marine" | Category == "Soil")
### BEGIN siber code
attach(subDatasetS1)
ngroups <- length(levels(Category))
# split the isotope data based on group
spx <- split(d15N, Category)
spy <- split(d18Ovsmow, Category)
# create some empty vectors for recording our metrics
SEA <- numeric(ngroups)
SEAc <- numeric(ngroups)
TA <- numeric(ngroups)
axisA <- numeric(ngroups)
axisB <- numeric(ngroups)
theta <- numeric(ngroups)
subellipses <- data.frame(d15N = NA, d18Ovsmow = NA, set = NA)
for (j in as.numeric(unique(Category))){
  SE <- standard.ellipse(spx[[j]],spy[[j]],steps = 1)
  SEA[j] <- SE$SEA
  SEAc[j] <- SE$SEAc
  axisA[j] <- SE$a
  axisB[j] <- SE$bc
  theta[j] <- SE$theta
  subellipses <- rbind(subellipses, data.frame(d15N = SE$xSEAc, d18Ovsmow = SE$ySEAc, set = j))
  CH <- convexhull(spx[[j]],spy[[j]])
  TA[j] <- CH$TA
}
detach(subDatasetS1)
subellipseMetrics <- data.frame(Category = levels(subDatasetS1$Category), SEA, SEAc, TA, axisA, axisB, theta, slope = tan(theta))
# Fix Category in the ellipses
subellipses$Category <- factor(subellipses$set, labels = c("Freshwater", "Marine", "Soil"))
### END siber code
###



###
# Need Continental values for Table 2
# Code for Continental calculations
# Continental created from subset above
ellipsesContinental <- data.frame(d15N = NA, d18Ovsmow = NA)
# Fit a standard ellipse to the data
SE <- standard.ellipse(Continental$d15N, Continental$d18Ovsmow, steps = 1)
SEA <- SE$SEA
SEAc <- SE$SEAc
axisA <- SE$a
axisB <- SE$bc
theta <- SE$theta
ellipsesContinental <- rbind(ellipsesContinental, data.frame(d15N = SE$xSEAc, d18Ovsmow = SE$ySEAc))
CH <- convexhull(Continental$d15N, Continental$d18Ovsmow)
TA <- CH$TA
# Stich together the metrics
ellipseMetricsContinental <- data.frame(Category = "Continental", SEA, SEAc, TA, axisA, axisB, theta, slope = tan(theta))
###


###
# TABLE 2: Descriptive metrics for the ellipses
table2a <- merge(subellipseMetrics, ddply(subDatasetS1, .(Category), summarise,
                                          N = length(d15N),
                                          meand15N = mean(d15N),
                                          sd15N = sd(d15N),
                                          meand18O = mean(d18Ovsmow),
                                          sd18O = sd(d18Ovsmow),
                                          cor = cor(d15N, d18Ovsmow)),
                 by="Category"
)
table2b <- merge(ellipseMetricsContinental, ddply(Continental, .(Category), summarise,
                                                  N = length(d15N),
                                                  meand15N = mean(d15N),
                                                  sd15N = sd(d15N),
                                                  meand18O = mean(d18Ovsmow),
                                                  sd18O = sd(d18Ovsmow),
                                                  cor = cor(d15N, d18Ovsmow)),
                 by = "Category"
)
table2 <- rbind(table2a, table2b)
print(table2, digits = 4)
#     Category    SEA   SEAc      TA  axisA axisB  theta   slope    N meand15N  sd15N meand18O  sd18O     cor
# 1  Freshwater 214.97 215.37 2014.65 12.546 5.459 0.7784  0.9861  527   -7.777  9.719    40.75  9.627 0.68210
# 2      Marine  21.95  22.32   78.61  3.623 1.945 1.5384 30.8679   62    5.139  1.931    44.76  3.621 0.04351
# 3        Soil 287.50 287.86 2694.27 13.028 7.029 0.6442  0.7512  794  -16.657 11.239    30.05  9.632 0.53409
# 4 Continental 298.50 298.73 2944.61 14.474 6.567 0.7484  0.9285 1321  -13.114 11.509    34.32 10.960 0.65772
# write.csv(table2, file = "table2.csv")
###



###
# Will make separate variables for the colours in each figure
# Need to subset some of the data to make plotting easier
subFMS.ell <- subset(subellipses, Category == "Freshwater" | Category == "Marine" | Category == "Soil")
###



###
# THIS WILL BE FIGURE 5
# Save at 3.27" x 3.27" and base_size = 8 and Arial for PLOS ONE
# Need the Supplementary Data Table 2
# Slightly larger axes by 10‰ in order to fit in error bars!
DatasetS2 <- read.csv('Dataset S2.csv', head = TRUE, na.strings=c("n/a"))
names(DatasetS2) <- c("Reference", "Category", "MeanFlux", "sdMeanFlux", "MedianFlux", "MinFlux", "MaxFlux", "d15N", "sdd15N", "SP", "sdSP", "d18Ovsmow", "sdd18Ovsmow", "n", "ErrorInfo", "Site")

cols.Fig5 <- brewer_pal(pal = "Dark2")(3)
cols.Fig5 <- c(cols.Fig5[c(3, 1, 2)], "grey","black","black")
ggplot(mapping = aes(x = d15N, y = d18Ovsmow, linetype = Category, shape = Category, colour = Category)) +
  geom_point(data = Stratosphere, size = 2, aes(colour = Category, shape = Category)) +
  geom_point(data = Troposphere, size = 2, aes(colour = Category, shape = Category)) +  
#   geom_path(data = subFMS.ell, size = 1, aes(colour = Category)) +
  geom_point(data = DatasetS2, size = 2, aes(colour = Category, shape = Category)) +
  geom_errorbar(data = DatasetS2, aes(x = d15N, ymin = d18Ovsmow - sdd18Ovsmow, ymax = d18Ovsmow + sdd18Ovsmow, width=0)) + 
  geom_errorbarh(data = DatasetS2, aes(xmin = d15N - sdd15N, xmax = d15N + sdd15N, height=0)) +
  scale_x_continuous(limits = c(-55,15)) +
  scale_y_continuous(limits = c(5,75)) +
  scale_colour_manual(values = cols.Fig5) +
  scale_shape_manual(values = c(16,17,15,3,8,17)) +
  scale_linetype_manual(values = c(1,1,1,0,0,0)) +
  labs(x = expression(paste("\u03b4" ^15, "N-N" [2], "O (\u2030 air N" [2], ")")), 
       y = expression(paste("\u03b4" ^18, "O-N" [2], "O (\u2030 VSMOW)")), 
       linetype = "", shape = "", colour = "") +
  coord_fixed() +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(legend.position = c(.22,0.88), legend.background = element_rect(), 
        legend.margin = unit(0.5,"cm"), legend.key.size = unit(0.25, "cm"), 
        legend.key = element_rect(colour = NA),
        plot.margin = unit(c(0,0,0,0), "mm"))
###

  
###
# THIS WILL BE FIGURE 6A (filtered data)
# Save at 3.27" x 3.27" and base_size = 8 and Arial for PLOS ONE
cols.Fig6A <- brewer_pal(pal = "Dark2")(3)
cols.Fig6A <- c(cols.Fig6A[c(3, 1, 2)], "grey","black")
ggplot(mapping = aes(x = d15N, y = d18Ovsmow, linetype = Category, shape = Category, colour = Category)) +
  geom_point(data = subFMS, size = 1, alpha = 1/3, aes(colour = Category, shape = Category)) +
  geom_point(data = Stratosphere, size = 1, aes(colour = Category, shape = Category)) +
  geom_point(data = Troposphere, size = 1, aes(colour = Category, shape = Category)) +
  geom_path(data = subFMS.ell, size = 1/2, aes(colour = Category)) +
  scale_x_continuous(limits = c(-40, 20)) +
  scale_y_continuous(limits = c(10, 70)) +
  scale_colour_manual(values = cols.Fig6A) +
  scale_shape_manual(values = c(16, 17, 15, 3, 8)) +
  scale_linetype_manual(values = c(1, 1, 1, 0, 0, 6)) +
  labs(x = expression(paste("\u03b4" ^15, "N-N" [2], "O (\u2030 air N" [2], ")")), 
       y = expression(paste("\u03b4" ^18, "O-N" [2], "O (\u2030 VSMOW)")), 
       linetype = "", shape = "", colour = "") +
  coord_fixed() +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(legend.position = c(.22, 0.88), legend.background = element_rect(fill = NA), 
        legend.margin = unit(0.5, "cm"), legend.key.size = unit(0.25, "cm"), 
        legend.key = element_rect(colour = NA),
        plot.margin = unit(c(0,0,0,0), "mm"))
###



###
# For Site Preference figures like Figure 6A
subDatasetS1SP <- subset.data.frame(DatasetS1SP, 
                                    subset = Filter == 'yes' & is.na(SP) == FALSE & 
                                      (Category == "Freshwater" | Category == "Marine" | Category == "Soil"))
### BEGIN siber code
attach(subDatasetS1SP)
ngroups <- length(levels(Category))
# split the isotope data based on group
spx <- split(d15N, Category)
spy <- split(SP, Category)
# create some empty vectors for recording our metrics
SEA <- numeric(ngroups)
SEAc <- numeric(ngroups)
TA <- numeric(ngroups)
axisA <- numeric(ngroups)
axisB <- numeric(ngroups)
theta <- numeric(ngroups)
subellipsesSPN <- data.frame(d15N = NA, SP = NA, set = NA)
for (j in as.numeric(unique(Category))){
  SE <- standard.ellipse(spx[[j]],spy[[j]],steps = 1)
  SEA[j] <- SE$SEA
  SEAc[j] <- SE$SEAc
  axisA[j] <- SE$a
  axisB[j] <- SE$bc
  theta[j] <- SE$theta
  subellipsesSPN <- rbind(subellipsesSPN, data.frame(d15N = SE$xSEAc, SP = SE$ySEAc, set = j))
  CH <- convexhull(spx[[j]],spy[[j]])
  TA[j] <- CH$TA
}
detach(subDatasetS1SP)
ellipseMetricsSPN <- data.frame(Category = levels(subDatasetS1SP$Category), SEA, SEAc, TA, axisA, axisB, theta, slope = tan(theta))
# Fix Category in the ellipses
subellipsesSPN$Category <- factor(subellipsesSPN$set, labels = c("Freshwater", "Marine", "Soil"))
### END siber code

### BEGIN siber code
attach(subDatasetS1SP)
ngroups <- length(levels(Category))
# split the isotope data based on group
spx <- split(d18Ovsmow, Category)
spy <- split(SP, Category)
# create some empty vectors for recording our metrics
SEA <- numeric(ngroups)
SEAc <- numeric(ngroups)
TA <- numeric(ngroups)
axisA <- numeric(ngroups)
axisB <- numeric(ngroups)
theta <- numeric(ngroups)
subellipsesSPO <- data.frame(d18Ovsmow = NA, SP = NA, set = NA)
for (j in as.numeric(unique(Category))){
  SE <- standard.ellipse(spx[[j]],spy[[j]],steps = 1)
  SEA[j] <- SE$SEA
  SEAc[j] <- SE$SEAc
  axisA[j] <- SE$a
  axisB[j] <- SE$bc
  theta[j] <- SE$theta
  subellipsesSPO <- rbind(subellipsesSPO, data.frame(d18Ovsmow = SE$xSEAc, SP = SE$ySEAc, set = j))
  CH <- convexhull(spx[[j]],spy[[j]])
  TA[j] <- CH$TA
}
detach(subDatasetS1SP)
ellipseMetricsSPO <- data.frame(Category = levels(subDatasetS1SP$Category), SEA, SEAc, TA, axisA, axisB, theta, slope = tan(theta))
# Fix Category in the ellipses
subellipsesSPO$Category <- factor(subellipsesSPO$set, labels = c("Freshwater", "Marine", "Soil"))
### END siber code



###
# THESE WILL BE FIGURES 6B and 6C
# It is the same as Figure 6A but SP
# Save at 3.27" x 3.27" and base_size = 8 and Arial for PLOS ONE
ggplot(mapping = aes(x = d15N, y = SP, linetype = Category, 
                     shape = Category, colour = Category)) +
  geom_point(data = subDatasetS1SP, size = 1, alpha = 1/2, aes(colour = Category, shape = Category)) +
  geom_point(data = StratoTropo, size = 1, aes(colour = Category, shape = Category)) +
  geom_path(data = subellipsesSPN, size = 1/2, aes(linetype = Category, colour = Category)) +
  scale_x_continuous(limits = c(-40, 20)) +
  scale_y_continuous(limits = c(-10, 50)) +
  scale_colour_manual(values = cols.Fig6A) +
  scale_shape_manual(values = c(16, 17, 15, 3, 8)) +
  scale_linetype_manual(values = c(1, 1, 1, 0, 0, 6)) +
  labs(x = expression(paste("\u03b4" ^15, "N-N" [2], "O (\u2030 air N" [2], ")")), 
       y = expression(paste("\u03b4" ^15, "N" ^{α}, "-N" [2], "O  \u2212  \u03b4" ^15, "N" ^{β}, "-N" [2], "O (\u2030 air N" [2], ")")), 
       linetype = "", shape = "", colour = "", fill = "") +
  coord_fixed() +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(legend.position = c(.22, 0.88), legend.background = element_rect(fill = NA), 
        legend.margin = unit(0.5, "cm"), legend.key.size = unit(0.25, "cm"), 
        legend.key = element_rect(colour = NA),
        plot.margin = unit(c(0,0,0,0), "mm"))

ggplot(mapping = aes(x = d18Ovsmow, y = SP, linetype = Category, 
                     shape = Category, colour = Category)) +
  geom_point(data = subDatasetS1SP, size = 1, alpha = 1/2, aes(colour = Category, shape = Category)) +
  geom_point(data = StratoTropo, size = 1, aes(colour = Category, shape = Category)) +
  geom_path(data = subellipsesSPO, size = 1/2, aes(linetype = Category, colour = Category)) +
  scale_x_continuous(limits = c(10, 70)) +
  scale_y_continuous(limits = c(-10, 50)) +
  scale_colour_manual(values = cols.Fig6A) +
  scale_shape_manual(values = c(16, 17, 15, 3, 8)) +
  scale_linetype_manual(values = c(1, 1, 1, 0, 0, 6)) +
  labs(x = expression(paste("\u03b4" ^18, "O-N" [2], "O (\u2030 VSMOW)")), 
       y = expression(paste("\u03b4" ^15, "N" ^{α}, "-N" [2], "O  \u2212  \u03b4" ^15, "N" ^{β}, "-N" [2], "O (\u2030 air N" [2], ")")), 
       linetype = "", shape = "", colour = "", fill = "") +
  coord_fixed() +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(legend.position = c(.22, 0.88), legend.background = element_rect(fill = NA), 
        legend.margin = unit(0.5, "cm"), legend.key.size = unit(0.25, "cm"), 
        legend.key = element_rect(colour = NA),
        plot.margin = unit(c(0,0,0,0), "mm"))
###


###
# Add the Global Model Solution data from literature values
globModSol <- read.csv('Global Model Solutions.csv')
names(globModSol) <- c("Reference", "Solution", "d15N", "s15N", "d18Ovsmow", "s18O", "Notes")
globModSol$legendName <- factor(paste(paste(as.numeric(globModSol$Reference)), "-", globModSol$Reference)) # for plotting
# separate Sowers since it is a box and not a point
globModSolBox <- subset(globModSol, subset = Notes == "box")
globModSol <- subset(globModSol, subset = Notes != "box")
levels(globModSol$Solution)[4] <- "Soil" # For plotting
###



###
# THIS WILL BE FIGURE 7 (the bottom-up top-down figure)
# Save at 7.20" x 6.83" and base_size = 12 and Arial for PLOS ONE
cols.Fig7 <- brewer_pal(pal = "Dark2")(5)
cols.Fig7 <- c("black", "grey", cols.Fig7[c(3, 1, 4, 2)], "grey", "black")
ggplot(mapping = aes(x = d15N, y = d18Ovsmow)) +
  geom_ribbon(data = globModSolBox, aes(ymin = 17, ymax = 26, shape = legendName, colour = Solution, fill = NA), show_guide = FALSE) +
  geom_path(data = subFMS.ell, size = 1, aes(colour = Category), show_guide = FALSE) +
  geom_point(data = Continental, size = 5, aes(x = mean(d15N), y = mean(d18Ovsmow), colour = Category), show_guide = FALSE) +
  geom_point(data = globModSol, size = 3, aes(colour = Solution, shape = legendName)) +
  geom_point(data = globModSol, size = 6, aes(colour = Solution), show_guide = FALSE) +
  geom_point(data = globModSol["10",], size = 3.5, colour = cols.Fig7[5], show_guide = FALSE) +
  geom_text(data = globModSol, size = 5, colour = "white", label = paste(as.numeric(globModSol$Reference)), show_guide = FALSE) +
  geom_text(data = data.frame(d15N = mean(globModSolBox$d15N), d18Ovsmow = mean(globModSolBox$d18Ovsmow)), size = 5, colour = "black", aes(label = "7")) +
  geom_point(data = Stratosphere, size = 3, shape = 3, aes(colour = Category), show_guide = FALSE) +
  geom_point(data = Troposphere, size = 3, shape = 8, aes(colour = Category), show_guide = FALSE) +
  #geom_point(data = StratoTropo, size = 3, aes(colour = Category)) +
  scale_x_continuous(limits = c(-30, 15)) +
  scale_y_continuous(limits = c(15, 65)) +
  scale_colour_manual(values = cols.Fig7) +
  scale_fill_manual(values = cols.Fig7) +
  scale_shape_manual(values = c(rep(NA,9))) +
  labs(x = expression(paste("\u03b4" ^15, "N-N" [2], "O (\u2030 air N" [2], ")")), 
       y = expression(paste("\u03b4" ^18, "O-N" [2], "O (\u2030 VSMOW)")), 
       linetype = "", shape = "       Studies", colour = "Sources", fill = "") +
  coord_fixed() +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(legend.position = c(0.35, 0.90), legend.background = element_rect(fill = NA),
        legend.margin = unit(0.4, "cm"), legend.key.size = unit(0.4, "cm"),
        legend.key = element_rect(colour = NA), legend.box.just = "top",
        legend.box = "horizontal",
        plot.margin = unit(c(0,0,0,0), "mm")) +
  guides(colour = guide_legend(order = 2, override.aes = list(shape = c(16,16,16,16,16,16,3,8))), 
         shape = guide_legend(order = 1))
###

# EOF 