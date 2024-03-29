library(sjPlot)
require(glmmTMB)
require(ggplot2)
require(lme4)
library(MASS)
library(reshape2) 
library(dplyr)
library(tidyr)
library(bbmle)
library(gridExtra)
library(viridis)
library(data.table) # for C val diff calculation

###### MAIN FIGURES ######

# Dataframe containing every plant within every population
master_boltbool <- read.csv(file = "CHDO_byPlant.csv", header = TRUE, sep = ",", na.strings = "NA")
# Dataframe containing every population
master_byPop <- read.csv(file = "CHDO_byPop.csv", header = TRUE, sep = ",", na.strings = "NA")
# Dataframe containing every estimate of genome size. Used for supplemental figures.
master_allCests <- read.csv(file = "CHDO_byGenome.csv", header = TRUE, sep = ",", na.strings = "NA")

# GLMER model using AHM as a predictor variable for bolting frequency
Bolting_Model <- glmer(Bolt_Bool ~ AHM + (1 | Population), data = master_boltbool, family=binomial)

master_nano <- with(master_boltbool, data.frame(Population,Bolt_Bool,AHM,AHM_ssp585))
masterfuture_nano <- with(master_boltbool, data.frame(Population,Bolt_Bool))
masterfuture_nano$AHM <- master_boltbool$AHM_ssp585

# Predict based on GLMER model contemporary and future bolting frequencies
real_predict_present <- predict(Bolting_Model, type= "response")
real_predict_present_fixef <- predict(Bolting_Model, re.form = NA,type= "response")
real_predict_future <- predict(Bolting_Model,newdata=masterfuture_nano, type= "response")
real_predict_future_fixef <- predict(Bolting_Model,newdata=masterfuture_nano, re.form = NA, type= "response")
master_nano$Predicted_F_Current <- real_predict_present
master_nano$Predicted_F_Current_fixef <- real_predict_present_fixef
master_nano$Predicted_F_ssp585 <- real_predict_future
master_nano$Predicted_F_ssp585_fixef <- real_predict_future_fixef

premelted_df <- with(master_nano, data.frame(Population,Predicted_F_Current,Predicted_F_ssp585,Predicted_F_Current_fixef,Predicted_F_ssp585_fixef,AHM,AHM_ssp585))
premelted_df_unique <- unique(premelted_df)
row.names(premelted_df_unique) <- NULL
premelted_df_unique$Frequency <- master_byPop$Frequency


### FIGURE 1 ###
# The observed frequency of bolting calculated for each of 95 Chaenactis douglasii 
# populations plotted against the predicted fixed effects values 
# (i.e., predicted frequency) from a generalized linear mixed model based on annual 
# heat moisture index as a climate predictor. A linear regression is represented 
# by the gray line. Pearson correlation coefficient (r) = 0.61, p < 0.0001.

cor.test(premelted_df_unique$Frequency,premelted_df_unique$Predicted_F_Current_fixef) # 0.6134
gg <- ggplot(data= premelted_df_unique, aes(Frequency, Predicted_F_Current_fixef)) +
  stat_smooth(method = "lm",se=FALSE,color="gray")+
  geom_point()+
  annotate("text",x=0.02,y=0.28,label=("r = 0.61")) +
  scale_x_continuous(breaks= seq(0, 1, 0.1)) +
  scale_y_continuous(breaks= seq(0, 1, 0.1)) +
  coord_cartesian(xlim=c(0,1),ylim=c(0,0.3)) +
  theme_bw() +
  xlab("Observed Frequency") + ylab("Predicted Frequency")
gg
ggsave("FreqObsVSPredicted_from_Bool_fixef+response.jpg", width = 9, height = 3)


### TABLE 1 ###
# General linear mixed-effects model output of presence and absence of bolting in 
# Chaenactis douglasii. Fixed effects were assigned to the annual heat moisture index (AHM).
# Confidence intervals (CI) and P values (p) are reported for fixed effects.
# Random effects are the assignment of individuals to seed source populations.
# σ2 = mean random effect variance, τ00 Population = variation between populations, 
# ICC = intraclass correlation coefficient (random effect variance / total variance).

tab_model(Bolting_Model)


### FIGURE 3 ###
# Predict from model bolting frequency based on contemporary and future predictors, as well as their fixed effects 

# Melt data to separate contemporary and future data into groups
melted_df <- pivot_longer(premelted_df_unique, cols=2:3,names_to="Data",values_to="Predicted")

# Rename values to appear more clearly on graph
melted_df <- melted_df %>%
  mutate(AHM_used = case_when(
    grepl("Current", Data) ~ AHM,
    grepl("ssp585", Data) ~ AHM_ssp585
  ))
melted_df <- melted_df %>%
  mutate(Data = case_when(
    grepl("Predicted_F_Current", Data) ~ "Current Predicted",
    grepl("Predicted_F_ssp585", Data) ~ "Mid-Century Predicted"
  ))

gg <- ggplot(melted_df, aes(x=AHM_used,y=Predicted)) +
  geom_line(aes(group = Population), alpha = 0.4) +
  geom_point(aes(color = Data)) +
  geom_line(aes(group = Data, color = Data), se = FALSE, stat='smooth', linetype = "dashed", alpha=1, size=0.7) +
  scale_x_continuous(breaks= seq(0, 100, 10)) +
  scale_y_continuous(breaks= seq(0, 1, 0.1)) +
  xlab('Annual Heat Moisture Index of Population')+ylab("Predicted Frequency")
gg
ggsave("AHM_VS_Current&Future_SSP585_Freq+smooth.jpg", width = 8, height = 8)


###### SUPPLEMENTAL FIGURES ######

### APPENDIX S1 ###
# A histogram of sample counts and 1C values (pg) for Chaenactis douglasii.
# Colors indicate putative cytotypes based on the genome size. 


# Apply rough estimates of ploidy
master_allCests <- master_allCests %>%
  mutate(Estimated_Cytotype = case_when(
    X1C.value < 3.4 ~ "2x",
    X1C.value >= 3.4 & X1C.value < 4.25 ~ "3x",
    X1C.value >= 4.25 & X1C.value < 6.3 ~ "4x",
    X1C.value >= 6.3 ~ "6x"
  ))

model_plot <- ggplot(master_allCests, aes(x=X1C.value,fill=Estimated_Cytotype)) +
  geom_histogram(binwidth=0.03) +
  scale_fill_manual(values=c("#E69F00","#009E73","#0072B2","#CC79A7")) +
  scale_x_continuous(breaks= seq(0, 9, 0.5)) +
  scale_y_continuous(breaks= seq(0, 12, 1)) +
  xlab('1C Value')+ylab("Count") +
  theme(legend.position = c(0.86, 0.78)) +
  guides(fill=guide_legend(title="Estimated Cytotype"))
model_plot
ggsave("1C_Histogram_withCytotype_USEDPOPS.jpg", width = 7, height = 4)


### APPENDIX S3 ###
# Heat map matrix depicting the Pearson's correlation coefficient between 
# significant predictor variables (correlated with bolting frequency at 
# r > 0.4 or r < -0.4). A white dot represents a co-linear relationship 
# (r > 0.7 or r < -0.7). Climate acronyms can be found in appendix S2.

# Restrict to relevant climate variables and melt correlation matrix
datapop2 = subset(master_byPop,select=c(Frequency,AHM_var,DD18_at,SHM_var,AHM,DD18,DD18_sm,SHM,MAT,DD1040,DD5,MWMT,CMD,DD_18_wt,DD_18,DD_18_at))
cor2 <- cor(datapop2)
melted_cor_df2 <- melt(cor2)
melted_cor_df2 <- melted_cor_df2 %>% 
  mutate(dot_bool = case_when(
    value <= -0.7 ~ TRUE,
    (value > -0.7) & (value < 0.7) ~ FALSE,
    value >= 0.7 ~ TRUE
  ))
ggplot(melted_cor_df2, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() + scale_fill_viridis(discrete = FALSE)+
  geom_tile() +
  geom_point(aes(size=ifelse(dot_bool, TRUE, NA)),show.legend=FALSE,color="white") +
  theme(axis.text.x = element_text(angle=90, hjust=1)) + 
  labs(x = "",
       y = "")
ggsave("Climate_Heatmap.jpg", width = 9.4, height = 8)


### APPENDIX S4 ###
# Predicted bolting frequency difference between current and mid-century for
# Chaenactis douglasii populations and the percent change, derived from 
# contemporary and mid-century annual heat moisture index values at each
# population location as a fixed effect predictor.

# Calculate change in frequency from contemporary to future predicted data
premelted_df_unique$Frequency_Change <- premelted_df_unique$Predicted_F_ssp585_fixef - premelted_df_unique$Predicted_F_Current_fixef
premelted_df_unique$Percent_Frequency_Change <- premelted_df_unique$Frequency_Change / premelted_df_unique$Predicted_F_Current_fixef
suppfigure_df <- with(premelted_df_unique, data.frame(Population, Frequency_Change, Percent_Frequency_Change))
grid.table(suppfigure_df)


### APPENDIX S5 ###
# Linear mixed-effects model results of genome size variation in Chaenactis douglasii.
# Fixed effects, climate predictors are variation in the number of frost-free days
# (NFFD var) and summer (June, July, and August) NFFD (NFFD sm). Confidence intervals
# (CI) and p-values (p) are reported for fixed effects. Random effects are the
#assignment of individuals to seed source populations.  σ2 = mean random effect
# variance, τ00 Population = variation between populations, 
# ICC = intraclass correlation coefficient (random effect variance / total variance).

# Dataframe containing every sample of the used populations whose genome size was measured
master_allCests <- read.csv(file = "Frankenstein1CwithClimate.csv", header = TRUE, sep = ",", na.strings = "NA")
Cmodel <- lmer(X1C.value ~ NFFD_var + NFFD_sm + (1 | Population), data = master_allCests)
tab_model(Cmodel)


### APPENDIX S6 ###
# General linear mixed-effects model results of presence and absence of bolting
# in Chaenactis douglasii, using yearly variation in Annual Heat Moisture from
# 1991 to 2020 as a predictor variable.

Bolting_Model_var <- glmer(Bolt_Bool ~ AHM_var + (1 | Population), data = master_boltbool, family=binomial)
tab_model(Bolting_Model_var)

