### 26/11/2024
### Ewan Harney

#################################################################
# Phenotypic analysis for InterChromaTE
#
# Traits:  CTMax in the F2
#          Egg number
#          Pupae number
#          Adult number
#          Egg-to-pupae viability
#          Egg-to-adult viability
#          Age at pupation
#          Age at Eclosion
#
# Split into the following sections:
#     1: CTMax in the F2
#     2: F3 phenotypic results (not including age)
#     3: F3 phenotypic results (age)
#     4: F6 phenotypic results (not including age)
#     5: F6 phenotypic results (age)
#
# Note that F3 actually refers to offspring of F3 (so F4)
# and F6 actually refers to offspring of F6 (so F7)
#
#################################################################

rm(list=ls(all=TRUE))

library(ggplot2) # for all plots
library(cowplot) # for some composite plots
library(MASS) # for dropterm
library(lme4) # for mixed effects models
library(emmeans) # for posthoc tests
library(car) # for better anovas
library(effects)
library(interactions)

#################################################
#################################################
# 1: CTMax data
#################################################
#################################################

# set working dir
workingDir = "C:/../../SuppMat/interchromate_phenotypic";
setwd(workingDir); 
getwd();

# read in data
ctmax_dat<-read.table("ctmax_data_akaa_manz.txt", header = TRUE)

# model for CTMax
lmer_ctm <- lmer(Ctmax ~ Population +(1|Batch), data = ctmax_dat)
summary(lmer_ctm)
dropterm(lmer_ctm, test="Chisq")
Anova(lmer_ctm, type = "II")
emmeans(lmer_ctm, pairwise ~ Population, adjust = "sidak")

# plot of ctmax
ggplot(data=ctmax_dat, aes(x=Population, y = Ctmax, fill = Population)) +
  geom_boxplot(alpha=0.8) +
  ylim(37.8,42.1) +
  ylab("CT Max (°C)")+
  scale_fill_manual(values = c( "#6295cd", "#cb4f42")) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size=18,margin = margin(t = 20)),
        axis.title.y = element_text(size=18,margin = margin(r = 10)),
        legend.position = "none")


#################################################
#################################################
# 2: F3 phenotypic results (not including age)
#################################################
#################################################

# read in data
dataf3<-read.table ("Phenotypic_F3_eggs-pupa-adults-viab.txt", header=TRUE)
dataf3$Pop_Treat_Day<-factor(paste(dataf3$Pop_Treat, dataf3$Day, sep="_"))
dataf3$Pop_Treat_Day <- factor(dataf3$Pop_Treat_Day, levels = c("Akaa_Ctrl_1_2", "Akaa_T37_1_2", "Manz_Ctrl_1_2", "Manz_T37_1_2",
                                                                  "Akaa_Ctrl_3_4", "Akaa_T37_3_4", "Manz_Ctrl_3_4", "Manz_T37_3_4",
                                                                  "Akaa_Ctrl_5_6", "Akaa_T37_5_6", "Manz_Ctrl_5_6", "Manz_T37_5_6",
                                                                  "Akaa_Ctrl_13_14", "Akaa_T37_13_14", "Manz_Ctrl_13_14", "Manz_T37_13_14"))

# Filter data sets
# Remove observations with zero eggs
dataf3no0<-na.omit(subset(dataf3, eggs > 0 )) # 369 / 400 observations remaining
# Remove observations where the number of miscounted eggs was 4 or more (i.e. pupal counts exceeded egg counts by 4 or more)
dataf3lomc<-na.omit(subset(dataf3no0, EMiscount < 4 )) # 347 / 369 obs. We will use this as our main testing data set
nrow(dataf3lomc[dataf3lomc$EMiscount==0,]) # 294 samples with no miscounts, 53 with miscounts between 1 and 3

# split into first cohort and cohorts 2-4
dataf3_1<-subset(dataf3lomc, dataf3lomc$Day== "1_2") # 92
dataf3_24<-subset(dataf3lomc, dataf3lomc$Day!= "1_2") # 255

#####################
# 2.1: F3 eggs

F3eggs1<-glm(eggs ~ Treat * Pop, data = dataf3_1, family = poisson)
dropterm(F3eggs1,test="Chisq")
F3eggs1a<-glm(eggs ~ Treat + Pop, data = dataf3_1, family = poisson)
dropterm(F3eggs1a,test="Chisq")
Anova(F3eggs1a, type = "II")

F3eggs2<-glmer(eggs ~ Treat * Pop  + (1|Day) + (1|ID), data = dataf3_24, family = poisson, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
dropterm(F3eggs2,test="Chisq")
F3eggs2b<-glmer(eggs ~ Treat + Pop  + (1|Day) + (1|ID), data = dataf3_24, family = poisson, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
dropterm(F3eggs2b,test="Chisq")
Anova(F3eggs2b, type = "III")

# Effect Plots
# ggplot theme for all figures
pheno_theme <- theme(panel.grid.major.x = element_blank(),
                     axis.title.x=element_text(vjust = -1, size = 16),
                     axis.title.y=element_text(size = 16),
                     axis.text=element_text(size = 12),
                     legend.position="none",
                     strip.text.x = element_text(size=12, color="black",face="bold.italic"),
                     strip.background = element_rect(colour = "grey90", fill="grey90"),
                     panel.border = element_rect(fill = NA, color="grey90", linewidth =1))

# Effects plot set up
F3eggs1_eff<- allEffects(F3eggs1)
F3eggs1_plot <- F3eggs1_eff[[1]]
F3eggs1_plot_df <- as.data.frame(F3eggs1_plot)
F3eggs1_plot_df$Treat_Pop<-factor(paste(F3eggs1_plot_df$Treat, F3eggs1_plot_df$Pop, sep="_"))

# plot
c1eggs<-ggplot(F3eggs1_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Eggs") +
  ylim(1,37) +
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme

# Effects plot set up
F3eggs2_eff<- allEffects(F3eggs2)
F3eggs2_plot <- F3eggs2_eff[[1]]
F3eggs2_plot_df <- as.data.frame(F3eggs2_plot)
F3eggs2_plot_df$Treat_Pop<-factor(paste(F3eggs2_plot_df$Treat, F3eggs2_plot_df$Pop, sep="_"))

# plot
c2eggs<-ggplot(F3eggs2_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Eggs") +
  ylim(1,37) +
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme

#####################
# 2.2: F3 pupae

F3pupa1<-glm(pupa ~ Treat * Pop, data = dataf3_1, family = poisson)
dropterm(F3pupa1,test="Chisq")
F3pupa1a<-glm(pupa ~ Treat + Pop, data = dataf3_1, family = poisson)
dropterm(F3pupa1a,test="Chisq")
Anova(F3pupa1a, type = "II")

F3pupa2<-glmer(pupa ~ Treat * Pop + (1|Day) + (1|ID), data = dataf3_24_a, family = poisson)
dropterm(F3pupa2,test="Chisq")
F3pupa2a<-glmer(pupa ~ Treat + Pop + (1|Day) + (1|ID), data = dataf3_24_a, family = poisson)
dropterm(F3pupa2a,test="Chisq")
Anova(F3pupa2a, type = "II")

#### Effect Plots
F3pupa1_eff<- allEffects(F3pupa1)
F3pupa1_plot <- F3pupa1_eff[[1]]
F3pupa1_plot_df <- as.data.frame(F3pupa1_plot)
F3pupa1_plot_df$Treat_Pop<-factor(paste(F3pupa1_plot_df$Treat, F3pupa1_plot_df$Pop, sep="_"))

c1pupa<-ggplot(F3pupa1_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Pupae") +
  ylim(1,37) +
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme

F3pupa2_eff<- allEffects(F3pupa2)
F3pupa2_plot <- F3pupa2_eff[[1]]
F3pupa2_plot_df <- as.data.frame(F3pupa2_plot)
F3pupa2_plot_df$Treat_Pop<-factor(paste(F3pupa2_plot_df$Treat, F3pupa2_plot_df$Pop, sep="_"))

c2pupa<-ggplot(F3pupa2_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Pupae") +
  ylim(1,37) +
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme

#####################
# 2.3: F3 adults

F3imag1<-glm(imag ~ Treat * Pop, data = dataf3_1, family = poisson)
dropterm(F3imag1,test="Chisq")
F3imag1a<-glm(imag ~ Treat + Pop, data = dataf3_1, family = poisson)
dropterm(F3imag1a,test="Chisq")
Anova(F3imag1a, type = "II")

F3imag2<-glmer(imag ~ Treat * Pop + (1|Day) + (1|ID), data = dataf3_24, family = poisson, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
dropterm(F3imag2,test="Chisq")
F3imag2a<-glmer(imag ~ Treat + Pop + (1|Day) + (1|ID), data = dataf3_24, family = poisson, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
dropterm(F3imag2a,test="Chisq")
Anova(F3imag2a, type = "II")

# effect plots

F3imag1_eff<- allEffects(F3imag1)
F3imag1_plot <- F3imag1_eff[[1]]
F3imag1_plot_df <- as.data.frame(F3imag1_plot)
F3imag1_plot_df$Treat_Pop<-factor(paste(F3imag1_plot_df$Treat, F3imag1_plot_df$Pop, sep="_"))

c1adults<-ggplot(F3imag1_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Imagoes") +
  ylim(1,37) +
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme

F3imag2_eff<- allEffects(F3imag2)
F3imag2_plot <- F3imag2_eff[[1]]
F3imag2_plot_df <- as.data.frame(F3imag2_plot)
F3imag2_plot_df$Treat_Pop<-factor(paste(F3imag2_plot_df$Treat, F3imag2_plot_df$Pop, sep="_"))

c2adults<-ggplot(F3imag2_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Imagoes") +
  ylim(1,37) +
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme

plot_grid(c1eggs,c1pupa,c1adults, ncol = 3)
plot_grid(c2eggs,c2pupa,c2adults, ncol = 3)

#####################
# 2.4: F3 Egg-to-pupae viability

F3pvia1 <- glm(cbind(pupa,eggs) ~ Treat * Pop , data = dataf3_1, family = binomial)
dropterm(F3pvia1,test="Chisq")
F3pvia1a <- glm(cbind(pupa,eggs) ~ Treat + Pop , data = dataf3_1, family = binomial)
dropterm(F3pvia1a,test="Chisq")
Anova(F3pvia1a, "II")

F3pvia2 <- glmer(cbind(pupa,eggs) ~ Treat * Pop +(1|Day) + (1|ID), data = dataf3_24, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
dropterm(F3pvia2,test="Chisq")
F3pvia2a <- glmer(cbind(pupa,eggs) ~ Treat + Pop +(1|Day) + (1|ID), data = dataf3_24, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
dropterm(F3pvia2a,test="Chisq")
Anova(F3pvia2a, "II")

## effects plots
F3pvia1_eff<- allEffects(F3pvia1)
F3pvia1_plot <- F3pvia1_eff[[1]]
F3pvia1_plot_df <- as.data.frame(F3pvia1_plot)
F3pvia1_plot_df$Treat_Pop<-factor(paste(F3pvia1_plot_df$Treat, F3pvia1_plot_df$Pop, sep="_"))

viapup_c1<-ggplot(F3pvia1_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Egg-to-pupa viability") +
  ylim(0.12,0.55) +
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme

F3pvia2_eff<- allEffects(F3pvia2)
F3pvia2_plot <- F3pvia2_eff[[1]]
F3pvia2_plot_df <- as.data.frame(F3pvia2_plot)
F3pvia2_plot_df$Treat_Pop<-factor(paste(F3pvia2_plot_df$Treat, F3pvia2_plot_df$Pop, sep="_"))

viapup_c2<-ggplot(F3pvia2_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Egg-to-pupa viability") +
  ylim(0.12,0.55) +
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme

#####################
# 2.5: F3 Egg-to-adult viability

F3ivia1a1 <- glm(cbind(imag,eggs) ~ Treat * Pop, data = dataf3_1, family = binomial)
dropterm(F3ivia1a1,test="Chisq")
F3ivia1a2 <- glm(cbind(imag,eggs) ~ Treat + Pop, data = dataf3_1, family = binomial)
dropterm(F3ivia1a2,test="Chisq")
Anova(F3ivia1a2, "II")

F3ivia2b1 <- glmer(cbind(imag,eggs) ~ Treat * Pop + (1|Day) + (1|ID), data = dataf3_24, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
dropterm(F3ivia2b1,test="Chisq")
F3ivia2b2 <- glmer(cbind(imag,eggs) ~ Treat + Pop + (1|Day) + (1|ID), data = dataf3_24, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
dropterm(F3ivia2b2,test="Chisq")
Anova(F3ivia2b2, "II")

### Effects plots
F3ivia1_eff<- allEffects(F3ivia1)
F3ivia1_plot <- F3ivia1_eff[[1]]
F3ivia1_plot_df <- as.data.frame(F3ivia1_plot)
F3ivia1_plot_df$Treat_Pop<-factor(paste(F3ivia1_plot_df$Treat, F3ivia1_plot_df$Pop, sep="_"))

ggplot(F3ivia1_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Egg-to-adult viability") +
  ylim(0.00,0.52) +
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme

F3ivia2_eff<- allEffects(F3ivia2)
F3ivia2_plot <- F3ivia2_eff[[1]]
F3ivia2_plot_df <- as.data.frame(F3ivia2_plot)
F3ivia2_plot_df$Treat_Pop<-factor(paste(F3ivia2_plot_df$Treat, F3ivia2_plot_df$Pop, sep="_"))

ggplot(F3ivia2_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Egg-to-adult viability") +
  ylim(0.00,0.52) +
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme

#################################################
#################################################
# 3: F3 phenotypic results (age)
#################################################
#################################################

# read in data
f3pupae1_4<-read.table ("f3_all_age-to-pupation.txt", header=TRUE)
f3imagoes1_4<-read.table ("f3_all_age-to-eclosion.txt", header=TRUE)

f3pupae1_4$Pop_Treat_Cohort<-factor(paste(f3pupae1_4$Pop_Treat, f3pupae1_4$Cohort, sep="_"))
f3pupae1_4$Pop_Treat_Cohort <- factor(f3pupae1_4$Pop_Treat_Cohort, levels = c("Akaa_Ctrl_1", "Akaa_T37_1", "Manz_Ctrl_1", "Manz_T37_1",
                                                                   "Akaa_Ctrl_2", "Akaa_T37_2", "Manz_Ctrl_2", "Manz_T37_2",
                                                                   "Akaa_Ctrl_3", "Akaa_T37_3", "Manz_Ctrl_3", "Manz_T37_3",
                                                                   "Akaa_Ctrl_4", "Akaa_T37_4", "Manz_Ctrl_4", "Manz_T37_4"))
f3pupae2_4<-subset(f3pupae1_4, f3pupae1_4$Cohort != 1)
f3pupae1_4suc<-subset(f3pupae1_4, Pupate ==1)
f3pupae2_4suc<-subset(f3pupae2_4, Pupate ==1)


f3imagoes1_4$Pop_Treat_Cohort<-factor(paste(f3imagoes1_4$Pop_Treat, f3imagoes1_4$Cohort, sep="_"))
f3imagoes1_4$Pop_Treat_Cohort <- factor(f3imagoes1_4$Pop_Treat_Cohort, levels = c("Akaa_Ctrl_1", "Akaa_T37_1", "Manz_Ctrl_1", "Manz_T37_1",
                                                                              "Akaa_Ctrl_2", "Akaa_T37_2", "Manz_Ctrl_2", "Manz_T37_2",
                                                                              "Akaa_Ctrl_3", "Akaa_T37_3", "Manz_Ctrl_3", "Manz_T37_3",
                                                                              "Akaa_Ctrl_4", "Akaa_T37_4", "Manz_Ctrl_4", "Manz_T37_4"))
f3imagoes2_4<-subset(f3imagoes1_4, f3imagoes1_4$Cohort != 1)
f3imagoes1_4suc<-subset(f3imagoes1_4, Eclose ==1)
f3imagoes2_4suc<-subset(f3imagoes2_4, Eclose ==1)

#####################
# 3.1: F3 pupae Age

F3age1<-lm(Age ~ Treat * Pop + Density, data = data1.1)
dropterm(F3age1,test="F")
Anova(F3age1, type = "III")
# emmeans can be slow to run
emmeans(F3age1,  pairwise ~ Treat|Pop, adjust = "sidak", pbkrtest.limit = 4900)

F3age2<-lmer(Age ~ Treat * Pop + Density +(1|Cohort) + (1|ID), data = f3pupae2_4suc)
dropterm(F3age2,test="Chisq")
Anova(F3age2, type ="III")
emmeans(F3age2,  pairwise ~ Treat|Pop , adjust = "sidak", pbkrtest.limit = 4900)

## Effects plots
F3age1_eff<- allEffects(F3age1)
F3age1_plot <- F3age1_eff[[2]]
F3age1_plot_df <- as.data.frame(F3age1_plot)
F3age1_plot_df$Treat_Pop<-factor(paste(F3age1_plot_df$Treat, F3age1_plot_df$Pop, sep="_"))

timepup_c1<-ggplot(F3age1_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Time to pupation (days)") +
  ylim(10,15) +
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme

F3age2_eff<- allEffects(F3age2)
F3age2_plot <- F3age2_eff[[2]]
F3age2_plot_df <- as.data.frame(F3age2_plot)
F3age2_plot_df$Treat_Pop<-factor(paste(F3age2_plot_df$Treat, F3age2_plot_df$Pop, sep="_"))

timepup_c2<-ggplot(F3age2_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Time to pupation (days)") +
  ylim(10,15) +
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme

plot_grid(viapup_c1,timepup_c1,
          viapup_c2,timepup_c2)

#####################
# 3.2: F3 imago Age

F3image1<-lm(Age ~ Treat * Pop + Density_pup, data = data1.2)
dropterm(F3image1,test="Chisq")
Anova(F3image1, type = "III")
emmeans(F3image1,  pairwise ~ Treat|Pop , adjust = "sidak", pbkrtest.limit = 4900)

F3image2<-lmer(Age ~ Treat * Pop + Density_pup +(1|Cohort) +  (1|ID), data = f3imagoes2_4suc)
dropterm(F3image2,test="Chisq")
Anova(F3image2, type = "III")
emmeans(F3image2,  pairwise ~ Treat|Pop , adjust = "sidak", pbkrtest.limit = 4900)

### Effects plots
F3image1_eff<- allEffects(F3image1)
F3image1_plot <- F3image1_eff[[2]]
F3image1_plot_df <- as.data.frame(F3image1_plot)
F3image1_plot_df$Treat_Pop<-factor(paste(F3image1_plot_df$Treat, F3image1_plot_df$Pop, sep="_"))

ggplot(F3image1_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Time to eclosion") +
  ylim(15.95,21) + 
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme

F3image2_eff<- allEffects(F3image2)
F3image2_plot <- F3image2_eff[[2]]
F3image2_plot_df <- as.data.frame(F3image2_plot)
F3image2_plot_df$Treat_Pop<-factor(paste(F3image2_plot_df$Treat, F3image2_plot_df$Pop, sep="_"))

ggplot(F3image2_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Time to eclosion") +
  ylim(15.95,21) + 
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme

#################################################
#################################################
# 4: F6 phenotypic results (not including age)
#################################################
#################################################

# read in data
dataf6<-read.table ("Phenotypic_F6_eggs-pupa-adults-viab.txt", header=TRUE)
dataf6$Pop_Treat_Day<-factor(paste(dataf6$Pop_Treat, dataf6$Day, sep="_"))
dataf6$Pop_Treat_Day <- factor(dataf6$Pop_Treat_Day, levels = c("Akaa_Ctrl_1_2", "Akaa_T37_1_2", "Manz_Ctrl_1_2", "Manz_T37_1_2",
                                                                "Akaa_Ctrl_3_4", "Akaa_T37_3_4", "Manz_Ctrl_3_4", "Manz_T37_3_4",
                                                                "Akaa_Ctrl_5_6", "Akaa_T37_5_6", "Manz_Ctrl_5_6", "Manz_T37_5_6",
                                                                "Akaa_Ctrl_13_14", "Akaa_T37_13_14", "Manz_Ctrl_13_14", "Manz_T37_13_14"))

#####################
# 4.1: F6 eggs

# Filtered data sets
# Remove observations with zero eggs
dataf6_all<-na.omit(subset(dataf6, eggs > 0 )) # 264 / 288
# Remove observations where the number of miscounted eggs was 4 or more (i.e. pupal counts exceeded egg counts by 4 or more)
dataf6_all<-na.omit(subset(dataf6_all, EMiscount < 4 )) # 261 / 264
nrow(dataf6_all[dataf6_all$EMiscount==0,]) # 260 samples with no miscounts, 1 with miscount between 1 and 3

# No need to split into cohorts (unlike in the F3 there was no strong effect of heat shock on the first cohort)

F6eggs1<-glmer(eggs ~ Treat * Pop + (1|Day) + (1|ID), data = dataf6_all, family = poisson, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
dropterm(F6eggs1,test="Chisq")
F6eggs1a<-glmer(eggs ~ Treat + Pop + (1|Day) + (1|ID), data = dataf6_all, family = poisson, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
dropterm(F6eggs1a,test="Chisq")
Anova(F6eggs1a, type = "II")

## Effects Plots

F6eggs1_eff<- allEffects(F6eggs1)
F6eggs1_plot <- F6eggs1_eff[[1]]
F6eggs1_plot_df <- as.data.frame(F6eggs1_plot)
F6eggs1_plot_df$Treat_Pop<-factor(paste(F6eggs1_plot_df$Treat, F6eggs1_plot_df$Pop, sep="_"))

f7_eggs<-ggplot(F6eggs1_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Eggs") +
  ylim(2,43) +
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme


#####################
# 4.2: F6 pupae

F6pupa1<-glmer(pupa ~ Treat * Pop +(1|Day) + (1|ID), data = dataf6_all, family = poisson, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
dropterm(F6pupa1,test="Chisq")
F6pupa1a<-glmer(pupa ~ Treat + Pop +(1|Day) + (1|ID), data = dataf6_all, family = poisson, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
dropterm(F6pupa1a,test="Chisq")
Anova(F6pupa1a, type = "II")

## Effects plot
F6pupa1_eff<- allEffects(F6pupa1)
F6pupa1_plot <- F6pupa1_eff[[1]]
F6pupa1_plot_df <- as.data.frame(F6pupa1_plot)
F6pupa1_plot_df$Treat_Pop<-factor(paste(F6pupa1_plot_df$Treat, F6pupa1_plot_df$Pop, sep="_"))

f7_pupae<-ggplot(F6pupa1_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Pupae") +
  ylim(2,43) +
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme

#####################
# 4.3: F6 adults

F6imag1<-glmer(imag ~ Treat * Pop+(1|Day) + (1|ID), data = dataf6_all, family = poisson, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
dropterm(F6imag1,test="Chisq")
F6imag1a<-glmer(imag ~ Treat + Pop+(1|Day) + (1|ID), data = dataf6_all, family = poisson, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
dropterm(F6imag1a,test="Chisq")
Anova(F6imag1a, type ="II")

## Effects plot
F6imag1_eff<- allEffects(F6imag1)
F6imag1_plot <- F6imag1_eff[[1]]
F6imag1_plot_df <- as.data.frame(F6imag1_plot)
F6imag1_plot_df$Treat_Pop<-factor(paste(F6imag1_plot_df$Treat, F6imag1_plot_df$Pop, sep="_"))

f7_adults<-ggplot(F6imag1_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Imagoes") +
  ylim(2,43) +
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme

plot_grid(f7_eggs,f7_pupae,f7_adults,ncol = 3)

#####################
# 4.4: F6 Egg-to-pupae viability

F6pvia1 <- glmer(cbind(pupa,eggs) ~ Treat * Pop + (1|Day) +  (1|ID), data = dataf6_all, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
dropterm(F6pvia1,test="Chisq") 
F6pvia1a <- glmer(cbind(pupa,eggs) ~ Treat + Pop + (1|Day) +  (1|ID), data = dataf6_all, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
dropterm(F6pvia1a,test="Chisq")
Anova(F6pvia1a, type ="II")

# Effects plot
F6pvia1_eff<- allEffects(F6pvia1)
F6pvia1_plot <- F6pvia1_eff[[1]]
F6pvia1_plot_df <- as.data.frame(F6pvia1_plot)
F6pvia1_plot_df$Treat_Pop<-factor(paste(F6pvia1_plot_df$Treat, F6pvia1_plot_df$Pop, sep="_"))

f7_pupvia<-ggplot(F6pvia1_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Egg-to-pupa viability") +
  ylim(0.0,0.5)+
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme

#####################
# 4.5: F6 Egg-to-adult viability

F6ivia1 <- glmer(cbind(imag,eggs) ~ Treat * Pop  + (1|Day) + (1|ID), data = dataf6_all, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
dropterm(F6ivia1,test="Chisq")
F6ivia1a <- glmer(cbind(imag,eggs) ~ Treat + Pop   + (1|Day) + (1|ID), data = dataf6_all, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
dropterm(F6ivia1a,test="Chisq")
Anova(F6ivia1a, type="II")

## Effects plot
F6ivia1_eff<- allEffects(F6ivia1)
F6ivia1_plot <- F6ivia1_eff[[1]]
F6ivia1_plot_df <- as.data.frame(F6ivia1_plot)
F6ivia1_plot_df$Treat_Pop<-factor(paste(F6ivia1_plot_df$Treat, F6ivia1_plot_df$Pop, sep="_"))

ggplot(F6ivia1_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Egg-to-adult viability") +
  ylim(0.0,0.52)+
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme

#################################################
#################################################
# 5: F6 phenotypic results (age)
#################################################
#################################################

# read in data
f6pupae1_4<-read.table("f6_all_age-to-pupation.txt", header=TRUE)
f6imagoes1_4<-read.table("f6_all_age-to-eclosion.txt", header=TRUE)

f6pupae1_4$Pop_Treat_Cohort<-factor(paste(f6pupae1_4$Pop_Treat, f6pupae1_4$Cohort, sep="_"))
f6pupae1_4$Pop_Treat_Cohort <- factor(f6pupae1_4$Pop_Treat_Cohort, levels = c("Akaa_Ctrl_1", "Akaa_T37_1", "Manz_Ctrl_1", "Manz_T37_1",
                                                                              "Akaa_Ctrl_2", "Akaa_T37_2", "Manz_Ctrl_2", "Manz_T37_2",
                                                                              "Akaa_Ctrl_3", "Akaa_T37_3", "Manz_Ctrl_3", "Manz_T37_3",
                                                                              "Akaa_Ctrl_4", "Akaa_T37_4", "Manz_Ctrl_4", "Manz_T37_4"))
f6pupae1_4suc<-subset(f6pupae1_4, Pupate ==1)

f6imagoes1_4$Pop_Treat_Cohort<-factor(paste(f6imagoes1_4$Pop_Treat, f6imagoes1_4$Cohort, sep="_"))
f6imagoes1_4$Pop_Treat_Cohort <- factor(f6imagoes1_4$Pop_Treat_Cohort, levels = c("Akaa_Ctrl_1", "Akaa_T37_1", "Manz_Ctrl_1", "Manz_T37_1",
                                                                                  "Akaa_Ctrl_2", "Akaa_T37_2", "Manz_Ctrl_2", "Manz_T37_2",
                                                                                  "Akaa_Ctrl_3", "Akaa_T37_3", "Manz_Ctrl_3", "Manz_T37_3",
                                                                                  "Akaa_Ctrl_4", "Akaa_T37_4", "Manz_Ctrl_4", "Manz_T37_4"))
f6imagoes1_4suc<-subset(f6imagoes1_4, Eclose ==1)

#####################
# 5.1: F6 pupae Age

F6age1<-lmer(Age ~ Treat * Pop + Density + (1|Cohort) + (1|ID) , data = f6pupae1_4suc)
dropterm(F6age1,test="Chisq")
F6age2<-lmer(Age ~ Treat + Pop + Density + (1|Cohort) + (1|ID) , data = f6pupae1_4suc)
dropterm(F6age2,test="Chisq")
Anova(F6age1, type ="III")

# Effects plot
F6page1_eff<- allEffects(F6age1)
F6page1_plot <- F6page1_eff[[2]]
F6page1_plot_df <- as.data.frame(F6page1_plot)
F6page1_plot_df$Treat_Pop<-factor(paste(F6page1_plot_df$Treat, F6page1_plot_df$Pop, sep="_"))

f7_pupage<-ggplot(F6page1_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Time to pupation (days)") +
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme

plot_grid(f7_pupvia,f7_pupage)

#####################
# 5.2: F6 imago Age

F6image1<-lmer(Age ~ Treat * Pop + Density_pup + (1|Cohort) + (1|ID), data = f6imagoes1_4suc)
dropterm(F6image1,test="Chisq")
Anova(F6image1, type ="III")
# emmeans very slow to run
#emmeans(F6image1,  pairwise ~ Treat|Pop , adjust = "sidak", pbkrtest.limit = 3370)

# Effects plot
F6image1_eff<- allEffects(F6image1)
F6image1_plot <- F6image1_eff[[2]]
F6image1_plot_df <- as.data.frame(F6image1_plot)
F6image1_plot_df$Treat_Pop<-factor(paste(F6image1_plot_df$Treat, F6image1_plot_df$Pop, sep="_"))

ggplot(F6image1_plot_df, aes(x=Treat, y=fit)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9))+
  geom_point(aes(colour = Treat_Pop), size = 7, shape =15) +
  scale_colour_manual(values = c("deepskyblue1", "red1", "deepskyblue4", "red4")) +
  xlab("Treatment") + ylab("Time to ecolsion") +
  ylim(15.95,21)+
  facet_grid(. ~Pop) +
  theme_minimal()+
  pheno_theme
