library(glmmTMB)
library(DHARMa)
library(emmeans)
library(ggplot2)
library(vegan)

setwd("~/Desktop/Github/healthy_preterm_analysis/")

# import and format data
md <- read.csv("Metadata.csv")
sp <- read.csv("Species_relab.csv", header = TRUE, row.names =1)

# samples as rows
sp <- t(sp)

# order levels of covariates
md$EBM_Summary <- factor(md$EBM_Summary, levels = c("Never", "Before", "During", "After"))
md$whichprob <- factor(md$whichprob, levels = c("None", "I", "LB2"))

# calculate Shannon diversity for each sample
md$Shannon <- diversity(sp, index = "shannon", MARGIN = 1, base = exp(1))

# check data is normally distributed (Shannon div. normally is)
hist(md$Shannon)
shapiro.test(md$Shannon)

# fit model, using gaussian distribution
summary(Shannon <-
          glmmTMB(Shannon ~
                    whichprob +
                    Abx_7d +
                    Dayfullfeed72hrs +
                    Birthweight +
                    Deliverymode +
                    Sex +
                    Season +
                    Weight_z_score_diff +
                    GA_exact +
                    EBM_Summary +
                    BMF_Summary +
                    Formula_Summary +
                    DOL +
                    (1|SERVISPatientNo),
                  data = md))

# test model fit and residuals, might need to re-adjust model accordingly
testResiduals(Shannon)

# get global p-values
car::Anova(Shannon)

# post-hoc tests (change formula to get values for different signficant covariates from model)
x <- emmeans(Shannon, specs = pairwise ~ whichprob)
x

plot(x) +
  theme_bw() +
  xlab("Emmean") +
  theme(text = element_text(size = 16),
        legend.title = element_blank())
