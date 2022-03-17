# This script can be used to fit multiple mixed effects logistic regression models (binomial)
# to see which clinical factors are associated with each preterm gut community type, and calculate aORs

library(glmmTMB)
library(DHARMa)
library(emmeans)

setwd("~/Desktop/Github/healthy_preterm_analysis")

data <- read.csv("Metadata.csv", header = TRUE, row.names = 1, sep = ',')

data$PGCT <- as.numeric(data$PGCT)

# reorder factor levels so that intuitive reference groups are used
data$EBM_Summary <- factor(data$EBM_Summary, levels = c("Never", "Before", "During", "After"))
data$Formula_Summary <- factor(data$Formula_Summary, levels = c("Never", "Before", "During", "After"))
data$BMF_Summary <- factor(data$BMF_Summary, levels = c("No", "Before", "During", "After"))
data$Deliverymode <- factor(data$Deliverymode, levels = c("Vaginal", "Caesarian"))
data$Sex <- factor(data$Sex, levels = c("Male", "Female"))
data$Abx_7d <- factor(data$Abx_7d, levels = c("No", "Yes"))
data$Weight_z_score_diff <- as.numeric(data$Weight_z_score_diff)
data$whichprob <- factor(data$whichprob, levels = c("None", "I", "LB2"))
data$Season <- factor(data$Season, levels = c("Winter", "Spring", "Summer", "Autumn"))


# sort response variable into 0 or 1
data.clst1 <- data
data.clst1$PGCT <- ifelse(data.clst1$PGCT == 1, data.clst1$PGCT <- 1, data.clst1$PGCT <- 0)

data.clst2 <- data
data.clst2$PGCT <- ifelse(data.clst2$PGCT == 2, data.clst2$PGCT <- 1, data.clst2$PGCT <- 0)

data.clst3 <- data
data.clst3$PGCT <- ifelse(data.clst3$PGCT == 3, data.clst3$PGCT <- 1, data.clst3$PGCT <- 0)

data.clst4 <- data
data.clst4$PGCT <- ifelse(data.clst4$PGCT == 4, data.clst4$PGCT <- 1, data.clst4$PGCT <- 0)

data.clst5 <- data
data.clst5$PGCT <- ifelse(data.clst5$PGCT == 5, data.clst5$PGCT <- 1, data.clst5$PGCT <- 0)

# fit model (binomial logistic regressions with mixed effects) for each cluster

data.input <- list(data.clst1, data.clst2, data.clst3, data.clst4, data.clst5)

clst.pval <- c()
ORs <- c()

for (i in 1:5){
    
  x <- data.frame(data.input[i])
  
  # fit each model, after checking residual plots etc. may need tweaking
  clst.model <-
          glmmTMB(PGCT ~ whichprob +
                    Birthweight +
                    Abx_7d +
                    EBM_Summary +
                    Formula_Summary +
                    BMF_Summary +
                    GA_exact +
                    Season +
                    Sex +
                    Deliverymode +
                    Weight_z_score_diff +
                    Dayfullfeed72hrs +
                    DOL +
                    (1|SERVISPatientNo),
                  data = x,
                  family = binomial)

summary(clst.model)

# get global p-values from the model (save in a list for all models)
clst.pval[[i]] <- car::Anova(clst.model)

# test the residuals from the model fit
DHARMa::testResiduals(clst.model)

# get OR and CI from the model (save in a list for all models)
ORs[[i]] <- exp(confint(clst.model, method = "wald", level = 0.95))

}
