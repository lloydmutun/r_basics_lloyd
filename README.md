######### day 4 ############
# GENERAL LINEAR MIXED-EFFECTS MODELS #

# CLEAR ENVIRONMENT
rm(list = ls())

# Productivity of Mauritius parakeets (from Tollington et al 2015 Jnl Anim Ecol) # 

# SETTING WORKING DIRECTORY # 
setwd("~/R_LEARNING_Material/Participants/data and code")

# LOAD PACKAGES #
library(readr)
dat_ti <- read_csv("dat.ti.csv")
View(dat_ti)

dat.ti <- dat_ti
str(dat.ti)


# clutch size : What variables do you think might affect clutch size? Our 
# main hypothesis here is that in some way, an outbreak of disease negatively 
# affected some element of reproductive fitness.

clutch.dat <- dat.ti[, c("hatched", "fledged", "eggs", "year", "damage", 
    "damSH","lay1stsep", "supp","outbreak","female")]
str(clutch.dat)

# Now, we need to create a few extra variables.
# We also need to transform the ‘dam age’ variable as lots of research shows that many variables 
# of reproductive fitness are related to the squared age of females.

clutch.dat$yearf <- as.factor(clutch.dat$year)
clutch.dat$damage2 <- (clutch.dat$damage)^2
str(clutch.dat)

# We now need to remove the NAs.

clutch.dat <- na.omit(clutch.dat)
str(clutch.dat)

# Zuur code that allows us to visualise our data and any relationships that exist.

source("Zuur_VIFunction.R")

# First, we need to create another separate data frame to perform the pair plot analysis, 
# which we restrict to explanatory variables only as it is collinearity amongst covariates that 
# we are interested in isolating.
library(tidyverse)
z <- clutch.dat %>%
  select(eggs, damage2, damSH, lay1stsep, year)

# Correlation functions require numeric data so, for the purposes of this operation, 
# we have omitted all non-numeric variables ("female", "supp" and "outbreak").
z <- as.data.frame(z)
pairs(z,lower.panel=panel.smooth2,upper.panel=panel.cor,diag.panel=panel.hist)

# we can check if any single variable is unacceptably collinear with the others.

corvif(z)

# Now we are satisfied that our predictors are free of potentially confounding collinearity, 
# let’s check our response variable to see what it looks like.

hist(clutch.dat$eggs)

# How should our model look? What is the response variable? What are our explanatory 
# variables and which are our random variables?

eggs~(damage2+damSH+lay1stsep+supp)*outbreak+(1|female)+(1|yearf)

# if we were wishing to specify random intercepts for females and yearf, and also random slopes 
# for the covariate damage2 and yearf, our expression would look like this:

eggs~(damage2+damSH+lay1stsep+supp)*outbreak+(1|female)+(1 + damage2|yearf)

# Installing the lme4 package to fit our model.
install.packages("lme4")
library(lme4)

# Fitting the model with the lme4 package. We will start with a simple model 
# and then build up to the full model.

clutch <- lmer(eggs~(damage2 + damSH + lay1stsep + supp) * outbreak + 
                 (1|female)+(1|yearf), clutch.dat)
summary(clutch)

# let look into our recent model

confint(clutch)

# We are going to start the model averaging procedure with our model for 
# clutch size.

options(na.action=na.fail)

# standardising all the variables in our model
install.packages("arm")
library(arm)

stdzmodel <- standardize(clutch,standardize.y=FALSE)
summary(stdzmodel)

# We can now use the dredge function to perform model selection and model averaging.

install.packages("MuMIn")
library(MuMIn)
modeldredge <- dredge(stdzmodel)
modeldredge

# arranging our model by AIC
topmodel <- get.models(modeldredge,subset=delta<4)
topmodel

model.sel(topmodel)

# Our final model is therefore:

clutchfinal <- lmer(eggs~damage2+lay1stsep+(1|female)+(1|yearf),clutch.dat)
summary(clutchfinal)

# lets see what we can infer from our final model

confint(clutchfinal)
# create two more similar models to our final one. 
# In each one we remove one of the fixed effects.
cf1 <- lmer(eggs~damage2+(1|female)+(1|yearf),clutch.dat)
cf2 <- lmer(eggs~lay1stsep+(1|female)+(1|yearf),clutch.dat)
summary(cf1)
summary(cf2)

# We can now compare the AIC of these three models to see if either of the two
# fixed effects in our final model are not contributing to the model fit.

anova(clutchfinal,cf1)
anova(clutchfinal,cf2)

# So removing ‘lay date’ from our model significantly weakens it. Removing 
#‘damage’ does not significantly weaken the model, but it does increase the 
#AIC by 2.5, which is not a huge increase but it is still an increase. 
#o we would be justified in keeping both of these fixed effects in our 
#final model.
anova(clutchfinal,cf2)
# produce a Marginal R2 which denotes the variance explained by the fixed 
# effects only and a Conditional R2 which accounts for the variance explained 
# by the fixed and random components of the model.
library(MuMIn)
r.squaredGLMM(clutchfinal)

# To reduce loading times when executing diagnostic functions 
# for GLMMs it is advisable to simulate residuals from the model object.
install.packages("DHARMa")
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = clutchfinal, plot = F)


### We can plot that relationship pretty simply
library(ggplot2)
ggplot(clutch.dat, aes(x=lay1stsep, y=eggs)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=F) +
  labs(x=expression(lay~1^st~September), y="No. eggs") +
  theme_bw()

# We’ll need it to make a summary of female age categories and the values of 
#clutchsize that go with it.
library(tidyverse)
clutchmeanage <- dat.ti %>% 
  group_by(damage) %>% 
  summarise(N = length(eggs), 
            clutch=mean(eggs, na.rm=T), 
            sd=sd(eggs, na.rm=T), 
            se=sd/sqrt(N), 
            CI=(1.96*se))
clutchmeanage
ggplot(clutchmeanage, aes(x=damage, y=clutch)) +
  geom_point(shape=1) +
  geom_errorbar(aes(ymin=clutch-se, ymax=clutch+se), width=.2) +
  labs(x="Damage", y="No. eggs") +
  theme_bw()

ggplot(clutchmeanage, aes(x=damage, y=clutch)) + 
  geom_point(shape=1) + 
  geom_line(aes(x=damage, y = clutch), size=0.5, color = "black") + 
  geom_line(aes(x=damage,y=(clutch-CI)), size=0.5,color = "black", linetype = "dashed") + 
  geom_line(aes(x=damage,y=(clutch+CI)), size=0.5,color = "black", linetype = "dashed") + 
  labs(x="Female age", y="Mean clutch size") + 
  scale_x_continuous(breaks=seq(1,13,1)) + 
  theme_bw()















































