# statistics for annabells project
setwd('E:/PhD/Undergrads/Annabell')

#look at distributions on masses, eggs, yolks
masses <- read.csv('egg.masses.byweek.csv') 
eggs <- read.csv('eggs.byweek.csv')
yolks <- read.csv('yolks.byweek.csv')
eggs.per <- read.csv('eggs.per.csv')

#comparisons are only informative prior to week 4 due to castration in the infected class
mass3 <- subset(masses, week <= 3)
egg3 <- subset(eggs, week <= 3)

mass3nox <- subset(mass3, Infection.Status != 1)
egg3nox <- subset(egg3, Infection.Status != 1)

library("ggpubr")
library('ggplot2')

egg.mass3 <- ggplot(mass3nox, aes(x = as.factor(week), y = egg.masses, fill = fact)) + 
  geom_violin(trim = T, scale = 'width', draw_quantiles = c(0.25, 0.5, 0.75), adjust = 1.5) +
  labs(x="Week", y = "Number of egg masses laid", hjust=10) +
  theme(text = element_text(size=50),
        axis.text.x = element_text(angle=90, hjust=1))
egg.mass3 + geom_dotplot(binaxis='y', stackdir='centerwhole', stackratio = 0.5, dotsize = 0.5, position=position_dodge(1)) +
  scale_fill_manual(values=c("#FC717F","#00BFC4")) + theme_classic() +
  theme(text = element_text(size=20), legend.title=element_blank())

egg.3 <- ggplot(egg3nox, aes(x = as.factor(week), y = eggs, fill = fact)) + 
  geom_violin(trim = T, draw_quantiles = c(0.25, 0.5, 0.75)) +
  labs(x="Week", y = "Number of eggs laid", hjust=10) +
  theme(text = element_text(size=50),
        axis.text.x = element_text(angle=90, hjust=1))#+ geom_point(position=position_jitterdodge())
egg.3 + geom_dotplot(binaxis='y', stackdir='centerwhole', stackratio = 0.5, dotsize = 0.5, position=position_dodge(1)) +
  scale_fill_manual(values=c("#FC717F","#00BFC4")) + theme_classic() +
  theme(text = element_text(size=20), legend.title=element_blank())

ggboxplot(mass3nox, x = 'week', y = 'egg.masses', color = 'fact')
ggboxplot(egg3nox, x = 'week', y = 'eggs', color = 'fact')


library(glmmTMB)
fit_hzip <- glmmTMB(egg.masses ~ fact1 * week + (1|SNAil.ID), data=mass3nox, ziformula=~., family=poisson(link="log"))
summary(fit_hzip)

fit_hzip2 <- glmmTMB(eggs ~ fact1 * week + (1|SNAil.ID), data=egg3nox, ziformula=~., family=poisson)
summary(fit_hzip2)


#check on egg masses
fit_zinbinom2 <- update(fit_hzip, family=nbinom2)
fit_zinbinom1 <- update(fit_hzip, family=nbinom1)
fit_hnbinom2 <- update(fit_zinbinom1, family=list(family="truncated_nbinom1",link="log"))

library(bbmle)
AICtab(fit_hzip, fit_zinbinom2, fit_zinbinom1, fit_hnbinom2)

#the best one
summary(fit_zinbinom1)

#check on eggs
fit_zinbinom2 <- update(fit_hzip2, family=nbinom2)
fit_zinbinom1 <- update(fit_hzip2, family=nbinom1)
fit_hnbinom2 <- update(fit_zinbinom1, family=list(family="truncated_nbinom1",link="log"))
AICtab(fit_hzip2, fit_zinbinom2, fit_zinbinom1, fit_hnbinom2)

#the best one
summary(fit_hnbinom2)


#data is continuous not count, so need to use a hurdle model or zero-inflated gamma model
#fit_zipoisson4 <- glmmTMB(eggs.egg.masses ~ fact1 + (1|SNAil.ID), data=eggs.per, ziformula=~1, family=poisson)
#summary(fit_zipoisson4)
noz <- subset(eggs.per, eggs.egg.masses > 0)
noz3 <- subset(noz, week <= 3)
noz3nox <- subset(noz3, Infection.Status != 1)

ggboxplot(noz3nox, x = 'week', y = 'eggs.egg.masses', color = 'fact')

library(lme4)
library(lattice)
library(car)
rep.noz <- lmer(eggs.egg.masses ~ fact1 * week + (1|SNAil.ID),
               data = noz3nox)
summary(rep.noz)
Betas  <- fixef(rep.noz)                  #Get the betas
SE     <-  sqrt(diag(vcov(rep.noz)))      #Get the SEs
pval   <- 2*pnorm(-abs(Betas  / SE)) #Z distribution
Output <- cbind(Betas,SE, pval)
print(Output, digits = 3) # P values for fixed effect
plot(rep.noz)
qqmath(rep.noz)
Anova(rep.noz, type="III")


#binomial glmms for hatching and growth
hatch <- read.csv('hatching.csv')
growth <- read.csv('growth.csv')

hatch.nox <- subset(hatch, Status != 1)
growth.nox <- subset(growth, Status != 1)

library(lme4)
hmm <- glmer(Hatched ~ Week + (1|ID) + (1|Jar) + (1|Status.1),
      data=hatch.nox,
      family=binomial)
library(car)
Anova(hmm)
summary(hmm)

hmm2 <- glmer(Hatched ~ Status.1 + (1|ID) + (1|Jar) + (1|Week),
             data=hatch.nox,
             family=binomial)
Anova(hmm2)
summary(hmm2)

library(emmeans)
multzip <- emmeans(hmm2, pairwise ~ Status.1)
test(multzip, adjust = 'mvt')


gmm <- glmer(Matured ~ Week + (1|ID) + (1|Jar) + (1|Status.1),
             data=growth.nox,
             family=binomial)
Anova(gmm)
summary(gmm)

gmm2 <- glmer(Matured ~ Status.1 + (1|ID) + (1|Jar) + (1|Week),
             data=growth.nox,
             family=binomial)
Anova(gmm2)
summary(gmm2)
multzip <- emmeans(gmm2, pairwise ~ Status.1)
test(multzip, adjust = 'mvt')



mist <- read.csv('logreg.errors.csv')
mmm <- glmer(error ~ group + (1|ID) + (1|week),
             data=mist,
             family=binomial)
summary(mmm)

library(emmeans)
multzip <- emmeans(mmm, list(pairwise ~ group), adjust = 'mvt')
multzip

