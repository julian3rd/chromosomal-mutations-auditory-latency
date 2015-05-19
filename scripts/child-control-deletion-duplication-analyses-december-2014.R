# clear all data
rm(list=ls())


# data input --------------------------------------------------------------

#set path and file

path <- '~/GitHub/chromosomal-mutations-auditory-latency/data/'

file  <- 
  file.path(path, 'simons-child-complete-data-long-format.csv')

#read in data
all.child.data <- read.csv(file, header = T)

# total recruited children
child.recruited <- subset(all.child.data, 
                     Chromosome == 'control' & Case == 'control' |
                       Chromosome == '16p' & Case == 'deletion' |
                       Chromosome == '16p' & Case == 'duplication')

child.recruited <- droplevels(child.recruited)

# remove artifacts, exlcusions and have only 16p
child.data <- subset(all.child.data, Exclusions == 'no')

child.data <- subset(child.data, 
                     Chromosome == 'control' & Case == 'control' |
                       Chromosome == '16p' & Case == 'deletion' |
                       Chromosome == '16p' & Case == 'duplication')

# only using validated data, so remove any missing diagnoses
child.data <- subset(child.data, 
                     ASD == 'FALSE' | ASD == 'TRUE' | ASD == 'control')

child.data <- droplevels(child.data)

# setting yearly age and 12yo age brackets and cutoff ---------------------

child.data$M100compCase <- 
  as.numeric(complete.cases(child.data$M100LatCorr))

child.data$M50compCase <- 
  as.numeric(complete.cases(child.data$M50LatCorr))

child.data$cutAge <-
  as.factor(ifelse(child.data$Age < 12,'under-12','12-and-over'))

min.age <- floor(min(child.data$Age_Calc))
max.age <- ceiling(max(child.data$Age_Calc))

child.data$breakAge <-
  cut(child.data$Age, breaks = seq(min.age, max.age, 1), include.lowest = T)


# case data subsets (control, deletion, duplication) ----------------------

child.control.data <- subset(child.data, Case == 'control')
child.deletion.data <- subset(child.data, Case == 'deletion')
child.duplication.data <- subset(child.data, Case == 'duplication')

child.control.data <- droplevels(child.control.data)
child.deletion.data <- droplevels(child.deletion.data)
child.duplication.data <- droplevels(child.duplication.data)

# data in wide format -----------------------------------------------------

#wide format
library(reshape2)
child.wide.analysis <-
  dcast(data = child.data, formula = Subject + Site + Age_Calc + Gender + 
          Handedness + ASD + NVIQ + VIQ + CELF.4 + SRS_parent + CTOPP + 
          Case + cutAge + breakAge ~ Hem + Cond,
        value.var = "M100LatCorr")

child.recruited.wide <-
  dcast(data = child.recruited, formula = Subject + Site + Age_Calc + Gender + 
          Handedness + ASD + NVIQ + VIQ + CELF.4 + SRS_parent + CTOPP + 
          Case + Exclusions ~ Hem + Cond,
        value.var = "M100LatCorr")


# color-blind friendly palettes -------------------------------------------

cbbPalette <- 
  c("#000000", "#E69F00", "#56B4E9", 
    "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbPalette <-
  c("#999999", "#E69F00", "#56B4E9", 
    "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



# Site LMM ----------------------------------------------------------------
#models
library(lme4)
library(ggplot2)

site.lmm <- lmer(M100LatCorr ~ Site + Age_Calc +(Cond + Hem||Subject), 
       data = child.data, na.action = na.exclude, REML = F)


# LMMs --------------------------------------------------------------------

child.fullmodel <- 
  lmer(M100LatCorr ~ Hem * Cond * Case + Age_Calc + Site + (Cond + Hem|Subject),
       data = child.data, na.action = na.exclude, REML = F)

child.addmodel <- 
  lmer(M100LatCorr ~ Hem + Cond + Case + Age_Calc + Site + (Cond + Hem|Subject),
       data = child.data, na.action = na.exclude, REML = F)

child.addmodel.lh <- 
  lmer(M100LatCorr ~ Cond + Case + Age_Calc + Site + (1|Subject),
       data = child.data, subset = Hem == '1-LH',
       na.action = na.exclude, REML = F)

child.addmodel.rh <- 
  lmer(M100LatCorr ~ Cond + Case + Age_Calc + Site + (1|Subject),
       data = child.data, subset = Hem == '2-RH',
       na.action = na.exclude, REML = F)

child.addmodel.case.age.interact <- 
  lmer(M100LatCorr ~ Hem + Cond + Case + Age_Calc + 
         Site + Case:Age +(Cond + Hem|Subject),
       data = child.data, na.action = na.exclude, REML = F)

child.addmodel.case.cond.interact <- 
  lmer(M100LatCorr ~ Hem + Cond + Case + Age_Calc + Site + 
         Cond:Case + (Cond + Hem|Subject),
       data = child.data, na.action = na.exclude, REML = F)

child.addmodel.case.cond.case.age.interact <- 
  lmer(M100LatCorr ~ Hem + Cond + Case + Age_Calc + Site + 
         Cond:Case + Case:Age + (Cond + Hem|Subject),
       data = child.data, na.action = na.exclude, REML = F)

# appending fitted and residuals ------------------------------------------

# appending fitted values and residuals to main effects models
child.addmodel.fortify <- fortify(child.addmodel)

child.addmodel.case.cond.interact.fortify <- 
  fortify(child.addmodel.case.cond.interact)

#replace NAs with population-level values
child.addmodel.fortify$Predicted <- 
  predict(child.addmodel, newdata = child.addmodel.fortify, 
          na.action = na.pass, allow.new.levels = T)

child.addmodel.case.cond.interact.fortify$Predicted <- 
  predict(child.addmodel.case.cond.interact, 
          newdata = child.addmodel.case.cond.interact.fortify, 
          na.action = na.pass, allow.new.levels = T)

# pscyh score covariate and gender LMMs -----------------------------------

# main effects models with psych score covariates and interactions
child.addmodel.nviq <- 
  lmer(M100LatCorr ~ Hem + Cond + Case + Age_Calc + Site + NVIQ + 
         (Cond + Hem|Subject), data = child.data, 
       na.action = na.exclude, REML = F)

child.addmodel.viq <- 
  lmer(M100LatCorr ~ Hem + Cond + Case + Age_Calc + Site + VIQ + 
         (Cond + Hem|Subject), data = child.data,
       na.action = na.exclude, REML = F)

child.addmodel.srs <- 
  lmer(M100LatCorr ~ Hem + Cond + Case + Age_Calc + Site + SRS_parent +  
         (Cond + Hem|Subject), data = child.data,
       na.action = na.exclude, REML = F)

child.addmodel.icv <- 
  lmer(M100LatCorr ~ Hem + Cond + Case + Age_Calc + Site + cmICV + 
         (Cond + Hem|Subject), data = child.data, na.action = na.exclude, REML = F)

child.addmodel.nviq.viq.srs <- 
  lmer(M100LatCorr ~ Hem + Cond + Case + Age_Calc + NVIQ + VIQ + 
         SRS_parent + Site +  (Cond + Hem|Subject),
       data = child.data, na.action = na.exclude, REML = F)

#gender and handedness
gender.lmm <-
  lmer(M100LatCorr ~ Gender + Age_Calc + Site + (Cond + Hem|Subject),
       data = child.data, na.action = na.exclude, REML = F)

hand.lmm <-
  lmer(M100LatCorr ~ Handedness + Age_Calc + Site + (Cond + Hem|Subject),
       data = child.data, na.action = na.exclude, REML = F)


# p-values and significance -----------------------------------------------

# p-values: use multcomp and car packages

# Wald type II Chi-square tests: significance of factors
library(car)
child.fullmodel.anova <- Anova(child.fullmodel)

child.addmodel.anova <- Anova(child.addmodel)

child.addmodel.case.cond.interact.anova <- 
  Anova(child.addmodel.case.cond.interact)

library(multcomp)
# multiple comparisons and effect sizes via Tukey HSD
child.fullmodel.hem.tukey <- 
  summary(glht(child.fullmodel,
               linfct = mcp(Hem = 'Tukey', covariate_average = T,
                            interaction_average = T)))

child.fullmodel.cond.tukey <- 
  summary(glht(child.fullmodel,
               linfct = mcp(Cond = 'Tukey', covariate_average = T,
                            interaction_average = T)))

child.fullmodel.case.tukey <- 
  summary(glht(child.fullmodel,
               linfct = mcp(Case = 'Tukey', covariate_average = T,
                            interaction_average = T)))

child.addmodel.hem.tukey <- 
  summary(glht(child.addmodel, linfct = mcp(Hem = 'Tukey', covariate_average = T)))

child.addmodel.cond.tukey <- 
  summary(glht(child.addmodel, linfct = mcp(Cond = 'Tukey', covariate_average = T)))

child.addmodel.case.tukey <- 
  summary(glht(child.addmodel, linfct = mcp(Case = 'Tukey', covariate_average = T)))

child.addmodel.case.cond.interact.hem.tukey <- 
  summary(glht(child.addmodel.case.cond.interact,
               linfct = mcp(Hem = 'Tukey',
                            covariate_average = T, interaction_average = T)))

child.addmodel.case.cond.interact.cond.tukey <- 
  summary(glht(child.addmodel.case.cond.interact,
               linfct = mcp(Cond = 'Tukey', 
                            covariate_average = T, interaction_average = T)))

child.addmodel.case.cond.interact.case.tukey <- 
  summary(glht(child.addmodel.case.cond.interact,
               linfct = mcp(Case = 'Tukey',
                            covariate_average = T, interaction_average = T)))



# success rate ------------------------------------------------------------

library(plyr)

# success rate by group and proportion tests
group.complete <-  ddply(child.data, .(Case), summarize,
        compCases = sum(complete.cases(M100LatCorr)),
        n = length(M100LatCorr),
        failCases = n - compCases)


# instead of two-sample proportion test, use analysis of deviance

# Q: should age and subject effects be considered?
group.complete.glm <- glm(M100compCase ~ Case, binomial, data = child.data)

group.complete.lmm <- 
  glmer(M100compCase ~ Case + Age_Calc + (1|Subject), 
        binomial, data = child.data)

# Tukey HSD for complete cases; 
# effect sizes must be back-transformed using plogis()
group.complete.tukey <- 
  summary(glht(group.complete.glm, linfct = mcp(Case = 'Tukey')))

# successs counts: Hem, Cond, Hem x Cond (and all three by group)
# added fail and total column to perform GLMs on proportions
hem.comp <- ddply(child.data, .(Hem), summarize,
        compCases = sum(complete.cases(M100LatCorr)),
        failCases = length(M100LatCorr) - sum(complete.cases(M100LatCorr)),
        totalCases = length(M100LatCorr))

case.comp <-
  ddply(child.data, .(Case), summarize,
        compCases = sum(complete.cases(M100LatCorr)),
        failCases = length(M100LatCorr) - sum(complete.cases(M100LatCorr)),
        totalCases = length(M100LatCorr))

cond.comp <- ddply(child.data, .(Cond), summarize,
        compCases = sum(complete.cases(M100LatCorr)),
        failCases = length(M100LatCorr) - sum(complete.cases(M100LatCorr)),
        totalCases = length(M100LatCorr))

hem.cond.comp <- ddply(child.data, .(Hem, Cond),summarize,
        compCases = sum(complete.cases(M100LatCorr)),
        failCases = length(M100LatCorr) - sum(complete.cases(M100LatCorr)),
        totalCases = length(M100LatCorr))


hem.case.comp <- ddply(child.data, .(Hem, Case), summarize,
        compCases = sum(complete.cases(M100LatCorr)),
        failCases = length(M100LatCorr) - sum(complete.cases(M100LatCorr)),
        totalCases = length(M100LatCorr))


cond.case.comp <- ddply(child.data, .(Cond, Case), summarize,
        compCases = sum(complete.cases(M100LatCorr)),
        failCases = length(M100LatCorr) - sum(complete.cases(M100LatCorr)),
        totalCases = length(M100LatCorr))


hem.cond.case.comp <- ddply(child.data, .(Hem, Case, Cond), summarize,
        compCases = sum(complete.cases(M100LatCorr)),
        failCases = length(M100LatCorr) - sum(complete.cases(M100LatCorr)),
        totalCases = length(M100LatCorr))


subject.comp <- ddply(child.data, .(Subject), summarize,
        compCases = sum(complete.cases(M100LatCorr)),
        failCases = length(M100LatCorr) - sum(complete.cases(M100LatCorr)),
        totalCases = length(M100LatCorr))

subject.case.comp <- ddply(child.data, .(Subject, Case), summarize,
        compCases = sum(complete.cases(M100LatCorr)),
        failCases = length(M100LatCorr) - sum(complete.cases(M100LatCorr)),
        totalCases = length(M100LatCorr))

# plots of complete observations bu subject and subject x case
# count of total number of subjects with given number of complete observations
subject.complete.plot <- ggplot(subject.comp, aes(x = factor(compCases))) + 
  geom_bar() + theme_bw() +
  labs(title = 'Complete observations count - overall',
        x = 'Complete Observations', y = 'Number of Subjects')+
  coord_cartesian(ylim = c(0, 45)) + scale_y_continuous(breaks = seq(5, 45, 5))

subject.case.complete.plot <- 
  ggplot(subject.case.comp, aes(x = factor(compCases))) +
  geom_bar(aes(fill = Case), position = 'dodge') + theme_bw() +
  labs(title = 'Complete observations count by Case',
       x = 'Complete Observations', y = 'Number of Subjects')+
  coord_cartesian(ylim = c(0, 25))+
  scale_y_continuous(breaks = seq(0, 25, 2)) +
  theme(strip.text = element_text(face = 'bold', size = rel(1.5)),
        strip.background=element_rect(fill = 'white', 
                                      colour = 'white', size = 1),
        legend.position = 'bottom',
        legend.text = element_text(face = 'bold', size = 16),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15, face = 'bold'),
        axis.text.y = element_text(size = 15,face = 'bold'),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 16, face = 'bold'),
        axis.title.y = element_text(size = 16, face = 'bold'),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = 20, face = 'bold'))+
  scale_fill_manual(values = cbPalette)




# least-squares means: full model -----------------------------------------

child.fullmodel.lsmeans.table <-
  doBy::LSmeans.lmerMod(child.fullmodel, effect = c('Hem', 'Cond', 'Case'))

child.fullmodel.lsmeans <- 
  cbind(child.fullmodel.lsmeans.table$grid, child.fullmodel.lsmeans.table$coef)

child.fullmodel.lsmeans$Hem <- as.factor(child.fullmodel.lsmeans$Hem)
child.fullmodel.lsmeans$Cond <- as.factor(child.fullmodel.lsmeans$Cond)
child.fullmodel.lsmeans$Case <- as.factor(child.fullmodel.lsmeans$Case)


dodge <- position_dodge(width = 0.9)

child.fullmodel.lsmeans$Hem <- 
  revalue(child.fullmodel.lsmeans$Hem, c('1-LH' = 'LH', '2-RH' = 'RH'))

child.fullmodel.lsmeans$Cond <- 
  revalue(child.fullmodel.lsmeans$Cond,
          c('1-200' = '200', '2-300' = '300', '3-500' = '500', '4-1000' = '1000'))


marginal.means.plot <-
  ggplot(child.fullmodel.lsmeans,  aes(x = Cond, y = estimate, fill = Case, group = Case)) +
  geom_bar(stat = 'identity', position = 'dodge', size = 1.3)+
  geom_errorbar(size = 1.2, colour = "black", 
                aes(ymin = estimate - se, ymax = estimate + se),
                position = dodge, width = 0.3)+
  facet_grid(. ~ Hem) + theme_bw() +
  labs(x = 'Stimulus condition (Hz)', y = 'Latency (ms)',
       title = 'M100 least squares means ',
       fill = 'Case', colour = 'Case', face = 'bold', size = 20)+
  coord_cartesian(ylim = c(90, 170)) +
  theme(strip.text = element_text(face = 'bold', size = rel(1.5)),
        strip.background=element_rect(fill = 'white', 
                                      colour = 'white', size = 1),
        legend.position = 'bottom',
        legend.text = element_text(face = 'bold', size = 16),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15, face = 'bold'),
        axis.text.y = element_text(size = 15,face = 'bold'),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 16, face = 'bold'),
        axis.title.y = element_text(size = 16, face = 'bold'),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = 20, face = 'bold'))+
  scale_fill_manual(values = cbPalette)



# least-squares means: condition ------------------------------------------

# means across Hemisphere

child.cond.lsmeans.table <-
  doBy::LSmeans.lmerMod(child.fullmodel, effect = c('Cond', 'Case'))

child.cond.lsmeans <- 
  cbind(child.cond.lsmeans.table$grid, child.cond.lsmeans.table$coef)

child.cond.lsmeans$Cond <- as.factor(child.cond.lsmeans$Cond)
child.cond.lsmeans$Case <- as.factor(child.cond.lsmeans$Case)

child.cond.lsmeans$Cond <- 
  revalue(child.cond.lsmeans$Cond, c('1-200' = '200', '2-300' = '300', '3-500' = '500', '4-1000' = '1000'))

cond.means.plot <-
  ggplot(child.cond.lsmeans, 
         aes(x = Cond, y = estimate, fill = Case, group = Case))+
  geom_bar(stat = 'identity', position = 'dodge', size = 1.3)+
  geom_errorbar(size = 1.2, colour = "black", 
                aes(ymin = estimate - se, ymax = estimate + se),
                position = dodge, width = 0.3)+
  theme_bw() +
  labs(x = 'Stimulus condition (Hz)', y = 'Latency (ms)',
       title = 'Condition means',
       fill = 'Case', colour = 'Case', face = 'bold', size = 20)+
  coord_cartesian(ylim = c(90, 170)) +
  theme(strip.text = element_text(face = 'bold', size = rel(1.5)),
        strip.background=element_rect(fill = 'white', 
                                      colour = 'white', size = 1),
        legend.position = 'none',
        legend.text = element_text(face = 'bold', size = 16),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15, face = 'bold'),
        axis.text.y = element_text(size = 15,face = 'bold'),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 16, face = 'bold'),
        axis.title.y = element_text(size = 16, face = 'bold'),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = 20, face = 'bold'))+
  scale_fill_manual(values = cbPalette) +
  annotate('text', x = 4, y = 160, label = 'b', size = 8)



# least-squares means: hemisphere -----------------------------------------


# means across condition

child.hem.lsmeans.table <-
  doBy::LSmeans.lmerMod(child.fullmodel, effect = c('Hem', 'Case'))

child.hem.lsmeans <- 
  cbind(child.hem.lsmeans.table$grid, child.hem.lsmeans.table$coef)

child.hem.lsmeans$Hem <- as.factor(child.hem.lsmeans$Hem)
child.hem.lsmeans$Case <- as.factor(child.hem.lsmeans$Case)

child.hem.lsmeans$Hem <- 
  revalue(child.hem.lsmeans$Hem, c('1-LH' = 'LH', '2-RH' = 'RH'))

hem.means.plot <-
  ggplot(child.hem.lsmeans, aes(x = Hem, y = estimate, fill = Case, group = Case))+
  geom_bar(stat = 'identity', position = 'dodge', size = 1.3)+
  geom_errorbar(size = 1.2, colour = "black", 
                aes(ymin = estimate - se, ymax = estimate + se),
                position = dodge, width = 0.3)+
  theme_bw() +
  labs(x = 'Hemisphere', y = 'Latency (ms)',
       title = 'Hemisphere means',
       fill = 'Case', colour = 'Case', face = 'bold', size = 20)+
  coord_cartesian(ylim = c(90, 170)) +
  theme(strip.text = element_text(face = 'bold', size = rel(1.5)),
        strip.background=element_rect(fill = 'white', 
                                      colour = 'white', size = 1),
        legend.position = 'none',
        legend.text = element_text(face = 'bold', size = 16),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15, face = 'bold'),
        axis.text.y = element_text(size = 15,face = 'bold'),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 16, face = 'bold'),
        axis.title.y = element_text(size = 16, face = 'bold'),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = 20, face = 'bold'))+
  scale_fill_manual(values = cbPalette) +
  annotate('text', x = 2, y = 160, label = 'c', size = 8)


# least-squares means: case -----------------------------------------------

child.case.lsmeans.table <- 
  doBy::LSmeans.lmerMod(child.addmodel.case.cond.interact, effect = c('Case'))

child.case.lsmeans <- 
  cbind(child.case.lsmeans.table$grid, child.case.lsmeans.table$coef)

child.case.lsmeans$Case <- as.factor(child.case.lsmeans$Case)

case.means.plot <-
  ggplot(child.case.lsmeans, 
         aes(x = Case, y = estimate, fill = Case, group = Case)) +
  geom_bar(stat = 'identity', position = 'dodge', size = 1.3)+
  geom_errorbar(size = 1.2, colour = "black", 
                aes(ymin = estimate - se, ymax = estimate + se),
                position = dodge, width = 0.3)+ 
  theme_bw() +
  labs(x = '', y = 'Latency (ms)', 
       title = 'Case means',
       fill = 'Case', face = 'bold', size = 20)+
  coord_cartesian(ylim = c(90, 170)) +
  theme(legend.position = 'none',
        legend.text = element_text(face = 'bold', size = 16),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15, face = 'bold'),
        axis.text.y = element_text(size = 15,face = 'bold'),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 16, face = 'bold'),
        axis.title.y = element_text(size = 16, face = 'bold'),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = 20, face = 'bold')) + 
  scale_fill_manual(values = cbPalette)

# RH 500 Hz LMM -----------------------------------------------------------

#RH, 500 Hz response only

rh.500.lmm <-
  lmer(M100LatCorr ~ Case + Age_Calc + Site + (1|Subject),
       data = child.data, subset = Hem == '2-RH' & Cond == '3-500',
       control = lmerControl(check.nobs.vs.nlev = 'ignore',
                             check.nobs.vs.nRE = 'ignore'),
       na.action = na.exclude, REML = F)

rh.500.anova <- Anova(rh.500.lmm)

rh.500.lmm.tukey <- 
  summary(glht(rh.500.lmm,linfct=mcp(Case='Tukey',covariate_average=T)))


# run residual analyses ---------------------------------------------------

path <- 
  '~/GitHub/chromosomal-mutations-auditory-latency/scripts/'

residual.analysis.file <- 
  file.path(path, 
        'child-control-deletion-duplication-residual-analysis-december-2014.R')

library(grid)
source(residual.analysis.file, echo = T)


# run ASD analyses --------------------------------------------------------

asd.analysis.file <- 
  file.path(path, 
        'child-control-deletion-duplication-asd-status-analysis-december-2014.R')

source(asd.analysis.file, echo = T)


# run auxiliary analyses --------------------------------------------------

auxiliary.analysis.file <- 
  file.path(path, 
        'child-control-deletion-duplication-auxiliary-analyses-december-2014.R')

source(auxiliary.analysis.file, echo = T)


# multipanel data plots ---------------------------------------------------

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# least-squres means + age plot -------------------------------------------

multiplot(fitted.vs.age.plot,cond.means.plot, hem.means.plot,
          layout = matrix(c(1, 1, 3, 4), nrow = 2, byrow = T), byrow = T)

print(cond.means.plot,
vp = viewport(layout.pos.row = 2, layout.pos.col = 1))

print(hem.means.plot,
vp = viewport(layout.pos.row = 2, layout.pos.col = 2))


# density estimate plot ---------------------------------------------------


# placing density estimate plots into a single file
grid.newpage()

pushViewport(viewport(layout = grid.layout(2, 1)))

print(case.means.plot,
      vp = viewport(layout.pos.row = 1,layout.pos.col = 1))

print(residual.plus.mean.density,
      vp = viewport(layout.pos.row = 2,layout.pos.col = 1))


# saving data -------------------------------------------------------------

# save analysis results (models, data, plots, etc.)
# 
# save.path <- 
#   '~/Dropbox/cerebral-cortex-submission/analysis-scripts/resubmission-analyses'
# 
# save.file <- 
#   'child-control-deletion-duplication-analyses-workspace-december-2014.rda'
# 
# save.data <- file.path(save.path, save.file)
# 
# save.image(save.data)
