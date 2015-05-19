# control and proband data w/o ASD ----------------------------------------

#controls vs. probands without ASD
# controls and probands w/o clinical ASD
child.data.no.asd <-
  subset(child.data, ASD == 'FALSE'|ASD == 'control')

child.data.no.asd <- droplevels(child.data.no.asd)


# probands with ASD data --------------------------------------------------

# ASD proband data
proband.asd.data <- subset(child.data, ASD == 'TRUE')
proband.asd.data <- droplevels(proband.asd.data)


# three-way ASD LMM -------------------------------------------------------

#three-way analysis: control, no clinical ASD, clinical ASD
asd.status.addmodel.interact <-
  lmer(M100LatCorr ~  Hem + Cond + ASD + Age + Age:ASD + 
         Site + (Cond + Hem|Subject), 
       data = child.data, na.action = na.exclude, REML = F)

asd.status.addmodel <-
  lmer(M100LatCorr ~ Hem + Cond + ASD + Age + Site + (Cond + Hem|Subject), 
       data = child.data, na.action = na.exclude, REML = F)


# hemisphere and condition, no case ---------------------------------------

# asd.status.nocase.addmodel <-
#   lmer(M100LatCorr ~ Hem + Cond + Age + Site + (Cond + Hem|Subject), 
#        data = child.data, na.action = na.exclude, REML = F)


# ASD only, with and without Age interaction ------------------------------

asd.status.only <-
  lmer(M100LatCorr ~ ASD * Age + Site + (Cond + Hem|Subject), 
       data = child.data, na.action = na.exclude, REML = F)

# ASD only, NO interaction with Age
asd.status.only.nointeract <-
  lmer(M100LatCorr ~ ASD + Age + Site + (Cond + Hem|Subject), 
       data = child.data, na.action = na.exclude, REML = F)


# model significance ------------------------------------------------------

# Wald Type II chi-square for factor significance
asd.status.addmodel.interact.anova <- Anova(asd.status.addmodel.interact)
asd.status.addmodel.anova <- Anova(asd.status.addmodel)
asd.status.only.anova <- Anova(asd.status.only)

# effect sizes and significance via Tukey HSD
asd.status.addmodel.interact.tukey <- 
  summary(glht(asd.status.addmodel.interact,
               linfct = mcp(ASD = 'Tukey', covariate_average = T, 
                            interaction_average = T)))

asd.status.only.tukey <- 
  summary(glht(asd.status.only,
               linfct = mcp(ASD = 'Tukey', covariate_average = T, interaction_average = T)))

asd.status.addmodel.tukey <- 
  summary(glht(asd.status.addmodel,
               linfct = mcp(ASD = 'Tukey', covariate_average = T)))


# controls and probands w/o ASD: LMM and significance ---------------------

#main effects and significance: controls and probands w/o ASD
child.no.asd.lmm <-
  lmer(M100LatCorr ~ Hem + Cond + ASD + Age + Site +  (Cond + Hem|Subject), 
       data = child.data.no.asd, na.action = na.exclude, REML = F)

child.no.asd.case.lmm <-
  lmer(M100LatCorr ~ Hem + Cond + Case + Age + Site + (Cond + Hem|Subject), 
       data = child.data.no.asd, na.action = na.exclude, REML = F)

child.no.asd.anova <- Anova(child.no.asd.lmm)


child.no.asd.hem.tukey <- 
  summary(glht(child.no.asd.lmm,
               linfct = mcp(Hem = 'Tukey', covariate_average = T)))

child.no.asd.cond.tukey <- 
  summary(glht(child.no.asd.lmm,
               linfct = mcp(Cond = 'Tukey', covariate_average = T)))

child.no.asd.asd.tukey <- 
  summary(glht(child.no.asd.lmm,
               linfct = mcp(ASD = 'Tukey',
              covariate_average = T, interaction_average = T)))


# within-proband ASD comparison -------------------------------------------

#within-proband ASD comparison (main effects model)
proband.asd.lmm <-
  lmer(M100LatCorr ~ Hem + Cond + Case + Age + Site + (Cond + Hem|Subject),
       data = proband.asd.data, na.action = na.exclude, REML = F)

proband.asd.anova <- Anova(proband.asd.lmm)

proband.asd.hem.tukey <- 
  summary(glht(proband.asd.lmm,
               linfct = mcp(Hem = 'Tukey', covariate_average = T)))

proband.asd.cond.tukey <- 
  summary(glht(proband.asd.lmm,
               linfct = mcp(Cond = 'Tukey', covariate_average = T)))

proband.asd.asd.tukey <- 
  summary(glht(proband.asd.lmm,
               linfct = mcp(Case = 'Tukey', covariate_average = T)))

# ASD status proportion test

proband.data <- 
  subset(child.data, Case == 'deletion'| Case == 'duplication')

proband.data <- droplevels(proband.data)

asd.prop.test <- prop.test(c(10,2), c(36, 13), correct = F, p = NULL)



# residuals plus means for all three ASD levels ---------------------------

# residuals plus means for all three groups (control, no ASD, ASD)
# use main effects model; add in means

# adding in fitted and residual values
asd.status.addmodel.fortify <- fortify(asd.status.addmodel)

asd.status.addmodel.fortify$Predicted <- 
  predict(asd.status.addmodel, newdata = asd.status.addmodel.fortify,
          na.action = na.pass, allow.new.levels = T)

# separating the data according to ASD status (control, TRUE, FALSE)
asd.status.addmodel.fortify.control <- 
  which(asd.status.addmodel.fortify$ASD == 'control')

asd.status.addmodel.fortify.false <- 
  which(asd.status.addmodel.fortify$ASD == 'FALSE')

asd.status.addmodel.fortify.true <- 
  which(asd.status.addmodel.fortify$ASD == 'TRUE')

# getting ASD status means
asd.status.means.table <- 
  doBy::LSmeans.lmerMod(asd.status.addmodel, effect = c('ASD'))

asd.means <- 
  cbind(asd.status.means.table$grid, asd.status.means.table$coef)

asd.means$ASD <- as.factor(asd.means$ASD)

# adding means to residuals

asd.status.addmodel.fortify$residPlusMean <- NA

asd.status.addmodel.fortify[asd.status.addmodel.fortify.control, 'residPlusMean'] <- 
  asd.status.addmodel.fortify[asd.status.addmodel.fortify.control, '.resid'] + 
  asd.means[1, 'estimate']

asd.status.addmodel.fortify[asd.status.addmodel.fortify.false, 'residPlusMean'] <- 
  asd.status.addmodel.fortify[asd.status.addmodel.fortify.false, '.resid'] + 
  asd.means[2, 'estimate']

asd.status.addmodel.fortify[asd.status.addmodel.fortify.true, 'residPlusMean'] <- 
  asd.status.addmodel.fortify[asd.status.addmodel.fortify.true, '.resid'] + 
  asd.means[3, 'estimate']


# density estimae of residuals plus means ---------------------------------

#density estiate of residuals plus mean by ASD status
asd.status.resid.plus.mean.density <-
  ggplot(asd.status.addmodel.fortify, aes(x = residPlusMean))+
  geom_density(aes(fill = ASD, colour = ASD, linetype = ASD),
               alpha = 0.7, size = 1.3) + theme_bw()+
  labs(x = 'Residuals + mean values (ms)', y = 'Density',
       title = 'LMM residual plus ASD status means',
       fill = 'ASD status', colour = 'ASD status', linetype = 'ASD status')+
  scale_fill_manual(labels = c('control', 'no ASD', 'ASD'),
                    values = c('black', 'green3', 'blue1'))+
  scale_colour_manual(labels = c('control', 'no ASD', 'ASD'),
                        values = c('black', 'green3', 'blue1'))+
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted'),
                          labels=c('control', 'no ASD', 'ASD'))+
  scale_x_continuous(breaks = seq(90, 170, 20))+
  scale_y_continuous(breaks = seq(0, 0.10, 0.02))+
  coord_cartesian(ylim = c(0, 0.10), xlim = c(90, 170))+
  theme(legend.position = 'bottom',
        legend.text = element_text(face = 'bold', size = 16),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15, face = 'bold'),
        axis.text.y = element_text(size = 15, face = 'bold'),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 20, face = 'bold'),
        axis.title.y = element_text(size = 20, face = 'bold'),
        plot.title = element_text(size = 20, face = 'bold'))


# collapsing fitted values by ASD status ----------------------------------

# take mean for each subject across hemisphere and condition
# keep demographic descriptors
# means: residual plus mean, fitted (from model), predicted values

asd.status.addmodel.fortify.mean.collapse <-
  ddply(asd.status.addmodel.fortify,
        .(Subject, ASD, Case, Copies, Age_Calc, NVIQ, CELF.4, SRS_parent, cmICV),
        summarize,
        compCases = sum(complete.cases(M100LatCorr)),
        meanFitted = mean(.fitted, na.rm = T),
        meanPredicted = mean(Predicted),
        medianFitted = median(.fitted, na.rm = T),
        medianPredicted = median(Predicted))


# fitted vs age for ASD status --------------------------------------------

# fitted vs. age by ASD status

fitted.vs.age.asd.plot <- 
  ggplot(asd.status.addmodel.fortify.mean.collapse,
       aes(y = meanPredicted, x = Age_Calc))+
  geom_point(aes(shape = ASD, colour = ASD), size = 3) + theme_bw()+
  labs(title = 'Mean predicted M100 vs. Age',
       y = 'Latency (ms)', x = 'Age (years)',
       shape = 'ASD status', colour = 'ASD status',
       fill = 'ASD status', linetype = 'ASD status')+
  geom_smooth(aes(linetype = ASD, colour = ASD),
              size = 1.3, method = 'lm', se = F)+
  scale_x_continuous(breaks = seq(7, 17, 2))+
  scale_y_continuous(breaks = seq(90, 190, 20))+
  coord_cartesian(ylim = c(90, 190))+
  theme(legend.position = 'bottom',
        legend.text = element_text(face = 'bold' ,size = 16),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15, face = 'bold'),
        axis.text.y = element_text(size = 15, face = 'bold'),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 20, face = 'bold'),
        axis.title.y = element_text(size = 20, face = 'bold'),
        plot.title = element_text(size = 20, face = 'bold'))+
  scale_shape_manual(values = c(16, 15, 18),
                     labels = c('control', 'no ASD', 'ASD')) +
  scale_colour_manual(values = c('black', 'green3', 'blue1'),
                      labels = c('control', 'no ASD', 'ASD')) +
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted'),
                        labels = c('control', 'no ASD', 'ASD'))