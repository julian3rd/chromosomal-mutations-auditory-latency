
# subjects w/M100 longer than 185ms ---------------------------------------

# subjects with M100s longer than 185ms (corrected)
long.m100.vec <- which(child.data$M100LatCorr > 185)
long.subjs <- child.data[long.m100.vec, ]
long.subjs <- droplevels(long.subjs)

long.subjs.wide <- 
  dcast(data = long.subjs,
        formula = Subject + Age_Calc + Case + ASD + NVIQ ~ Hem + Cond,
        value.var = 'M100LatCorr')


# subjects with at least one measured M100 --------------------------------

# subjects with at least one measured M100
child.data.real.m100 <- subset(child.data, M100LatCorr != 'NA')
child.data.real.m100 <- droplevels(child.data.real.m100)

child.data.real.m100.wide <- 
  dcast(data = child.data.real.m100, 
        formula = Subject + Site + Handedness + Age_Calc + Case + ASD + NVIQ + VIQ + 
          Gender + CELF.4 + SRS_parent + CTOPP + cmICV + dB.SL ~ Hem + Cond, 
        value.var = 'M100LatCorr')

control.data.real <- which(child.data.real.m100.wide$Case == 'control')
deletion.data.real <- which(child.data.real.m100.wide$Case == 'deletion')
duplication.data.real <- which(child.data.real.m100.wide$Case == 'duplication')


# ANOVA: Age, psych for at least M100 -------------------------------------

age.aov <- aov(Age_Calc ~ Case, data = child.data.real.m100.wide)
nviq.aov <- aov(NVIQ ~ Case, data = child.data.real.m100.wide)
viq.aov <- aov(VIQ ~ Case, data = child.data.real.m100.wide)
celf.aov <- aov(CELF.4 ~ Case, data = child.data.real.m100.wide)
srs.aov <- aov(SRS_parent ~ Case, data = child.data.real.m100.wide)
ctopp.aov <- aov(CTOPP ~ Case, data = child.data.real.m100.wide)


# psych score distributions -----------------------------------------------

# all data
child.psych.score.summary <- 
  ddply(child.wide.analysis, .(Case), summarize,
      n = length(Case),
      meanAge = mean(Age_Calc, na.rm = T), sdAge = sd(Age_Calc, na.rm = T),
      meanNVIQ = mean(NVIQ, na.rm = T), sdNVIQ = sd(NVIQ, na.rm = T),
      meanVIQ = mean(VIQ, na.rm = T), sdVIQ = sd(VIQ, na.rm = T),
      meanCELF = mean(CELF.4, na.rm = T), sdCELF = sd(CELF.4, na.rm = T),
      meanSRS = mean(SRS_parent, na.rm = T), sdSRS = sd(SRS_parent, na.rm = T),
      meanCTOPP = mean(CTOPP, na.rm = T), sdCTOPP = sd(CTOPP, na.rm = T))

# with at least one M100
child.real.psych.score.summary <- 
  ddply(child.data.real.m100.wide, .(Case), summarize,
        n = length(Case),
        meanAge = mean(Age_Calc, na.rm = T), sdAge = sd(Age_Calc, na.rm = T),
        meanNVIQ = mean(NVIQ, na.rm = T), sdNVIQ = sd(NVIQ, na.rm = T),
        meanVIQ = mean(VIQ, na.rm = T), sdVIQ = sd(VIQ, na.rm = T),
        meanCELF = mean(CELF.4, na.rm = T), sdCELF = sd(CELF.4, na.rm = T),
        meanSRS = mean(SRS_parent, na.rm = T), sdSRS = sd(SRS_parent, na.rm = T),
        meanCTOPP = mean(CTOPP, na.rm = T), sdCTOPP = sd(CTOPP, na.rm = T))


# Case-ASD analysis -------------------------------------------------------

# Case-ASD analysis
child.data2 <- child.data

child.data2$Case.ASD <- 
  as.factor(paste(child.data2$Case, child.data2$ASD, sep = '-'))

child.case.asd.addmodel <- 
  lmer(M100LatCorr ~ Hem + Cond + Case.ASD + Age + Site + (Cond + Hem|Subject), 
       data = child.data2, na.action = na.exclude, REML = F)

child.case.asd.addmodel.fortify <- fortify(child.case.asd.addmodel)

child.case.asd.addmodel.fortify$residPlusMean <- NA

control.vec <- 
  which(child.case.asd.addmodel.fortify$Case.ASD == 'control-control')

deletion.false.vec <- 
  which(child.case.asd.addmodel.fortify$Case.ASD == 'deletion-FALSE')

deletion.true.vec <- 
  which(child.case.asd.addmodel.fortify$Case.ASD == 'deletion-TRUE')

duplication.false.vec <- 
  which(child.case.asd.addmodel.fortify$Case.ASD == 'duplication-FALSE')

duplication.true.vec <- 
  which(child.case.asd.addmodel.fortify$Case.ASD == 'duplication-TRUE')


# Case-ASD means and plot -------------------------------------------------

# case-ASD means and plot
child.case.asd.means.table <- 
  doBy::LSmeans.lmerMod(child.case.asd.addmodel, effect = c('Case.ASD'))

case.asd.means <-
  cbind(child.case.asd.means.table$grid, child.case.asd.means.table$coef)

case.asd.means$Case.ASD <- as.factor(case.asd.means$Case.ASD)

case.asd.means$Case.ASD <- 
  revalue(case.asd.means$Case.ASD, c('control-control' = 'control'))

case.asd.means$Case <- case.asd.means$Case.ASD
case.asd.means$ASD <- case.asd.means$Case.ASD

case.asd.means$Case <- 
  revalue(case.asd.means$Case, 
          c('deletion-FALSE' = 'deletion', 'deletion-TRUE' = 'deletion',
          'duplication-TRUE' = 'duplication', 'duplication-FALSE' = 'duplication'))

case.asd.means$ASD <- 
  revalue(case.asd.means$ASD, c('deletion-FALSE' = 'FALSE', 'deletion-TRUE' = 'TRUE',
          'duplication-TRUE' = 'TRUE', 'duplication-FALSE' = 'FALSE'))

case.asd.means.colorblind.plot <- 
  ggplot(case.asd.means, 
         aes(x = Case.ASD, y = estimate, fill = Case, group = Case))+
  geom_bar(stat = 'identity', position = 'dodge', size = 1.3, aes(fill = Case))+
  geom_errorbar(size = 1.2, colour = "black", 
                aes(ymin = estimate - se, ymax = estimate + se),
                position = dodge, width = 0.3)+
  theme_bw() +
  labs(x = '', y = 'Latency (ms)',
       title = 'M100 least-squares means\npopulation subdivsion',
       fill = 'Case', colour = 'Case', face = 'bold', size = 20)+
  coord_cartesian(ylim = c(90, 170)) +
  scale_y_continuous(breaks = seq(90, 170, 15)) +
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
        plot.title = element_text(size = 20, 
                                  face = 'bold')) +
  scale_fill_manual(values = cbPalette) +
  scale_x_discrete(labels = c('control', 'deletion', 'deletion',
                              'duplication', 'duplication')) +
  annotate('text', x = c(3, 5), y = c(140, 130), 
           label = 'ASD', face = 'bold', size = 8)


# Case-ASD residuals plus means -------------------------------------------

child.case.asd.addmodel.fortify[control.vec, 'residPlusMean'] <- 
  child.case.asd.addmodel.fortify[control.vec, '.resid'] + 
  case.asd.means[1, 'estimate']

child.case.asd.addmodel.fortify[deletion.false.vec, 'residPlusMean'] <- 
  child.case.asd.addmodel.fortify[deletion.false.vec, '.resid'] + 
  case.asd.means[2, 'estimate']

child.case.asd.addmodel.fortify[deletion.true.vec, 'residPlusMean'] <- 
  child.case.asd.addmodel.fortify[deletion.true.vec, '.resid'] + 
  case.asd.means[3, 'estimate']

child.case.asd.addmodel.fortify[duplication.false.vec, 'residPlusMean'] <- 
  child.case.asd.addmodel.fortify[duplication.false.vec, '.resid'] + 
  case.asd.means[4, 'estimate']

child.case.asd.addmodel.fortify[duplication.true.vec, 'residPlusMean'] <- 
  child.case.asd.addmodel.fortify[duplication.true.vec, '.resid'] + 
  case.asd.means[5, 'estimate']

case.asd.residual.plus.mean.plot <- 
  ggplot(child.case.asd.addmodel.fortify, aes(x = residPlusMean)) +
  geom_density(aes(fill = Case, colour = Case, linetype = Case.ASD),
               size = 1.2, alpha = 0.5) + theme_bw() +
  labs(x = 'Residual + mean values (ms)', y = 'Density',
       title = 'LMM residuals plus means\npopulation subdivision',
       fill = 'Case', colour = 'Case', linetype = 'ASD status')+
  scale_fill_manual(values = cbPalette,
                    labels = c('control', 'deletion', 'duplication'))+
  scale_colour_manual(values = cbPalette,
                      labels = c('control', 'deletion', 'duplication'))+
  scale_linetype_manual(values = c('solid', 'dotted', 'dashed', 
                                  'dotted', 'dashed'),
                        labels = c('control', 'deletion', 'deletion\nASD',
                                   'duplication', 'duplication\nASD'))+
  scale_x_continuous(breaks = seq(90, 170, 20))+
  scale_y_continuous(breaks = seq(0, 0.12, 0.02))+
  coord_cartesian(ylim = c(0, 0.12), xlim = c(90, 170))+
  theme(legend.position = 'bottom',
        legend.text = element_text(face = 'bold', size = 16),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15, face = 'bold'),
        axis.text.y = element_text(size = 15, face = 'bold'),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 20, face = 'bold'),
        axis.title.y = element_text(size = 20, face = 'bold'),
        plot.title = element_text(size = 20, face = 'bold'))


# testing to see if Case-ASD residuals variances differ -------------------

case.asd.levene <- 
  leveneTest(residPlusMean ~ Case.ASD,  
             data = child.case.asd.addmodel.fortify, center = 'median')

# Case-ASD residual plus means KDEs ---------------------------------------

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))

print(case.asd.means.colorblind.plot, 
      vp = viewport(layout.pos.row = 1,layout.pos.col = 1))

print(case.asd.residual.plus.mean.plot, 
      vp = viewport(layout.pos.row = 2, layout.pos.col = 1))



# dB SL: population differences, as covariate -----------------------------


# LMMS: done on participants with evaluable data 
dB.SL.lmm <- 
  lmer(dB.SL ~ Hem + Case + Age_Calc + Site + (1|Subject), 
       data = child.data.real.m100, na.action = na.exclude, REML = F)

dB.SL.aov <- 
  aov(dB.SL ~ Case * Age_Calc + Site, data = child.data.real.m100.wide)

dB.SL.lmm.anova <- Anova(dB.SL.lmm)

dB.SL.lmm.hem.tukey <- 
  summary(glht(dB.SL.lmm, linfct = mcp(Hem = 'Tukey', covariate_average = T)))

dB.SL.lmm.case.tukey <- 
  summary(glht(dB.SL.lmm, linfct = mcp(Case = 'Tukey', covariate_average = T)))

dB.SL.lmm.site.tukey <- 
  summary(glht(dB.SL.lmm, linfct = mcp(Site = 'Tukey', covariate_average = T)))

# LMM: dB.SL as M100 covariate

child.addmodel.dB.SL <- 
  lmer(M100LatCorr ~ Hem + Cond + Case + Age_Calc + Site + 
         dB.SL + (Cond + Hem|Subject),
       data = child.data, na.action = na.exclude, REML = F)


child.addmodel.dB.SL.interact <- 
  lmer(M100LatCorr ~ Hem + Cond + Case + Age_Calc + Site + 
         dB.SL + dB.SL:Case + (Cond + Hem|Subject),
       data = child.data, na.action = na.exclude, REML = F)

child.addmodel.dB.SL.anova <- Anova(child.addmodel.dB.SL)
child.addmodel.dB.SL.interact.anova <- Anova(child.addmodel.dB.SL.interact)

# classification/group separation -----------------------------------------

library(caret)
library(kernlab)

child.addmodel.fortify.mean.collapse.2 <- 
  na.omit(child.addmodel.fortify.mean.collapse)

# separate into train and test sets (50/50 split)
RNGkind('Mersenne-Twister')

set.seed(847508904)

train.rows <- 
  createDataPartition(child.addmodel.fortify.mean.collapse.2$Case, 
                      p = 0.5, list = F)

train.set <- child.addmodel.fortify.mean.collapse.2[train.rows, ]
test.set <- child.addmodel.fortify.mean.collapse.2[-train.rows, ]

# control parameters
cv.ctrl <- 
  trainControl(method = 'repeatedcv', summaryFunction = defaultSummary,
               repeats = 5,  classProbs = T)


# random forest model
rf.grid <- data.frame(.mtry = c(2, 3)) # use sqrt(variables), i.e., 5
rf.grid.no.m100 <- data.frame(.mtry = c(2, 3))
rf.grid.m100.only <- data.frame(.mtry = c(1, 2 ))


set.seed(768505036)

case.rf.tune <- 
  train(Case ~ meanPredicted + NVIQ + VIQ + CELF.4 + SRS_parent,
        data = train.set, method = 'rf', metric = 'Accuracy', 
        tuneGrid = rf.grid, trControl = cv.ctrl)

case.rf.tune.no.m100 <- 
  train(Case ~ NVIQ + VIQ + CELF.4 + SRS_parent,
        data = train.set, method = 'rf', metric = 'Accuracy', 
        tuneGrid = rf.grid.no.m100, trControl = cv.ctrl)

case.rf.tune.m100.only <- 
  train(Case ~ meanPredicted,
        data = train.set, method = 'rf', metric = 'Accuracy', 
        tuneGrid = rf.grid.m100.only, trControl = cv.ctrl)

case.rf.tune

plot(case.rf.tune)

# SVM model
set.seed(768505036)

case.svm.tune <- 
  train(Case ~ meanPredicted + NVIQ + VIQ + CELF.4 + SRS_parent,
        data = train.set, method = 'svmRadial', metric = 'Accuracy',
        tuneLength = 9, preProcess = c('center', 'scale'),
        trControl = cv.ctrl)

case.svm.tune.no.m100 <- 
  train(Case ~ NVIQ + VIQ + CELF.4 + SRS_parent,
        data = train.set, method = 'svmRadial', metric = 'Accuracy',
        tuneLength = 9, preProcess = c('center', 'scale'),
        trControl = cv.ctrl)

case.svm.tune.m100.only <- 
  train(Case ~ meanPredicted,
        data = train.set, method = 'svmRadial', metric = 'Accuracy',
        tuneLength = 9, preProcess = c('center', 'scale'),
        trControl = cv.ctrl)


case.svm.tune

plot(case.svm.tune)

# confusion matrices
rf.predict <- predict(case.rf.tune, test.set)
confusionMatrix(rf.predict, test.set$Case)

rf.predict.no.m100 <- predict(case.rf.tune.no.m100, test.set)
confusionMatrix(rf.predict.no.m100, test.set$Case)

rf.predict.m100.only <- predict(case.rf.tune.m100.only, test.set)
confusionMatrix(rf.predict.m100.only, test.set$Case)


svm.predict <- predict(case.svm.tune, test.set)
confusionMatrix(svm.predict, test.set$Case)

svm.predict.no.m100 <- predict(case.svm.tune.no.m100, test.set)
confusionMatrix(svm.predict.no.m100, test.set$Case)

svm.predict.m100.only <- predict(case.svm.tune.m100.only, test.set)
confusionMatrix(svm.predict.m100.only, test.set$Case)

