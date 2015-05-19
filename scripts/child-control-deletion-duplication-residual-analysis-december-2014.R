
# adding case means to residuals ------------------------------------------

child.addmodel.fortify$residPlusMean <- NA

child.addmodel.fortify.control <- 
  which(child.addmodel.fortify$Case == 'control')

child.addmodel.fortify.deletion <- 
  which(child.addmodel.fortify$Case == 'deletion')

child.addmodel.fortify.duplication <- 
  which(child.addmodel.fortify$Case == 'duplication')

child.addmodel.fortify[child.addmodel.fortify.control, 'residPlusMean'] <- 
  child.addmodel.fortify[child.addmodel.fortify.control, '.resid'] + 
  child.case.lsmeans[1, 'estimate']

child.addmodel.fortify[child.addmodel.fortify.deletion, 'residPlusMean'] <- 
  child.addmodel.fortify[child.addmodel.fortify.deletion, '.resid'] + 
  child.case.lsmeans[2, 'estimate']

child.addmodel.fortify[child.addmodel.fortify.duplication, 'residPlusMean'] <- 
  child.addmodel.fortify[child.addmodel.fortify.duplication, '.resid'] + 
  child.case.lsmeans[3, 'estimate']


# density estimate of residuals plus means and variance test --------------

residual.plus.mean.density <-
  ggplot(child.addmodel.fortify, aes(x = residPlusMean)) +
  geom_density(aes(fill = Case, colour = Case, linetype = Case),
               alpha = 0.7, size = 1.3) + theme_bw() +
  labs(x = 'Residuals + mean values (ms)', y = 'Density',
       title = 'LMM residuals plus Case means',
       fill = 'Case', colour = 'Case', linetype = 'Case')+
  scale_x_continuous(breaks = seq(90, 170, 20))+
  scale_y_continuous(breaks = seq(0, 0.1, 0.02))+
  coord_cartesian(ylim = c (0, 0.11), xlim = c(90, 170))+
  theme(legend.position = 'bottom',
        legend.text = element_text(face = 'bold', size = 16),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15, face = 'bold'),
        axis.text.y = element_text(size = 15, face = 'bold'),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 20, face = 'bold'),
        axis.title.y = element_text(size = 20, face = 'bold'),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = 20, face = 'bold')) + 
  scale_fill_manual(values = cbPalette) + 
  scale_colour_manual(values = cbPalette) + 
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted'))

case.levene <- 
  leveneTest(residPlusMean ~ Case, 
             data = child.addmodel.fortify, center = 'median')


# K-S tests on residuals plus means (Case pairwise) -----------------------

control.resid.plus.mean <- subset(child.addmodel.fortify, Case == 'control')
deletion.resid.plus.mean <- subset(child.addmodel.fortify, Case == 'deletion')
duplication.resid.plus.mean <- subset(child.addmodel.fortify, Case == 'duplication')

control.resid.plus.mean <- droplevels(control.resid.plus.mean)
deletion.resid.plus.mean <- droplevels(deletion.resid.plus.mean)
duplication.resid.plus.mean <- droplevels(duplication.resid.plus.mean)

control.vs.deletion.ks.test <- 
  ks.test(control.resid.plus.mean$residPlusMean, 
        deletion.resid.plus.mean$residPlusMean, exact = F)

control.vs.duplication.ks.test <-
  ks.test(control.resid.plus.mean$residPlusMean, 
        duplication.resid.plus.mean$residPlusMean, exact = F)

deletion.vs.duplication.ks.test <-
  ks.test(deletion.resid.plus.mean$residPlusMean, 
        duplication.resid.plus.mean$residPlusMean, exact = F)



# duplication power calculation (uses residual standard error) ------------

model.effect.sizes <- 
  glht(child.addmodel, linfct = mcp(Case = 'Tukey', covariate_average = T))

model.effect.sizes.fortify <- fortify(model.effect.sizes, aes(lhs, estimate))

residual.moments <- 
  ddply(child.addmodel.case.cond.interact.fortify, .(Case), 
        summarize, residSD = sd(.resid, na.rm = T))

control.duplication.power <- 
  power.t.test(delta = model.effect.sizes.fortify[2, 'estimate'], 
               sd = residual.moments[3, 'residSD'], power = 0.8)

# residual plus mean collapsed --------------------------------------------

# residual plus mean vs. Age
# take mean for each subject across hemisphere and condition
# keep demographic descriptors
# means: residual plus mean, fitted (from model), predicted values

child.addmodel.fortify.mean.collapse <-
  ddply(child.addmodel.fortify,
        .(Subject, Case, Age_Calc, 
          ASD, NVIQ, VIQ, CELF.4, SRS_parent, CTOPP, cmICV),
        summarize, 
        compCases = sum(complete.cases(M100LatCorr)),
        meanFitted = mean(.fitted, na.rm = T),
        meanPredicted = mean(Predicted),
        medianFitted = median(.fitted, na.rm = T),
        medianPredicted = median(Predicted))


# fitted vs age -----------------------------------------------------------


fitted.vs.age.plot <-
  ggplot(child.addmodel.fortify.mean.collapse,
         aes(y = meanPredicted, x = Age_Calc, shape = Case))+
  geom_point(size = 3, aes(colour = Case)) + theme_bw() +
  labs(y = 'Latency (ms)', x = 'Age (years)',
       title = 'Mean predicted M100 vs. age',
       shape = 'Case', linetype = 'Case')+
  geom_smooth(aes(linetype = Case, colour = Case),
              size = 1.3, method = 'lm', se = F)+
  scale_x_continuous(breaks = seq(5, 17, 2))+
  scale_y_continuous(breaks = seq(90, 215, 20))+
  coord_cartesian(ylim = c(90, 190))+
  theme(legend.position = 'bottom',
        legend.text = element_text(face = 'bold', size = 16),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15, face = 'bold'),
        axis.text.y = element_text(size = 15, face = 'bold'),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 20, face = 'bold'),
        axis.title.y = element_text(size = 20, face = 'bold'),
        plot.title = element_text(size = 20, face = 'bold')) + 
  scale_fill_manual(values = cbPalette) + 
  scale_colour_manual(values = cbPalette) +
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted')) +
  annotate('text', x = 6, y = 160, label = 'a', size = 8)

# regression slopes and significances

by(child.addmodel.fortify.mean.collapse,
   child.addmodel.fortify.mean.collapse[, 'Case'],
   function(x) summary(lm(meanPredicted ~ Age_Calc, data = x)))


by(child.addmodel.fortify.mean.collapse,
   child.addmodel.fortify.mean.collapse[, 'Case'],
   function(x) cor.test(~ meanPredicted + Age_Calc, method = 's', data = x))

m100.age.dist.case.aov <- 
  aov(meanPredicted ~ Age_Calc * Case, 
      data = child.addmodel.fortify.mean.collapse)

# fitted vs NVIQ ----------------------------------------------------------

fitted.vs.nviq.plot <- 
  ggplot(child.addmodel.fortify.mean.collapse,
       aes(y = meanPredicted, x = NVIQ, shape = Case)) +
  geom_point(size = 3, aes(colour = Case)) + theme_bw() +
  labs(title = 'Mean predicted M100 vs. NVIQ', y = 'Latency (ms)', x = 'Case',
       shape = 'Case', linetype = 'Case', colour = 'Case') +
  geom_smooth(aes(linetype = Case, colour = Case),
              size = 1.3, method = 'lm', se = F) +
  scale_x_continuous(breaks = seq(50, 130, 20))+
  scale_y_continuous(breaks = seq(90, 190, 20))+
  coord_cartesian(xlim = c(50, 140))+
  theme(legend.position = 'bottom',
        legend.text = element_text(face = 'bold', size = 8),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15, face = 'bold'),
        axis.text.y = element_text(size = 15, face = 'bold'),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 15, face = 'bold'),
        axis.title.y = element_text(size = 15, face = 'bold'),
        plot.title = element_text(size = 15, face = 'bold')) + 
  scale_fill_manual(values = cbPalette) + 
  scale_colour_manual(values = cbPalette) +
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted'))

# linear model statistics for M100 vs. NVIQ

by(child.addmodel.fortify.mean.collapse,
   child.addmodel.fortify.mean.collapse[, 'Case'],
   function(x) summary(lm(meanPredicted ~ NVIQ, data = x)))

by(child.addmodel.fortify.mean.collapse,
   child.addmodel.fortify.mean.collapse[, 'Case'],
   function(x) cor.test(~ meanPredicted + NVIQ, method = 's', data = x))


# fitted vs VIQ -----------------------------------------------------------

fitted.vs.viq.plot <- 
  ggplot(child.addmodel.fortify.mean.collapse,
       aes(y = meanPredicted, x = VIQ, shape = Case)) +
  geom_point(size = 3, aes(colour = Case)) + theme_bw() +
  labs(title = 'Mean predicted M100 vs. VIQ', y = 'Latency (ms)', x = 'Case',
       shape = 'Case', linetype = 'Case', colour = 'Case')+
  geom_smooth(aes(linetype = Case, colour = Case),
              size = 1.3, method = 'lm', se = F)+
  scale_x_continuous(breaks = seq(50, 150, 20))+
  scale_y_continuous(breaks = seq(90, 190, 20))+
  coord_cartesian(xlim = c(50, 150))+
  theme(legend.position = 'bottom',
        legend.text = element_text(face = 'bold', size = 15),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15, face = 'bold'),
        axis.text.y = element_text(size = 15, face = 'bold'),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 15, face = 'bold'),
        axis.title.y = element_text(size = 15, face = 'bold'),
        plot.title = element_text(size = 15, face = 'bold')) + 
  scale_fill_manual(values = cbPalette) + 
  scale_colour_manual(values = cbPalette) +
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted'))

# linear model statistics for M100 vs. VIQ
by(child.addmodel.fortify.mean.collapse,
   child.addmodel.fortify.mean.collapse[, 'Case'],
   function(x) summary(lm(meanPredicted ~ VIQ, data = x)))

by(child.addmodel.fortify.mean.collapse,
   child.addmodel.fortify.mean.collapse[, 'Case'],
   function(x) cor.test(~ meanPredicted + VIQ, method = 's', data = x))

# fitted vs SRS -----------------------------------------------------------

fitted.vs.srs.plot <-
  ggplot(child.addmodel.fortify.mean.collapse,
         aes(y = meanPredicted, x = SRS_parent)) +
  geom_point(size = 3, aes(shape = Case, colour = Case)) + theme_bw() +
  labs(title = 'Mean predicted M100 vs. SRS',
       y = 'Latency (ms)', x = 'SRS', shape = 'Case', colour = 'Case')+
  stat_smooth(size = 1.3, method = 'lm', se = F,
              aes(linetype = Case, colour = Case))+
  scale_x_continuous(breaks = seq(0, 168, 24)) +
  scale_y_continuous(breaks = seq(90, 190, 20)) +
  coord_cartesian(ylim = c(90, 190), xlim = c(-4, 150))+
  theme(legend.position = 'bottom',
        legend.text = element_text(face = 'bold', size = 8),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15, face = 'bold'),
        axis.text.y = element_text(size = 15, face = 'bold'),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 15, face = 'bold'),
        axis.title.y = element_text(size = 15, face = 'bold'),
        plot.title = element_text(size = 15, face = 'bold')) + 
  scale_fill_manual(values = cbPalette) + 
  scale_colour_manual(values = cbPalette) +
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted'))

fitted.vs.sqrt.srs.plot <-
  ggplot(child.addmodel.fortify.mean.collapse,
         aes(y = meanPredicted, x = sqrt(SRS_parent))) +
  geom_point(size = 3, aes(shape = Case, colour = Case)) + theme_bw() +
  labs(title = 'Mean predicted M100 vs. sqrt(SRS)',
       y = 'Latency (ms)', x = 'sqrt(SRS)', shape = 'Case', colour = 'Case')+
  stat_smooth(size = 1.3, method = 'lm', se = F,
              aes(linetype = Case, colour = Case))+
  scale_x_continuous(breaks = seq(0, 14, 2)) +
  scale_y_continuous(breaks = seq(90, 190, 20)) +
  coord_cartesian(ylim = c(90, 190), xlim = c(-1, 14))+
  theme(legend.position = 'bottom',
        legend.text = element_text(face = 'bold', size = 12),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15, face = 'bold'),
        axis.text.y = element_text(size = 15, face = 'bold'),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 15, face = 'bold'),
        axis.title.y = element_text(size = 15, face = 'bold'),
        plot.title = element_text(size = 15, face = 'bold')) + 
  scale_fill_manual(values = cbPalette) + 
  scale_colour_manual(values = cbPalette) + 
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted'))

# linear model statistics for M100 vs. SRS

# untransformed SRS
by(child.addmodel.fortify.mean.collapse,
   child.addmodel.fortify.mean.collapse[, 'Case'],
   function(x) summary(lm(meanPredicted ~ SRS_parent, data = x)))

by(child.addmodel.fortify.mean.collapse,
   child.addmodel.fortify.mean.collapse[, 'Case'],
   function(x) cor.test(~ meanPredicted + SRS_parent, method = 's', data = x))

# sqrt(SRS)
by(child.addmodel.fortify.mean.collapse,
   child.addmodel.fortify.mean.collapse[, 'Case'],
   function(x) summary(lm(meanPredicted ~ sqrt(SRS_parent), data = x)))

by(child.addmodel.fortify.mean.collapse,
   child.addmodel.fortify.mean.collapse[, 'Case'],
   function(x) cor.test(~ meanPredicted + sqrt(SRS_parent), method = 's', data = x))

# fitted vs CTOPP ---------------------------------------------------------

fitted.vs.ctopp.plot <- 
  ggplot(child.addmodel.fortify.mean.collapse,
         aes(y = meanPredicted, x = CTOPP)) +
  geom_point(size = 3, aes(shape = Case, colour = Case)) + theme_bw()+
  labs(title = 'Mean predicted M100 vs. CTOPP',
       y = 'Latency (ms)', x = 'CTOPP (standard score)',
       shape = 'Case', colour = 'Case')+
  stat_smooth(size = 1.3, method = 'lm', se = F,
              aes(linetype = Case, colour = Case))+
  scale_x_continuous(breaks = seq(0, 12, 2))+
  scale_y_continuous(breaks = seq(90, 215, 20))+
  theme(legend.position = 'none',
        legend.text = element_text(face = 'bold', size = 15),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15, face = 'bold'),
        axis.text.y = element_text(size = 15, face = 'bold'),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 15, face = 'bold'),
        axis.title.y = element_text(size = 15, face = 'bold'),
        plot.title = element_text(size = 15, face = 'bold')) + 
  scale_fill_manual(values = cbPalette) + 
  scale_colour_manual(values = cbPalette) +
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted'))

# linear model statistics for M100 vs. CTOPP
by(child.addmodel.fortify.mean.collapse,
   child.addmodel.fortify.mean.collapse[, 'Case'],
   function(x) summary(lm(meanPredicted ~ CTOPP, data = x)))

by(child.addmodel.fortify.mean.collapse,
   child.addmodel.fortify.mean.collapse[, 'Case'],
   function(x) cor.test(~ meanPredicted + CTOPP, method = 's', data = x))

# fitted VS ICV -----------------------------------------------------------


fitted.vs.icv.overall.plot <- 
  ggplot(child.addmodel.fortify.mean.collapse,
         aes(y = meanPredicted, x = cmICV))+
  geom_point(size = 3, aes(shape = Case, colour = Case)) + theme_bw()+
  labs(title = 'Mean M100 vs. ICV',
       y = 'Latency (ms)', x = 'ICV (cm cubed)',
       shape = 'Case')+
  stat_smooth(size = 1.3, method = 'lm', se = F, colour = 'black')+
  scale_x_continuous(breaks = seq(1e3, 2e3, 0.1e3))+
  scale_y_continuous(breaks = seq(90, 215, 20))+
  theme(legend.position = 'bottom',
        legend.text = element_text(face = 'bold', size = 15),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15, face = 'bold'),
        axis.text.y = element_text(size = 15, face = 'bold'),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 15, face = 'bold'),
        axis.title.y = element_text(size = 15, face = 'bold'),
        plot.title = element_text(size = 15, face = 'bold')) + 
  scale_colour_manual(values = cbPalette)

fitted.vs.icv.copies.plot <- 
  ggplot(child.addmodel.fortify.mean.collapse,
         aes(y = meanPredicted, x = cmICV))+
  geom_point(size = 3, aes(shape = Case, colour = Case)) + theme_bw()+
  labs(title = 'Mean M100 vs. ICV',
       y = 'Latency (ms)', x = 'ICV (cm cubed)',
       shape = 'Case', linetype = 'Case', colour = 'Case')+
  stat_smooth(size = 1.3, method = 'lm', se = F,
              aes(linetype = Case, colour = Case))+
  scale_x_continuous(breaks = seq(1e3, 2e3, 0.2e3))+
  scale_y_continuous(breaks = seq(90, 215, 20))+
  theme(legend.position = 'bottom',
        legend.text = element_text(face = 'bold', size = 15),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15, face = 'bold'),
        axis.text.y = element_text(size = 15, face = 'bold'),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 15, face = 'bold'),
        axis.title.y = element_text(size = 15, face = 'bold'),
        plot.title = element_text(size = 15, face = 'bold')) + 
  scale_fill_manual(values = cbPalette) + 
  scale_colour_manual(values = cbPalette) +
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted'))

# linear model statistics for M100 vs. ICV
by(child.addmodel.fortify.mean.collapse,
   child.addmodel.fortify.mean.collapse[, 'Case'],
   function(x) summary(lm(meanPredicted ~ cmICV, data = x)))

by(child.addmodel.fortify.mean.collapse,
   child.addmodel.fortify.mean.collapse[, 'Case'],
   function(x) cor.test(~ meanPredicted + cmICV, data = x))

# fitted vs CELF-4 --------------------------------------------------------

fitted.vs.celf4.plot <- 
  ggplot(child.addmodel.fortify.mean.collapse,
         aes(y = meanPredicted, x = CELF.4, shape = Case)) +
  geom_point(size = 3, aes(colour = Case)) + theme_bw() +
  labs(title = 'Mean M100 vs. CELF-4', y = 'Latency (ms)', x = 'CELF-4',
       shape = 'CELF-4', linetype = 'CELF-4', colour = 'CELF-4') +
  geom_smooth(aes(linetype = Case, colour = Case),
              size = 1.3, method = 'lm', se = F) +
  scale_x_continuous(breaks = seq(40, 140, 20))+
  scale_y_continuous(breaks = seq(90, 190, 20))+
  coord_cartesian(xlim = c(30, 140))+
  theme(legend.position = 'none',
        legend.text = element_text(face = 'bold', size = 15),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15, face = 'bold'),
        axis.text.y = element_text(size = 15, face = 'bold'),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 15, face = 'bold'),
        axis.title.y = element_text(size = 15, face = 'bold'),
        plot.title = element_text(size = 15, face = 'bold')) + 
  scale_fill_manual(values = cbPalette) + 
  scale_colour_manual(values = cbPalette) +
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted'))

# linear model statistics for M100 vs. CELF-4
by(child.addmodel.fortify.mean.collapse,
   child.addmodel.fortify.mean.collapse[, 'Case'],
   function(x) summary(lm(meanPredicted ~ CELF.4, data = x)))

by(child.addmodel.fortify.mean.collapse,
   child.addmodel.fortify.mean.collapse[, 'Case'],
   function(x) cor.test(~ meanPredicted + CELF.4, data = x))


# fitted vs. NVIQ, SRS, CELF-4, CTOPP regressions -------------------------

# plotting NVIQ, SRS, CELF-4, CTOPP regressions
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))

print(fitted.vs.nviq.plot, 
      vp = viewport(layout.pos.row = 1,layout.pos.col = 1))

print(fitted.vs.srs.plot, 
      vp = viewport(layout.pos.row = 1, layout.pos.col = 2))

print(fitted.vs.celf4.plot, 
      vp = viewport(layout.pos.row = 2,layout.pos.col = 1))

print(fitted.vs.ctopp.plot, 
      vp = viewport(layout.pos.row = 2,layout.pos.col = 2))