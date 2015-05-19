#Simons child M100 analysis final version README

analyses done in R 3.1.2

##Packages used:
ggplot2  
lme4  
car  
doBy/lmerTest  
multcomp  
plyr  
reshape2  
caret  
kernlab  


##Files for analysis:

1. child-control-deletion-duplication-analyses-december-2014.R  
2. child-control-deletion-duplication-residual-analysis-december-2014.R  
3. child-control-deletion-duplication-asd-status-december-2014.R  
4. child-control-deletion-duplication-auxiliary-analyses.R  


Each file/script is separated into code sections that can be easily accessed using RStudio. File descriptions list what occurs in each code section; code sections have the format of ‘# description ------‘.  

All objects (LMM calls, tests, plots, etc.) are named so as to describe what is going on as best as possible.  

Analysis results are saved in the workspace named ‘child-control-deletion-duplication-analyses-worspace-december-2014.rda’  


Description of files:  

##child-control-deletion-duplication-analyses-december-2014.R

This is the main analysis file. It should be noted that in the multiplot() section to create Fig 2, there is a bug that results in a Traceback error. If you run everything up until that point (line 644) and then run everything afterward, output will be fine.  

Contents (in order):  
▪	Subset of data to analyze (only validated data)  
▪	setting age cutoff (over/under 12 years, yearly brackets)  
▪	Case subsets (control, deletion duplication)  
▪	recruited and analyzed data converted to wide format  
▪	color-blind friendly palette definition  
▪	loading ggplot2, lme4; Site comparison LMM 
▪	LMMs: all the various models considered in the analysis; primary models of interest are child.fullmodel and child.addmodel (the former used for the condition and hemisphere mean plots; the latter used in the reporting of the main results)  
▪	appending fitted, residuals and predicted values to the main data frame  
▪	psychological scores as covariates in LMMs  
▪	p-values and significance tables  
▪	success rate calculation and evaluation  
▪	least-squares means and plot for full model  
▪	least-squares means and plot for each condition  
▪	least-squares means and plot for each hemisphere  
▪	least-squares means and plot for each case  
▪	RH 500 Hz LMM  
▪	command to run residual analyze script  
▪	command to run ASD analysis script  
▪	command to run auxiliary analysis script  
▪	function stup to plot Fig 2 in manuscript  
▪	least-saures mean + age plot (Fig 2)  
▪	density estimate plot  
▪	saving data  

##child-control-deletion-duplication-residual-analysis-december-2014.R

This files contains all the analyses related to the residuals  

Contents (in order)  
▪	addition of Case means to residuals  
▪	residual + mean density plot and Brown-Forsythe test of those distributions between Cases  
▪	K-S tests on residuals plus means  
▪	duplication power test (ultimately not used since pooled variance was used in final manuscript)  
▪	collapsing M100 to single value for each Subject (mean and median predicted and fitted values are calculated for each subject)  
▪	predicted vs age plot (says fitted, but y = meanPredicted in the code)  
▪	regressions lopes and significance, ANOVA to test differences between slopes for each Case  
▪	predicted vs NVIQ and significance  
▪	predicted vs VIQ and significance  
▪	predicted vs SRS and significance (transformed and untransformed)  
▪	predicted vs CTOPP and significance  
▪	predicted vs ICV and significance (ICV in centimeters cubed)  
▪	predicted vs CELF-4 and significance  
▪	plots of predicted vs. psych scores on a grid  


##child-control-deletion-duplication-asd-status-analysis-december-2014.R

File contains all the analyses related to clinical ASD status  

Contents (in order)  
▪	creating dataset of controls and pro bands without clinical ASD  
▪	creating dataset of pro bands with ASD  
▪	three way ASD status main effects LMM (one with Age:ASD interaction, one without)  
▪	hemisphere and condition, no case model (commented out and not used)  
▪	ASD status alone LMMs (one with Age:ASD interaction, one without)  
▪	model significance for the above  
▪	main effects and significance for controls and pro bands w/o ASD  
▪	within pro band ASD comparison (main effects model and significances) and ASD status proportion test  
▪	residuals plus means and plot for all three ASD levels  
▪	single subject M100 (taken from ASD status LMM) calculation  
▪	predicted vs age by ASD status plot  

##child-control-deletion-duplication-auxiliary-analyses-december-2014.R

This file contains analyses that were not part of the planned comparisons but were informative or added to the text to make the descriptions and analyses more complete.  

Contents (in order)  
▪	Subjects with M100 latency (corrected) greater than 185 ms (long and wide format)  
▪	subjects with at least one observed M100 (long and wide format); entire dataset and separated out by Case  
▪	Age, psych score ANOVAs examining differences in distribution between each Case  
▪	psych score distributions for each Case (entire analyzed dataset, subjects with at least one M100)  
▪	Case-ASD analysis (models, least-squares means and plots)  
▪	Case-ASD residuals plus means and plot  
▪	CASE-ASD Brown-Forsythe test  
▪	Case-ASD residual plus means least-shares mean and KDE plots on grid  
▪	dB SL analyses; differences between Case groups and as a covariate for M100  
▪	classification using 50/50 split of collapsed (i.e., one M100 per subject) data  