##R code for leachate project plotting and data analysis##
##Created by Jim Junker, 10-14-2014##

#install and load necessary packages#
#only run once if you don't have the packages
#install.packages(c("reshape2", "plyr", "ggplot2", "segmented", "Rcpp"))

#### Load in packages ####
#run everytime
library(plyr)
library(tidyverse)
library(reshape2)
library(segmented)
library(Rcpp)
library(gridExtra)
library(GGally)
library(data.table)
library(corrplot)
library(PerformanceAnalytics)
library(gtable)
library(ggrepel)
library(grid)
library(egg)
library(quantmod)
library(nlme)
library(lme4)

#set plotting themes and color palette#
cbbPalette <- c("#006600", "#FF9900", "#990099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#####
##load in data sheet for DOC timeseries and cumulative measures ####

doc<-read.csv(file="./data-files/Leachate_TimeSeries1.csv",T)
colnames(doc)[2] = "dTime" #rename time between samples column to 'dTime'
#separate dfs for each leachate type
pine<-doc[which(doc$Ctype=="Pine"),]
leaves<-doc[which(doc$Ctype=="Leaves"),]
grass<-doc[which(doc$Ctype=="Grasses"),]
#load in respiration files
resp2 <- read.csv(file = "./data-files/Leachate_Respiration2.csv", T)
resp <-read.csv(file="./data-files/Leachate_Respiration1_mod.csv",T)
#o2 consumption files
o2_resp <- read.csv(file = "./data-files/2018-06-30_O2-consump.csv", T)
#data table of initial metrics for leachates.
doc_tot <- read.csv(file="./data-files/Leachate_CumulativeJan2018.csv",T)

## standardizes to % of carbon loss
doc$log.Per <- log(doc$PerDOCrem)
#doc$time.mod <- doc$Time + 1#not sure we need this
doc$N <- doc$nitrite + doc$nitrate
N_fix = which(doc$N == 0)
doc[N_fix, 'N'] = 0.01
P_fix = which(doc$phosphate == 0)
doc[P_fix, 'phosphate'] = 0.05
doc$NP = (doc$N/doc$phosphate)*(30.974/14.007)
doc$DOC.DIN = (doc$DOC/doc$N) *(14.007/12.011)
doc$DOC.SRP = (doc$DOC/doc$phosphate) *(30.974/12.011)

#NP_fix = which(doc$NP) == T) 
#doc[NP_fix, 'NP'] = -Inf

#DIN_fix = which(is.nan(doc$DOC.DIN) == T) 
#doc[DIN_fix, 'DOC.DIN'] = -Inf

#SRP_fix = which(is.nan(doc$DOC.SRP) == T) 
#doc[SRP_fix, 'DOC.SRP'] = -Inf

doc$doc.cell <- doc$PerdDOC/(doc$Ncells/1000000)/doc$dTime
labile_avg = c()
for(i in 2:nrow(doc)){
  print(doc$labile[i])
  print(doc$labile[i-1])
  labile_avg[i-1] = mean(c(doc$labile[i],doc$labile[i-1]))
  print(labile_avg[i-1])
}
labile_avg
labile_avg = c(NA, labile_avg)

doc$labile_mean = labile_avg
rm(labile_avg)
doc_gr2 <- doc[which(doc$group == 2),]

doc.change <- doc %>% group_by(Ctype) %>% mutate(dDOC = c(NA, diff(PerDOCrem)))
doc.change <- doc.change %>% group_by(Ctype) %>% mutate(U = abs((log(PerDOCrem)-log(shift(PerDOCrem,n=1,type = "lag")))/dTime))
doc.change <- doc.change %>% mutate(Ub = U/(Ncells))

#Figure 1 for publication Respiration and cell count latticed, with [DOC] inset #######
###Analyses on respiration data##
#2018-07-24 NEW ORDER and new size
# - [doc]
# - rates
# - CO2 respiration
# - microbial abundances
resp <- resp[which(resp$Ctype != "Pine2"),]
theme_set(theme_bw(16))
#This is the final plot
F_resp<- ggplot(data=resp, aes(x=Days, y=CO2.accum, group=Ctype, colour=Ctype)) +
  geom_rect(xmin = 6.31, xmax = 7.66, ymin = -50, ymax = 850, fill = "grey",alpha = 0.9, colour = NA) +
  geom_vline(xintercept = 7, linetype = "dashed", size = 1.5) +
  geom_errorbar(aes(ymin = CO2.accum - stdev.CO2.accum, ymax = CO2.accum + stdev.CO2.accum), width = NA, size = 1.5) +
  geom_line(size=1.5) + 
  geom_point(aes(shape = Ctype, fill = Ctype), colour = "black", size = 4, stroke = 1.1) +
  ylab(bquote(~CO[2]~ 'accumulation ('*mu*'L)')) +
  xlab(expression(paste("Time (days)"))) +
  scale_shape_manual(name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                                         labels = c("Grasses", "Leaves", "Pine"), values = c(21,24,22)) +
  scale_fill_manual( name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                     labels = c("Grasses", "Leaves", "Pine"), values = cbbPalette) +
  scale_colour_manual( name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                       labels = c("Grasses", "Leaves", "Pine"), values = cbbPalette) +
  scale_x_continuous(limits = c(0,45), expand = c(0.008,0)) +
  scale_y_continuous(limits = c(-5,810), expand = c(0.04,0)) + coord_cartesian(ylim = c(0,810))+
  #annotate("text", x = 2, y = 770, label = "(b)")+
  theme(legend.position = "none", axis.line.x  = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1), axis.text.x = element_blank(), axis.text.y = element_text(), legend.title = element_text(vjust = 2),
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), panel.border = element_blank());F_resp

##Linearlize and plot the decay over time. 
#this is the final plot for [DOC]
F.lndoc<-ggplot(data=doc, aes(x=Time, y=log(PerDOCrem), group=Ctype, colour=Ctype)) +
  geom_rect(xmin = 6.31, xmax = 7.66, ymin = -50, ymax = 850, fill = "grey",alpha = 0.9, colour = NA) +
  geom_vline(xintercept = 7, linetype = "dashed", size = 1.5) +
  geom_line(size=1.2) + 
  geom_point (aes(fill = Ctype, shape = Ctype),size = 4, colour = "black", stroke = 1.1) +
  ylab(expression(paste("Ln(DOC remaining [%])"))) +
  xlab(expression(paste("Time [days]"))) +
  scale_shape_manual(name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                    labels = c("Grasses", "Leaves", "Pine"), values = c(21,24,22)) +
  scale_colour_manual(name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                    labels = c("Grasses", "Leaves", "Pine"), values = cbbPalette) +
  scale_fill_manual(name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                    labels = c("Grasses", "Leaves", "Pine"), values = cbbPalette) +
  #annotate("text", x = 2, y = 4.8, label = "(a)") +
  scale_x_continuous(limits = c(0,45), expand = c(0.008,0)) +
  scale_y_continuous(limits = c(2.6,4.9))+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.line = element_line(size = 1), axis.title.y = element_text(margin = margin(0,7,0,0)),
      axis.title.x = element_blank(), axis.line.x  = element_blank(),axis.text.x = element_blank(),
      axis.ticks.x = element_blank(), panel.border = element_blank());F.lndoc

####

#DOC uptake 
Uplot = ggplot(doc.change, aes(y = log(U), x = PerDOCrem, group = Ctype)) + 
  geom_path(aes(colour = Ctype), size = 1.5) +
  geom_point(aes(shape = Ctype, fill = Ctype), size = 2.5, stroke = 1)  +
  ylab(bquote('Ln(U [' ~d^-1* '])')) +
  xlab(expression(paste("DOC remaining (%)"))) +
  scale_shape_manual(name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                     labels = c("Grasses", "Leaves", "Pine"), values = c(21,24,22)) +
  scale_fill_manual( name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                     labels = c("Grasses", "Leaves", "Pine"), values = cbbPalette) +
  scale_colour_manual( name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                       labels = c("Grasses", "Leaves", "Pine"), values = cbbPalette) +
  scale_x_continuous(limits = c(0,100), expand = c(0.009,0)) +
  scale_y_continuous(limits = c(-7,-0.01), expand = c(0.009, 0)) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_text(size = 14),        
        plot.background = element_rect( fill = "transparent", colour = NA), axis.title.y = element_text(size = 14,vjust = 2));Uplot

#### This sets a custom annotation with the plot that allows to save as an object
g = ggplotGrob(Uplot)

F1 = F.lndoc + annotation_custom(grob = g, xmin = 18, xmax = 43.8,
                                ymin = 3.45, ymax = 5);F1
print(F1)

## Cells over time multiple looks
Ncell.fin<- ggplot( data = doc, aes( x = Time, y = log10(Ncells), group = Ctype, colour = Ctype)) +
  geom_rect(xmin = 6.31, xmax = 7.66, ymin = -50, ymax = 850, fill = "grey",alpha = 0.9, colour = NA) +
  geom_vline(xintercept = 7, linetype = "dashed", size = 1.5) +
  geom_line (size = 1.5) +
  geom_point (aes(fill = Ctype, shape = Ctype),size = 4, colour = "black", stroke = 1.1) +
  ylab(bquote(~Log[10]*'(Cells ['~mL^-1* '])')) +	
  xlab("Time (days)")+
  scale_shape_manual(name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                       labels = c("Grasses", "Leaves", "Pine"), values = c(21,24,22)) +
  scale_fill_manual( name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                     labels = c("Grasses", "Leaves", "Pine"), values = cbbPalette) +
  scale_colour_manual( name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                       labels = c("Grasses", "Leaves", "Pine"), values = cbbPalette) +
  ylab(bquote(~Log[10]*' (cells [' ~mL^-1* ']) ')) +
  xlab("Time (days)") +
  #annotate("text", x = 2, y = 5.2, label = "(c)") +
  scale_x_continuous(limits = c(0,45), expand = c(0.008,0)) +
  scale_y_continuous(limits = c(5,8.3), expand = c(0.008, 0)) +
theme(legend.position = c(0.85,0.8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(), axis.text.y = element_text(), axis.title.y = element_text(margin = margin(0,15,0,0)),
      panel.border = element_blank(), axis.line = element_line(size = 0.75));Ncell.fin

##### build final plot using ggarrange######
#ggarrange(F1, F_resp, Ncell.fin,ncol = 1)

tiff("Figure1.tiff", res = 600, height = 234, width = 174, units = "mm", compression = "lzw")
ggarrange(F1, F_resp, Ncell.fin, ncol = 1)
dev.off()
####### Done with Figure 1 #####

#################playing with efficiency of DOC respiration#################
#2018-07-24 New size
doc.change_gr2 <- doc.change[which(doc.change$group == 2),]
theme_set(theme_bw(14))
set.seed(101)
Feff <- ggplot( data = doc.change, aes( x = X.labile, y = Ub*1e06)) +
  stat_smooth(data = doc.change_gr2, method = "lm", se = F, alpha = 0.6, size = 1.5, colour = "black", linetype = 1) +
  geom_point(size = 3.5, aes(shape = Ctype, fill = Ncells), colour = "black", stroke = 1.1) +
  scale_shape_manual(name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                     labels = c("Grasses", "Leaves", "Pine"), values = c(21,24,22)) +
  scale_fill_gradient(low = "blue", high = 'red', name = "Cell count") +
  xlab(expression('%'~Lability[EEMs]*'')) +
  ylim(c(-10,20)) + coord_cartesian(ylim = c(0,0.125))+
  #ylab(expression(atop("DOM Processing Efficiency", paste("(U/"*10^6~"cells)")))) +
  ylab(expression("DOM processing efficiency (U/"*10^6~"cells)"))+
  scale_x_reverse() +
  geom_rug(sides = 't')+
  geom_text_repel(data = doc, aes(x = X.labile, y = 0.125, label = Time), size = 2.2, ylim = c(0.115, 0.128), direction = "y", point.padding = NA, segment.colour = "grey", force = 1.4) +
  theme(legend.position = c(0.15,0.53),panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12), 
        axis.title.x = element_text(), axis.title.y = element_text(margin = margin(0,0,0,0)),
        legend.title = element_text(size = 11), legend.text = element_text(size = 10),
        legend.background = element_rect( fill = "transparent", colour = NA));Feff
tiff("Figure3.tiff", width = 129, height = 129, units = "mm", res = 600, compression = "lzw")
Feff
dev.off()
cor.test(doc.change_gr2$X.labile, doc.change_gr2$Ub*1e06)
#######  End cell DOC efficiency #####

###### Compare breakpoint in [DOC] #####
plot(log(PerDOCrem) ~ Time, data = doc[doc$Ctype == "Pine",]) #plot the data

pine.lm <- glm(log(PerDOCrem ) ~  0 + Time, offset = rep(4.60517, length(Time)),data = doc[doc$Ctype == "Pine",])
pine_offset.lm <- glm(log(PerDOCrem ) ~ Time,offset = rep(4.60517, length(Time)),data = doc[doc$Ctype == "Pine",]);summary(pine_offset.lm)
plot(pine.lm);summary(pine.lm)  #plotting the linear model. Look @ standardized residuals and fitted values to estimate time of possible breakpoints
set.seed(123);seg.pine1 = segmented(pine.lm, seg.Z = ~Time, psi = 1)
summary(seg.pine1)
seg.pine1$psi[[2]]
######  Quick test to check sensitivity of breakpoint to starting value #####
seg.pine3 = segmented(pine.lm, seg.Z = ~Time, psi = 3)
summary(seg.pine3)
seg.pine5 = segmented(pine.lm, seg.Z = ~Time, psi = 5)
summary(seg.pine5)
seg.pine15 = segmented(pine.lm, seg.Z = ~Time, psi = 15)
summary(seg.pine15)
seg.pine40 = segmented(pine.lm, seg.Z = ~Time, psi = 40)
summary(seg.pine40)

start = c(1,3,5,15,40)
est = c(seg.pine1$psi[[2]], seg.pine3$psi[[2]], seg.pine5$psi[[2]], seg.pine15$psi[[2]],
        seg.pine40$psi[[2]])
est.se = c(seg.pine1$psi[[3]], seg.pine3$psi[[3]], seg.pine5$psi[[3]], seg.pine15$psi[[3]],
           seg.pine40$psi[[3]])
df = data.frame(start, est, est.se)
limits = c(ymin = est - est.se, ymax = est + est.se)
ggplot(df, aes(x = start, y = est)) + geom_errorbar(aes(ymin = est - est.se, ymax = est + est.se), width = 0, size = 1) +geom_point(size = 3) +
  ylim(limits =c(0,20))
##### End test ######
slope(seg.pine1)
plot(seg.pine1)
points(log(PerDOCrem)-4.60517 ~ Time, data = doc[doc$Ctype == "Pine",])
AICcmodavg::AICc(pine.lm, seg.pine1)
#Running with  grasses
plot(log(PerDOCrem) ~ Time, data = doc[doc$Ctype == "Grasses",]) #plot the data
 #set for running the segmented model 
grass.lm <- glm(log(PerDOCrem) ~ 0 + Time, offset = rep(4.60517, length(Time)), data = doc[doc$Ctype == "Grasses",]) #the Linear model to use for breakpoint estimation
grass_offset.lm <- glm(log(PerDOCrem ) ~ Time,offset = rep(4.60517, length(Time)),data = doc[doc$Ctype == "Grasses",]);summary(grass_offset.lm)
plot(grass.lm);summary(grass.lm)  #plotting the linear model. Look @ standardized residuals and fitted values to estimate time of possible breakpoints
set.seed(123);seg.grass1 = segmented(grass.lm, seg.Z = ~Time, psi = 10) #Piecewise with single breakpoint
summary(seg.grass1)
slope(seg.grass1)

plot(seg.grass1)
points(log(PerDOCrem)-4.60517 ~ Time, data = doc[doc$Ctype == "Grasses",])
AIC(grass.lm, seg.grass1)
##Running the leaves
plot(log(PerDOCrem) ~ Time, data = doc[doc$Ctype == "Leaves",]) #plot the data
 #set for running the segmented model 
leaves.lm <- glm(log(PerDOCrem) ~ 0 + Time, offset = rep(4.60517, length(Time)),data = doc[doc$Ctype == "Leaves",]) #the Linear model to use for breakpoint estimation
leaves_offset.lm <- glm(log(PerDOCrem ) ~ Time,offset = rep(4.60517, length(Time)),data = doc[doc$Ctype == "Leaves",]);summary(leaves_offset.lm)
HEAD
plot(leaves.lm);summary(leaves.lm)  #plotting the linear model. Look @ standardized residuals and fitted values to estimate time of possible breakpoints
755be2983a2adfe9092d06a83579b880f601a9e1
set.seed(123);seg.leaves1 = segmented(leaves.lm, seg.Z = ~Time, psi = 10) #Piecewise with single breakpoint
summary(seg.leaves1)

plot(seg.leaves1)
points(log(PerDOCrem)-4.60517 ~ Time, data = doc[doc$Ctype == "Leaves",])
AIC(leaves.lm, seg.leaves1)
slope(seg.leaves1)
#### DOC psi est ####
seg.pine1$psi[2]
seg.grass1$psi[2]
seg.leaves1$psi[2]

mean(c(seg.pine1$psi[2],seg.grass1$psi[2],seg.leaves1$psi[2]))
##### 
#### DOC k coeffs ####
seg.pine1$coefficients[1]
100-exp(4.60517+seg.pine1$coefficients[1]*seg.pine1$psi[2])
seg.grass1$coefficients[1]
100-exp(4.60517+seg.grass1$coefficients[1]*seg.grass1$psi[2])
seg.leaves1$coefficients[1] 
100-exp(4.60517+seg.leaves1$coefficients[1]*seg.leaves1$psi[2])

#### look at  
short.DOC<-ggplot(data=doc, aes(x=Time, y=log(PerDOCrem), group=Ctype, colour=Ctype)) +
  geom_line(size=1.2) + 
  geom_point (aes(fill = Ctype, shape = Ctype),size = 4, colour = "black", stroke = 1.1) +
  geom_vline(xintercept = 7, linetype = "dashed", size = 1.5) +
  ylab(expression(paste("Ln(DOC remaining [%])"))) +
  xlab(expression(paste("Time [days]"))) +
  scale_shape_manual(name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                     labels = c("Grasses", "Leaves", "Pine"), values = c(21,24,22)) +
  scale_colour_manual(name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                      labels = c("Grasses", "Leaves", "Pine"), values = cbbPalette) +
  scale_fill_manual(name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                    labels = c("Grasses", "Leaves", "Pine"), values = cbbPalette) +
  #annotate("text", x = 43, y = 4.55, label = "(c)") +
  scale_x_continuous(limits = c(0,45), expand = c(0.008,0)) + coord_cartesian(xlim = c(0,8))+
  scale_y_continuous(limits = c(2.6,4.9))+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1), axis.title.y = element_text(margin = margin(0,7,0,0)),
        axis.title.x = element_blank(), axis.line.x  = element_blank(),axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), panel.border = element_blank());short.DOC

###### End [DOC] Breakpoint comparisons #####

#respiration estimates at short term##
f <- function(u) which.min(abs(as.numeric(resp2$Time) - as.numeric(u)))
resp2 %>%
  group_by(Ctype) %>%
  spread(key = Rep, CO2.accum) -> resp.w

resp.pine = resp.w[which(resp.w$Ctype == "Pine"),]
resp.grass = resp.w[which(resp.w$Ctype == "Grasses"),]
resp.leaves = resp.w[which(resp.w$Ctype == "Leaves"),]

resp.pine[f(seg.pine1$psi[2]),];mean(c(467,467,476))
resp.grass[f(seg.grass1$psi[2]),];mean(c(399,359,399))
resp.leaves[f(seg.leaves1$psi[2]),];mean(c(507,556,502))

### Dissolved stoic figure for supplemental ####
CV = function(x){sd(x, na.rm = T)/mean(x, na.rm = T)}
NP_doc = doc[apply(doc[c('Time', 'NP')], 1, function(x) all(is.finite(x))),c('Time', 'Ctype','NP')]
aggregate(NP~Ctype, data = doc, FUN = CV)

NP.plot = ggplot(data = doc, aes(x = Time, y = NP, group = Ctype, colour = Ctype)) +
  geom_path(size = 1.2) +geom_point(aes(shape = Ctype, fill = Ctype), size = 6, colour = "black") + 
  ylab("DIN:SRP (molar)") + xlab("Time (days)") + scale_x_continuous(limits = c(0,44), expand = c(0.01,0)) +
  scale_shape_manual(name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                     labels = c("Grasses", "Leaves", "Pine"), values = c(21,24,22)) +
  scale_fill_manual( name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                     labels = c("Grasses", "Leaves", "Pine"), values = cbbPalette) +
  scale_colour_manual( name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                       labels = c("Grasses", "Leaves", "Pine"), values = cbbPalette) +
  theme(legend.position = c(0.85, 0.85), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.text.y = element_text(size = 22), axis.text.x = element_text(size = 22), 
        axis.title.x = element_text(size = 24), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 0, l = 0, b = 0)));NP.plot

tiff("NP_plot.tiff", res = 600, height = 8, width = 8, units = "in")
NP.plot
dev.off()
N_doc = doc[apply(doc[c('Time', 'DOC.DIN')], 1, function(x) all(is.finite(x))), c('Time', 'Ctype', 'DOC.DIN')]
aggregate(DOC.DIN~Ctype, data = doc, FUN = CV)

DOC_N.plot = ggplot(doc, aes(x = Time, y = DOC.DIN, group = Ctype, colour = Ctype)) +
  geom_point(aes(shape = Ctype, fill = Ctype), size = 6, colour = "black") + geom_path(size = 1.2) +
  ylab("DOC:DIN (molar)") + xlab("Time (days)") +
  scale_shape_manual(name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                     labels = c("Grasses", "Leaves", "Pine"), values = c(21,24,22)) +
  scale_fill_manual( name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                     labels = c("Grasses", "Leaves", "Pine"), values = cbbPalette) +
  scale_colour_manual( name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                       labels = c("Grasses", "Leaves", "Pine"), values = cbbPalette) +
  theme(legend.position = c(0.9, 0.75), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.text.y = element_text(size = 22), axis.text.x = element_blank(), 
        axis.title.x = element_blank(), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 0, l = 0, b = 0)));DOC_N.plot

P_doc = doc[apply(doc[c('Time', 'DOC.SRP')], 1, function(x) all(is.finite(x))), c('Time', 'Ctype', 'DOC.SRP')]
aggregate(DOC.SRP~Ctype, data = doc, FUN = CV)

DOC_P.plot = ggplot(doc, aes(x = Time, y = DOC.SRP, group = Ctype, colour = Ctype)) +
  geom_point(aes(shape = Ctype, fill = Ctype), size = 6, colour = "black") + geom_path(size = 1.2) +
  ylab("DOC:SRP (molar)") + xlab("Time (days)") +
  scale_shape_manual(name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                     labels = c("Grasses", "Leaves", "Pine"), values = c(21,24,22)) +
  scale_fill_manual( name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                     labels = c("Grasses", "Leaves", "Pine"), values = cbbPalette) +
  scale_colour_manual( name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                       labels = c("Grasses", "Leaves", "Pine"), values = cbbPalette) +
  theme(legend.position = 'none', panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.text.y = element_text(size = 22), axis.text.x = element_blank(), 
        axis.title.x = element_blank(), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 0, l = 0, b = 0)));DOC_P.plot

tiff("Dissolved_stoic.tiff", res = 600, height = 15, width = 10, unit = "in")
ggarrange(DOC_N.plot, DOC_P.plot, NP.plot, ncol = 1)
dev.off()

ggarrange(DOC_N.plot, DOC_P.plot, NP.plot, ncol = 1)
grid.arrange(DOC_N.plot, DOC_P.plot, NP.plot, ncol = 1)

###### End Stoic supplemental figure #####

####### Exploratory analysis section #######
###loading in full respiration data #####
str(o2_resp)
#source("./segmented.lme/source-segmented.lme.r")
o2_resp.long = tidyr::gather(o2_resp, Ctype, O2.consum, Grasses:Pine, factor_key = T)
#Build models for each Carbon type
plot(O2.consum ~ Days, data = o2_resp.long[o2_resp.long$Ctype == 'Pine',])
plot(CO2.accum ~ Time, data = resp2[resp2$Ctype == "Pine",]) #plot the data
set.seed(123) #set for running the segmented model
O2.glm = gls(O2.consum ~ Days * Ctype, correlation = corAR1(0.8, form = ~1|Ctype), data = o2_resp.long)
pine.glm1 = gls(O2.consum ~ Days, correlation = corAR1(0.8), data = o2_resp.long[which(o2_resp.long$Ctype == "Pine"),])
pine.glm2 = gls(O2.consum ~ Days, data = o2_resp.long[which(o2_resp.long$Ctype == "Pine"),])
pine.lm <- lm(O2.consum ~ 0 + Days, data = o2_resp.long[o2_resp.long$Ctype == "Pine",])
#pine.lm <- lm(CO2.accum ~ Time, data = resp2[resp2$Ctype == "Pine",])
plot(pine.glm1);summary(pine.glm1)  #plotting the linear model. Look @ standardized residuals and fitted values to estimate time of possible breakpoints
plot(pine.glm2);summary(pine.glm2)
plot(CO2.accum~Time, data = resp2[resp2$Ctype == 'Pine',])
abline(pine.lm)
######  Quick test to check sensitivity of breakpoint to starting value #####
seg.pine1 = segmented(pine.glm1, seg.Z = ~Days, psi = 10)
summary(seg.pine1)
seg.pine1$psi
slope(seg.pine1)
seg.pine2 = segmented(pine.glm1, seg.Z = ~Days, psi = c(1, 40))
summary(seg.pine2)
seg.pine2$psi
slope(seg.pine2)
seg.pine3 = segmented(pine.glm1, seg.Z = ~Days, psi = c(1,3,20))
summary(seg.pine3)
seg.pine3$psi
slope(seg.pine3)
seg.pine5 = segmented(pine.glm1, seg.Z =  ~Days, psi = c(1,5,15,20,37))
seg.pine5$psi
slope(seg.pine5)

start = c(1,2.5,5,10,40)
est = c(seg.pine1$psi[[2]], seg.pine2$psi[[2]], seg.pine5$psi[[2]], seg.pine10$psi[[2]],
        seg.pine40$psi[[2]])
est.se = c(seg.pine1$psi[[3]], seg.pine2$psi[[3]], seg.pine5$psi[[3]], seg.pine10$psi[[3]],
           seg.pine40$psi[[3]])
df = data.frame(start, est, est.se)
limits = c(ymin = est - est.se, ymax = est + est.se)
ggplot(df, aes(x = start, y = est)) + geom_errorbar(aes(ymin = est - est.se, ymax = est + est.se)) +geom_point(size = 1)
##### End test ######

slope(seg.pine1)
plot(seg.pine1)
plot(seg.pine2)
slope(seg.pine2)
points(O2.consum ~ Days, data = o2_resp.long[o2_resp.long$Ctype == "Pine",])
pacf(resid(seg.pine1))
exp(3.016)
AIC(pine.glm, seg.pine1, seg.pine2, seg.pine3, seg.pine5)
#Running with  grasses
plot(O2.consum ~ Days, data = o2_resp.long[o2_resp.long$Ctype == 'Grasses',])
plot(CO2.accum ~ Time, data = resp2[resp2$Ctype == "Grasses",]) #plot the data
set.seed(123) #set for running the segmented model 
grass.lm <- lm(O2.consum ~ 0+Days, data = o2_resp.long[o2_resp.long$Ctype == 'Grasses',]) #the Linear model to use for breakpoint estimation
#grass.lm <- lm(CO2.accum ~ 0 + Time, data = resp2[resp2$Ctype == "Grasses",]) #the Linear model to use for breakpoint estimation
plot(grass.lm);summary(grass.lm)  #plotting the linear model. Look @ standardized residuals and fitted values to estimate time of possible breakpoints

seg.grass1 = segmented(grass.lm, seg.Z = ~Days, psi = 1) #Piecewise with single breakpoint
summary(seg.grass1)
slope(seg.grass1)
seg.grass2 = segmented(grass.lm, seg.Z = ~Days, psi = c(1, 40))
summary(seg.grass2)
seg.grass3 = segmented(grass.lm, seg.Z = ~Days, psi = c(.5,1,3))
summary(seg.grass3)
seg.grass5 = segmented(grass.lm, seg.Z =  ~Days, psi = c(1,5,15,20,37))

plot(seg.grass1)
points(O2.consum ~ Days, data = o2_resp.long[o2_resp.long$Ctype == 'Grasses',])
AIC(grass.lm, seg.grass1, seg.grass2, seg.grass3, seg.grass5)
plot(resid(seg.grass1))
plot(resid(seg.grass2))
plot(resid(seg.grass3))
plot(resid(seg.grass5))

##Running the leaves
plot(CO2.accum ~ Time, data = resp2[resp2$Ctype == "Leaves",]) #plot the data
set.seed(123) #set for running the segmented model 
leaves.lm <- lm(CO2.accum ~ 0+ Time, data = resp2[resp2$Ctype == "Leaves",]) #the Linear model to use for breakpoint estimation
plot(leaves.glm);summary(leaves.glm)  #plotting the linear model. Look @ standardized residuals and fitted values to estimate time of possible breakpoints
seg.leaves1 = segmented(leaves.lm, seg.Z = ~Time, psi = 1) #Piecewise with single breakpoint
summary(seg.leaves1)
slope(seg.leaves1)
exp(2.920)
plot(seg.leaves1)
points(CO2.accum ~ log(Time), data = resp[resp$Ctype == "Leaves",])
pacf(resid(seg.leaves1))

#compare the slope estimates for each regression
slope(seg.pine1)
slope(seg.grass1)
slope(seg.leaves1)
seg.pine1$coefficients[1]

exp(1)^seg.pine1$coefficients[1] -1
exp(1)^seg.grass1$coefficients[1] -1
exp(1)^seg.leaves1$coefficients[1] -1

###### End segmented respiration model #####
###### attempting to build a time step model of michaelis menten kinetics of DOC uptake #####
# Hyperbolic model formula #
# d[DOC]/dt = Vmax[DOC]/Km + [DOC] 

# Can we fit a model that explains d[DOC]/dt with a consistent hyperbolic model?
# fitting a MM model to  DOC loss?
# Just fitting a MM model to DOC change over a day with PerDOC
doc.change <- doc %>% group_by(Ctype) %>% mutate(dDOC = c(NA, diff(PerDOCrem)))
doc.change <- doc.change %>% group_by(Ctype) %>% mutate(U = abs(dDOC/(PerDOCrem*dTime)))
doc.change <- doc.change %>% mutate(Ub = U/(Ncells))
#doc.change = rbind(doc.change, doc.change)
ggplot(doc.change, aes(y = U, x = Time, group = Ctype)) + geom_point(aes(colour = Ctype), size = 4)  +geom_path(aes(colour = Ctype))
ggplot(doc.change, aes(y = U, x = log10(Ncells), group = Ctype)) + geom_point(aes(colour = Ctype), size = 4) + geom_path(aes(colour = Ctype))
ggplot(doc.change, aes(y = U, x = PerDOCrem, group = Ctype)) + geom_point(aes(colour = Ctype), size = 4) + geom_path(aes(colour = Ctype))
ggplot(doc.change, aes(y = log10(Ncells), x = Time, group = Ctype)) + geom_point(aes(colour = Ctype), size = 4) + geom_path(aes(colour = Ctype))
ggplot(doc.change, aes(x = PerDOCrem, y = Ub, group = Ctype)) + geom_point(aes(colour = Ctype), size = 4)
ggplot(doc.change, aes(x = Time, y = Ub, group = Ctype)) + geom_point(aes(colour = Ctype), size = 4) +geom_path(aes(colour = Ctype))
ggplot(doc.change, aes(y = log10(Ncells), x = PerDOCrem, group = Ctype)) + geom_point(aes(colour = Ctype), size = 4)

# fit nls to doc decay
library(nlme)
doc.change$c = NA
pine_fix = which(doc.change$Ctype == "Pine")
doc.change[pine_fix, "c"] = 1
leaves_fix = which(doc.change$Ctype == "Leaves")
doc.change[leaves_fix, "c"] = 2
grass_fix = which(doc.change$Ctype == "Grasses")
doc.change[grass_fix, "c"] = 3

lin_form = formula(U~a*log(Time)*Ncells)#constant uptake
preview(lin_form, data = doc.change, start = )
lin_nls = nls(formula = lin_form, data = doc.change, trace = T, start = list(a = 10))
summary(lin_nls)
mm_form = formula(U~((a*PerDOCrem)/(b+PerDOCrem)))#Uptake rate scaled to DOC
mmcellc_nls = nls(formula = mmcellc_form, data = doc.change, trace = T, start = list(a = a_start, b = b_start))
mmcell_form = formula(log(U)~((a*PerDOCrem)/(b+PerDOCrem)*log10(Ncells)))#Uptake rate scaled to DOC
mmcell_nls = nls(formula = mmcell_form, data = doc.change, trace = T, start = list(a = a_start, b = b_start))
summary(mmcell_nls)
PerDOCrem = seq(10,100, 1)
Ncells = seq(min(doc.change$Ncells), max(doc.change$Ncells), max(doc.change$Ncells)-min(doc.change$Ncells)/91)
predict_mmcell = predict(mmcell_nls, newdata = data.frame(PerDOCrem = PerDOCrem, Ncells = Ncells))
pred_df = data.frame(PerDOCrem, Ncells, predict_mmcell)

#Plotting the data and model
plot(U~PerDOCrem, data = Leaves, pch = 21, col = "black", bg = "orange", cex = 2.5)
points(U~PerDOCrem, data = Pine, pch = 21, col = "black", bg = "purple", cex = 2.5)
points(U~PerDOCrem, data = Grasses, pch = 21, col = "black", bg = "dark green", cex = 2.5 )
lines(U~PerDOCrem, data = Leaves, lwd = 3, col = "orange")
lines(U~PerDOCrem, data = Pine, lwd = 3, col = "purple")
lines(U~PerDOCrem, data = Grasses, lwd = 3, col = "dark green")

lines(predict_mmcell~PerDOCrem, data = pred_df, col = "black", lwd = 3)


mm_form = formula(U~((a*PerDOCrem/(b+PerDOCrem))))#Saturating Uptake scaled to DOC


mmcelltyp_form = formula(U~(isTrue(Ctype == "Pine"))*(a*PerDOCrem/(b+PerDOCrem))*log10(Ncells) +
                           isTrue(Ctype == "Leaves")*(a*PerDOCrem/(b+PerDOCrem))*log10(Ncells) +
                           isTrue(Ctype == "Grasses")*(a*PerDOCrem/(b+PerDOCrem))*log10(Ncells))
mmcelltyp_nls = nls(mmcelltyp_form, data= doc.change, trace = T)
a_start = 0.5
b_start = -20

fit <- nlsLM(U ~ ((a*PerDOCrem)/(b+PerDOCrem)*log10(Ncells)) | c, data=doc.change, start = list(a=1,b=1))

require(broom)
fit_data <- augment(fit)

plot(.fitted~rev, data=fit_data)

mmcellc_form = mle2(U~((a*PerDOCrem)/(b+PerDOCrem)))

plot(doc.change$U[which(doc.change$Ctype == "Leaves")]~doc.change$PerDOCrem[which(doc.change$Ctype == "Leaves")], pch = 21,bg = "orange", xlim = c(10,80), ylim = c(0,0.8), xlab = "PerDOCrem", ylab = "U")
points(doc.change$U[which(doc.change$Ctype == "Pine")]~doc.change$PerDOCrem[which(doc.change$Ctype == "Pine")], pch = 21, bg = "purple")
points(doc.change$U[which(doc.change$Ctype == "Grasses")]~doc.change$PerDOCrem[which(doc.change$Ctype == "Grasses")], pch = 21, bg = "dark green")

set.seed(123)
m = nlsList(formula = mm_form, data = doc.change, start = list(a = a_start, b = b_start))
m = nls(formula = lin_form, data = doc.change, start = list(a = a_start, b = b_start))
m = nls(mmcell_form, data = doc.change, start = list(a = a_start, b = b_start))

summary(m)
logLik(m)
#Doc of 16+ sample
doc.late = doc.change[which(doc.change$Time > 15),]
m1 = nls(formula = mm_form, data = doc.late)
summary(m1)
logLik(m1)
plot(doc.late$U[which(doc.late$Ctype == "Leaves")]~doc.late$PerDOCrem[which(doc.late$Ctype == "Leaves")], pch = 21,bg = "orange", xlim = c(13,40), ylim = c(0,0.15), xlab = "PerDOCrem", ylab = "U")
points(doc.late$U[which(doc.late$Ctype == "Pine")]~doc.late$PerDOCrem[which(doc.late$Ctype == "Pine")], pch = 21, bg = "purple")
points(doc.late$U[which(doc.late$Ctype == "Grasses")]~doc.late$PerDOCrem[which(doc.late$Ctype == "Grasses")], pch = 21, bg = "dark green")

Pine = doc.late[which(doc.late$Ctype == "Pine"),]
Leaves = doc.late[which(doc.late$Ctype == "Leaves"),]
Grasses = doc.late[which(doc.late$Ctype == "Grasses"),]

m_grasses = nls(U~SSasymp(PerDOCrem, Asym = 0.1, R0 = 1, lrc = -8), data = Grasses, trace = T);summary(m_pine)
getInitial(U~SSasymp(U, Asym = 0.2, R0 = 1, lrc = -6.2), data = Pine)

m_leaves = nls(mm_form, data = Leaves, trace = T);summary(m_leaves)
m_grasses = nls(mm_form, data = Grasses, trace = T);summary(m_grasses)

#######
Pine = doc.change[which(doc.change$Ctype == "Pine"),]
Leaves = doc.change[which(doc.change$Ctype == "Leaves"),]
Grasses = doc.change[which(doc.change$Ctype == "Grasses"),]

set.seed(123)
m_pine = nls(Ub~PerDOC^(b*PerDOCrem*Ncells), data = Pine)
lm_Pine = lm(Ub~PerDOCrem, data = Pine);summary(m_Pine)
set.seed(123)
m_Leaves = nls(formula = mm_form, data = Leaves[-3,], start = list( a = 0.5, b= 59), trace = T);summary(m_Leaves)
lm_Leaves = lm(Ub~PerDOCrem, data = Leaves);summary(lm_Leaves)
set.seed(123)
exp_Grasses = nls(Ub~a+b^PerDOCrem, data = Grasses, trace = T);summary(exp_Grasses)
lm_Grasses = lm(Ub~PerDOCrem, data = Grasses);summary(lm_Grasses)

PerDOCrem = seq(10,100, 1)
predict_leaves = predict(m_Leaves, newdata = data.frame(PerDOCrem = PerDOCrem))
predict_pine = predict(lm_Pine, newdata = data.frame(PerDOCrem = PerDOCrem))
predict_grasses = predict(exp_Grasses, newdata = data.frame(PerDOCrem = PerDOCrem))
leaves_df = data.frame(PerDOCrem, predict_leaves)
pine_df = data.frame(PerDOCrem, predict_pine)
grass_df = data.frame(PerDOCrem, predict_grasses)
#Plotting the data and model
plot(U~PerDOCrem, data = Leaves, pch = 21, col = "black", bg = "orange", cex = 2.5)
points(U~PerDOCrem, data = Pine, pch = 21, col = "black", bg = "purple", cex = 2.5)
points(U~PerDOCrem, data = Grasses, pch = 21, col = "black", bg = "dark green", cex = 2.5 )
lines(U~PerDOCrem, data = Leaves, lwd = 3, col = "orange")
lines(U~PerDOCrem, data = Pine, lwd = 3, col = "purple")
lines(U~PerDOCrem, data = Grasses, lwd = 3, col = "dark green")

lines(predict_leaves~PerDOCrem, data = leaves_df, col = "orange", lwd = 3)
lines(predict_pine~PerDOCrem, data = pine_df, col = "purple", lwd = 3)
lines(predict_grasses~PerDOCrem, data = grass_df, col = "dark green", lwd = 3)

ggplot(doc.change, aes(x = Time, y = PerDOCrem, group = Ctype)) + geom_point(aes(colour = Ctype), size = 4)
predict_leaves
#

plot(doc.change$absdDOC~doc.change$PerDOCrem)
lines(doc.change$PerDOCrem[which(doc.change$PerDOCrem != 100)], predict(m))
length(doc.change$PerDOCrem)
length(predict(m))
is.na(doc.change$absdDOC)
##### END MICHAELIS-MENTEN KINETICS PREDICTION ######

###Correlation  tables of DOC and other variables for each litter type
## Correlation matrix with p-values. See http://goo.gl/nahmV for documentation of this function

cor.prob <- function (X, dfr = nrow(X) - 2) {
  R <- cor(X, use="pairwise.complete.obs", method = "spearman")
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr/(1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  R[row(R) == col(R)] <- NA
  R
}

## Use this to dump the cor.prob output to a 4 column matrix
## with row/column indices, correlation, and p-value.
## See StackOverflow question: http://goo.gl/fCUcQ
flattenSquareMatrix <- function(m) {
  if( (class(m) != "matrix") | (nrow(m) != ncol(m))) stop("Must be a square matrix.") 
  if(!identical(rownames(m), colnames(m))) stop("Row and column names must be equal.")
  ut <- upper.tri(m)
  data.frame(i = rownames(m)[row(m)[ut]],
             j = rownames(m)[col(m)[ut]],
             cor=t(m)[ut],
             p=m[ut])
}

doctot_cor <- doc_tot[,c(2:21)]
#correlation matrix
cor(doctot_cor)
#correlation matrix with p-values
cor.prob(doctot_cor)
#'flatten' that table
flattenSquareMatrix(cor.prob(doctot_cor))
#plot the data
doctot.cor <- chart.Correlation(doctot_cor)

jpeg("DOC_cor_all.jpg", height = 6000, width = 6000, quality = 80, pointsize = 30)
chart.Correlation(doctot_cor, pch = 19)
dev.off()


###### End exploratory analysis ####
####Don't know what this stuff is below here. Clean later #####

Ncell.time.plot <- ggplot( data = doc, aes( x = PerDOCrem, y = log10(Ncells), group = Ctype, colour = Ctype)) +
  geom_line( size = 2.5)
Ncell.time.plot
##statistical test of slopes of logged data

doc.ln.lm <- lm(log.Per ~ log(time.mod) * Ctype, data = doc)
summary(doc.ln.lm)
##change in lability over time (EEMS)

eems.plot<-ggplot(data=doc, aes(x = Time, y = X.labile, group=Ctype, colour=Ctype)) +
	geom_line(size=2.5) + 
	geom_point (aes(fill = Ctype),size = 5.5, colour = "black", shape = 21) +
	ylab(expression(paste("% Lability"))) +
	ylim(c(25,90)) +
	xlab(expression(paste("Time (days)")))
eems.plot

#This is the final plot
dev.new()
F2 = eems.plot + scale_colour_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	scale_fill_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size = 2),
		axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
		axis.title.y = element_text(vjust = 1.2, size = 25), plot.margin = unit(rep(0.5, 4), "lines"), panel.border = element_blank(), 
		plot.background = element_rect(fill = "transparent", colour = NA), axis.line.y = element_line( size = 2), axis.text.y = element_text(size = 15, face = "bold"))





#Estimate the greatest change between days of each of the Ctype

library(quantmod)
resp.sum <- ddply(resp, c("Ctype", "Time"), summarize, CO2.mean = mean(CO2.accum, na.rm = T))
resp.se <- ddply(resp, c("Ctype", "Time"), summarize, CO2.se = 1.96 * (sd(CO2.accum, na.rm = T)/sqrt(3)))
resp.sum <- ddply(resp.sum, c("Ctype"), transform, CO2.del = round(Delt(CO2.mean),3))
head(resp.sum)
resp.sum$Delt.1.arithmetic[is.na(resp.sum$Delt.1.arithmetic)] <- 0
resp.sum$Delt.1.arithmetic[which(resp.sum$Delt.1.arithmetic == "Inf")] <- 0
del.max = ddply(resp.sum, c("Ctype"), summarize, max(Delt.1.arithmetic))

## 
## Estimate when (time) Ctype differ in CO2 accum and when/if (time) they converge

## Pair-wise function

resp = read.csv(file.choose(),T) #only need to run if you have not yet

group = unique(resp$Time)

pairwise.Ttest <- function(x, group, p.adjust.m){
  co = as.matrix(combn(unique(x$Ctype),2))
  pairs = c()
  p.value = c()
  for(i in group){
    Xi = x[which(x$Time == i),];
    pw.t = pairwise.t.test(Xi$CO2.accum, Xi$Ctype, p.adjust.method = "none"); #runs the pairwise t.test
    p.value = cbind(p.value, c(pw.t$p.value[[1]], pw.t$p.value[[2]], pw.t$p.value[[4]])) 
  }
  p.adjusted = matrix(p.adjust(p.value, method = p.adjust.m), nrow = 3)
  #pairs = c(pairs, paste(as.character(co[1,1]), 'vs', as.character(co[2,1])), paste(as.character(co[1,2]), 'vs', as.character(co[2,2])), paste(as.character(co[1,3]), 'vs', as.character(co[2,3])))
  pairs = c(paste(as.character(co[1,1]), 'vs', as.character(co[2,1])), paste(as.character(co[1,2]), 'vs', as.character(co[2,2])), paste(as.character(co[1,3]), 'vs', as.character(co[2,3])))
  pairwise = data.frame(pairs, p.adjusted)
  #pairwise = data.frame(pairs, p.value)
  
  colnames(pairwise) = c("Pairs", paste(as.character(unique(x$Time))))
  print(paste("P value adjustment method:", p.adjust.m))
  return(pairwise)
  write.table(pairwise, file = "Resp_time_Ttest.txt", sep = "\t", quote = F, append = F)
}
pairwise.Ttest(resp, group, p.adjust.m = "bonf")
pairwise.Ttest(resp, group, p.adjust.m = "holm")
##########################
#######################################
cor.test(doc_gr2$X.labile, doc_gr2$doc.cell)


### Build cell count over time to 

cell.plot <- cell.plot +
  geom_line(data = doc, aes( x = Time, y = log10(Ncells/10000000), group = Ctype, colour = Ctype), size = 1.5)

cell.plot <- cell.plot +
  geom_line(data = doc, aes( x = Time, y = log10(DOC), group = Ctype, colour = Ctype), linetype = "dashed", size = 1.5)
cell.plot

####Test plot for groups
cell.plot5 <- ggplot( data = doc, aes( x = X.labile, y = doc.cell)) +
	geom_point( size = 5.5) +
	stat_smooth(data = doc_gr2, method = "lm", se = F) +	
	#scale_colour_gradient(low = "red") +
	xlab(expression('%'~Lability[EEMs]*'')) +
	ylim(c(-0.5,10)) +
	ylab(expression(atop("DOC Processing efficiency", paste("("*Delta*"DOC/"*10^6~"cells/day)"))))
cell.plot5




fluor.plot <- ggplot(data = doc, aes( x = Time, y = fluoride, group = Ctype, colour = Ctype, shape = Ctype)) +
	geom_line ( size = 2.5) +
	geom_point (aes(fill = Ctype), size = 5.5, colour = "black", shape = 21, na.rm = T) +
	geom_errorbar ( aes(ymin = fluoride - fluor.se, ymax = fluoride + fluor.se)) +
	ylab(bquote('Fluoride (mg/L)')) 
fluor.plot

fluor.plot + scale_colour_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	scale_fill_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	theme(legend.position = c(0.93, 0.9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
		axis.title.y = element_text(vjust = 2), panel.border = element_blank(), axis.line.y = element_line())

Chlor.plot <- ggplot(data = doc, aes( x = Time, y = chloride, group = Ctype, colour = Ctype)) +
	geom_line ( size = 2.5) +
	geom_point (aes(fill = Ctype), size = 5.5, colour = "black", shape = 21, na.rm = T) +
	geom_errorbar ( aes(ymin = chloride - cl.se, ymax = chloride + cl.se)) +
	ylab(bquote('Chloride (mg/L)')) 

Chlor.plot + scale_colour_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	scale_fill_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	theme(legend.position = c(0.93, 0.9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
		axis.title.y = element_text(vjust = 2), panel.border = element_blank(), axis.line.y = element_line())

nit.plot <- ggplot(data = doc, aes( x = Time, y = N, group = Ctype, colour = Ctype)) +
	geom_line ( size = 2.5) +
	geom_point (aes(fill = Ctype), size = 5.5, colour = "black", shape = 21, na.rm = T) +
	geom_errorbar ( aes(ymin = N - (no2.se + no3.se), ymax = N + (no2.se + no3.se))) +
	ylab(bquote('Nitrite + Nitrate (mg/L)')) 

F4 = nit.plot + scale_colour_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	scale_fill_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line = element_line(size = 2), 
		axis.line.y = element_line(), axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
		axis.title.x = element_blank(), axis.title.y = element_text(vjust = 1, size = 25),  plot.margin = unit(rep(0.5, 4), "lines"),
		plot.background = element_rect(fill = "transparent", colour = NA), axis.text.y = element_text(size = 15, face = "bold"))

Phos.plot <- ggplot(data  = doc, aes( x = Time, y = phosphate, group = Ctype, colour = Ctype)) +
	geom_line ( size = 2.5) +
	geom_point (aes(fill = Ctype), size = 5.5, colour = "black", shape = 21, na.rm = T) +
	geom_errorbar ( aes(ymin = phosphate - po4.se, ymax = phosphate + po4.se)) +
	ylab(bquote('Phosphate (mg/L)')) +
	xlab(expression(paste("Time (days)"))) 

F5 = Phos.plot + scale_colour_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	scale_fill_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		axis.line = element_line(size = 2), axis.text.x = element_text(size = 15, face = "bold"), axis.ticks.x = element_line(), axis.title.x = element_text(size = 20, face = "bold"),
		axis.title.y = element_text(vjust = 1, size = 25),  plot.margin = unit(rep(0.5, 4), "lines"), panel.border = element_blank(), axis.line.y = element_line(),
		plot.background = element_rect(fill = "transparent", colour = NA), axis.text.y = element_text(size = 15, face = "bold"))


plots <- list(F1, F2, F3, F4, F5)
grobs <- list()
widths <- list()

for (i in 1:length(plots)){
    grobs[[i]] <- ggplotGrob(plots[[i]])
    widths[[i]] <- grobs[[i]]$widths[2:5]
}
dev.off()

pdf("final_plot.pdf", height = 20, width = 15, bg = "white")


maxwidth <- do.call(grid::unit.pmax, widths)

for (i in 1:length(grobs)){
     grobs[[i]]$widths[2:5] <- as.list(maxwidth)
}

do.call("grid.arrange", c(grobs, ncol = 1))

dev.off()


##Leachate cumulative and quality measures

doc_tot <- read.csv(file="Leachate_Cumulative.csv",T)

CN.plot <- ggplot(data  = doc_tot, aes( x = CN, y = TotCRes, group = Ctype, colour = Ctype)) +
	geom_point (aes(fill = Ctype), size = 5.5, colour = "black", shape = 21, na.rm = T) +
	ylab(bquote('% DOC loss')) +
	xlab(expression(paste("C:N (molar)")))
CN.plot 

 

F6 = CN.plot + scale_colour_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	scale_fill_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	theme(legend.position = c(0.20, 0.90), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		axis.line = element_line(size = 2), axis.text.x = element_text(size = 15, face = "bold"), axis.ticks.x = element_line(), axis.title.x = element_text(size = 20, face = "bold"),
		axis.title.y = element_text(vjust = 1, size = 25),  plot.margin = unit(rep(0.5, 4), "lines"), panel.border = element_rect(), axis.line.y = element_line(),
		plot.background = element_rect(fill = "transparent", colour = NA), axis.text.y = element_text(size = 15, face = "bold"))


CP.plot <- ggplot(data  = doc_tot, aes( x = CP, y = TotCRes, group = Ctype, colour = Ctype)) +
	geom_point (aes(fill = Ctype), size = 5.5, colour = "black", shape = 21, na.rm = T) +
	ylab(bquote('% DOC loss')) +
	xlab(expression(paste("C:P (molar)")))
CP.plot 

F7 = CP.plot + scale_colour_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	scale_fill_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		axis.line = element_line(size = 2), axis.text.x = element_text(size = 15, face = "bold"), axis.ticks.x = element_line(), axis.title.x = element_text(size = 20, face = "bold"),
		axis.title.y = element_text(vjust = 1, size = 25),  plot.margin = unit(rep(0.5, 4), "lines"), panel.border = element_rect(), axis.line.y = element_line(),
		plot.background = element_rect(fill = "transparent", colour = NA), axis.text.y = element_text(size = 15, face = "bold"))

NP.plot <- ggplot(data  = doc_tot, aes( x = NP, y = TotCRes, group = Ctype, colour = Ctype)) +
	geom_point (aes(fill = Ctype), size = 5.5, colour = "black", shape = 21, na.rm = T) +
	ylab(bquote('% DOC loss')) +
	xlab(expression(paste("N:P (molar)"))) 

F8 = NP.plot + scale_colour_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	scale_fill_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		axis.line = element_line(size = 2), axis.text.x = element_text(size = 15, face = "bold"), axis.ticks.x = element_line(), axis.title.x = element_text(size = 20, face = "bold"),
		axis.title.y = element_text(vjust = 1, size = 25),  plot.margin = unit(rep(0.5, 4), "lines"), panel.border = element_rect(), axis.line.y = element_line(),
		plot.background = element_rect(fill = "transparent", colour = NA), axis.text.y = element_text(size = 15, face = "bold"))

MLB.plot <- ggplot(data  = doc_tot, aes( x = MLB, y = TotCRes, group = Ctype, colour = Ctype)) +
	geom_point (aes(fill = Ctype), size = 5.5, colour = "black", shape = 21, na.rm = T) +
	ylab(bquote('% DOC loss')) +
	xlab(expression(paste("MLB"))) 

F9 = MLB.plot + scale_colour_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	scale_fill_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		axis.line = element_line(size = 2), axis.text.x = element_text(size = 15, face = "bold"), axis.ticks.x = element_line(), axis.title.x = element_text(size = 20, face = "bold"),
		axis.title.y = element_text(vjust = 1, size = 25),  plot.margin = unit(rep(0.5, 4), "lines"), panel.border = element_rect(), axis.line.y = element_line(),
		plot.background = element_rect(fill = "transparent", colour = NA), axis.text.y = element_text(size = 15, face = "bold"))

MLBw.plot <- ggplot(data  = doc_tot, aes( x = MLBw, y = TotCRes, group = Ctype, colour = Ctype)) +
	geom_point (aes(fill = Ctype), size = 5.5, colour = "black", shape = 21, na.rm = T) +
	ylab(bquote('% DOC loss')) +
	xlab(expression(paste("MLBw"))) 

F10 = MLBw.plot + scale_colour_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	scale_fill_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		axis.line = element_line(size = 2), axis.text.x = element_text(size = 15, face = "bold"), axis.ticks.x = element_line(), axis.title.x = element_text(size = 20, face = "bold"),
		axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(rep(0.5, 4), "lines"), panel.border = element_rect(), axis.line.y = element_blank(),
		plot.background = element_rect(fill = "transparent", colour = NA), axis.text.y = element_blank())


EEMS.plot <- ggplot(data  = doc_tot, aes( x = EEMS, y = TotCRes, group = Ctype, colour = Ctype)) +
	geom_point (aes(fill = Ctype), size = 5.5, colour = "black", shape = 21, na.rm = T) +
	ylab(bquote('% DOC loss')) +
	xlab(expression(paste("EEMS"))) 

F11 = EEMS.plot + scale_colour_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	scale_fill_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		axis.line = element_line(size = 2), axis.text.x = element_text(size = 15, face = "bold"), axis.ticks.x = element_line(), axis.title.x = element_text(size = 20, face = "bold"),
		axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(rep(0.5, 4), "lines"), panel.border = element_rect(), axis.line.y = element_blank(),
		plot.background = element_rect(fill = "transparent", colour = NA), axis.text.y = element_blank())

pdf("lability_plot.pdf", width = 15)
grid.arrange(F6, F11, F10, ncol = 3, nrow = 1)
dev.off()

#############

stoic.plot2 <- ggplot(data  = doc_tot, aes( x = CN, y = co2accum, group = Ctype, colour = Ctype)) +
	geom_point (aes(fill = Ctype), size = 5.5, colour = "black", pch = 21, na.rm = T) +
	geom_point(data = doc_tot, aes(x = CP, y = co2accum, group = Ctype, colour = Ctype, fill = Ctype), size = 5.5, pch = 24, colour = "black") +
	ylim(c(600,800)) +
	ylab(bquote(~CO[2]~ 'Accumulation (' *mu* 'L)')) +
	scale_x_log10() +
	xlab(expression(paste("C:X (molar)"))) 
stoic.plot2


F12 = stoic.plot2 + scale_colour_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	scale_fill_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		axis.line = element_line(size = 2), axis.text.x = element_text(size = 15, face = "bold"), axis.ticks.x = element_line(), axis.title.x = element_text(size = 20, face = "bold"),
		axis.title.y = element_text(vjust = 1, size = 25),  plot.margin = unit(rep(0.5, 4), "lines"), panel.border = element_rect(), axis.line.y = element_line(),
		plot.background = element_rect(fill = "transparent", colour = NA), axis.text.y = element_text(size = 15, face = "bold"))
F12


EEMS.plot2 <- ggplot(data  = doc_tot, aes( x = EEMS, y = co2accum, group = Ctype, colour = Ctype)) +
	geom_point (aes(fill = Ctype), size = 5.5, colour = "black", shape = 21, na.rm = T) +
	ylim(c(600,800)) +
	ylab(bquote(~CO[2]~ 'Accumulation (' ~mu* 'g)')) +
	xlab(expression(paste("EEMS")))
EEMS.plot2
 

F13 = EEMS.plot2 + scale_colour_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	scale_fill_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		axis.line = element_line(size = 2), axis.text.x = element_text(size = 15, face = "bold"), axis.ticks.x = element_line(), axis.title.x = element_text(size = 20, face = "bold"),
		axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(rep(0.5, 4), "lines"), panel.border = element_rect(), axis.line.y = element_blank(),
		plot.background = element_rect(fill = "transparent", colour = NA), axis.text.y = element_blank())


MLBw.plot2 <- ggplot(data  = doc_tot, aes( x = MLBw, y = co2accum, group = Ctype, colour = Ctype)) +
	geom_point (aes(fill = Ctype), size = 5.5, colour = "black", shape = 21, na.rm = T) +
	ylim(c(600,800)) +
	ylab(bquote(~CO[2]~ 'Accumulation (' ~mu* 'g)')) +
	xlab(expression(paste("MLBw"))) 

F14 = MLBw.plot2 + scale_colour_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	scale_fill_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = cbbPalette) +
	theme(legend.position = c(0.85,0.9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		axis.line = element_line(size = 2), axis.text.x = element_text(size = 15, face = "bold"), axis.ticks.x = element_line(), axis.title.x = element_text(size = 20, face = "bold"),
		axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(rep(0.5, 4), "lines"), panel.border = element_rect(), axis.line.y = element_blank(),
		plot.background = element_rect(fill = "transparent", colour = NA), axis.text.y = element_blank())

pdf("lability_plot2.pdf", width = 15, height = 7)
grid.arrange(F12, F13, F14, ncol = 3, nrow = 1)
dev.off()

######
nut.plot4 <- ggplot(data = doc, aes( x = Time, y = phosphate, group = Ctype, colour = Ctype)) +
	geom_line( size = 1.5)
nut.plot4

CN.plot <- ggplot(data = doc, aes (x = PerDOCrem, y = nitrate, group = Ctype, colour = Ctype)) +
	geom_line( size = 1.5)
CN.plot

NP.plot <- ggplot( data = doc, aes ( x = nitrate, y = phosphate, group = Ctype, colour = Ctype)) +
	geom_line( size = 1.5)
NP.plot


##Linear models of the decay of each litter type

lm.grass <- lm(log.Per ~ log(time.mod), data = grass); summary(lm.grass) 
lm.pine <- lm(log.Per ~ log(time.mod), data = pine); summary(lm.pine) 
lm.leaves <- lm(log.Per ~ log(time.mod), data = leaves); summary(lm.leaves) 





#This is the final respiration plot (greyscale)

F_resp = ggplot(data = resp, aes(x = Time, y = CO2.accum, group = Ctype, colour = Ctype, shape = Ctype)) +
	geom_line(size = 2) +	
	geom_errorbar(aes(ymin = CO2.accum - stdev.CO2.accum, ymax = CO2.accum + stdev.CO2.accum), colour = "black") +
	geom_point(aes(fill = Ctype, shape = Ctype), colour = "black", size = 3, stroke = 1) + 
	scale_colour_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = c("grey0", "grey", "grey40")) +
	scale_fill_manual(name="Litter type", breaks=c("Pine", "Leaves", "Grasses"),
	labels = c("Pine", "Leaves", "Grasses"), values = c("grey0", "grey", "grey40")) +
	scale_shape_manual(name = "Litter type", breaks=c("Pine", "Leaves", "Grasses"), 
	labels = c("Pine", "Leaves", "Grasses"), values = c(21,22,25)) +
	ylab(bquote(~CO[2]~ 'accumulation ('*mu*'L)')) +
	xlab("Time (days)") +
	scale_x_continuous(limits = c(0,50), expand = c(0.01, 0)) +
	#scale_y_continuous(limits = c(5,9), expand = c(0,0)) +
	theme(axis.line = element_line(size = 2),	panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_text(vjust = 1.2, size = 25), plot.margin = unit(rep(0.5, 4), "lines"), panel.border = element_blank(), 
		plot.background = element_rect(fill = "transparent", colour = NA), axis.line.y = element_line( size = 2), axis.text.y = element_text(size = 22, face = "bold"), axis.text.x = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(size = 24), legend.justification = c(1,0), legend.position = c(1,0))
F_resp

tiff("resp_plot.tiff", height = 10, width = 15, units = "in", compression = "lzw", res = 300)
F_resp
dev.off()


#################  Old code #####################
#old cell respiration efficiency plot code
cell.plot1 <- ggplot( data = doc, aes( x = labile, y = doc.cell, group = Time, colour = Time)) +
  geom_point( size = 5.5) +
  scale_colour_gradient(low = "red") +
  xlab(bquote('%'~Lability[T-1]*'')) +
  ylab(bquote('Change in [DOC]'));cell.plot1

cell.plot2 <- ggplot( data = doc, aes( x = X.labile, y = doc.cell, group = Ncells, colour = Ncells)) +
  geom_point( size = 5.5) +
  scale_x_continuous() +
  scale_colour_gradient(low = "red") +
  xlab(expression('%'~Lability[T]*'')) +
  ylab(bquote('Change in [DOC]'));cell.plot2

pdf("cell_plot.pdf")
grid.arrange(cell.plot1, cell.plot2, ncol = 1, nrow = 2)
dev.off()

cell.plot3 <- ggplot( data = doc, aes( x = labile, y = doc.cell, group = Ncells, colour = Ncells)) +
  geom_point( size = 5.5) +
  scale_colour_gradient(low = "red") +
  xlab(bquote('%'~Lability[T-1]*'')) +
  ylab(bquote('Change in [DOC]'));cell.plot3

set.seed(101)
Feff <- ggplot( data = doc, aes( x = X.labile, y = doc.cell)) +
  stat_smooth(data = doc_gr2, method = "lm", se = T, alpha = 0.6, size = 1.5, colour = "black", linetype = 1) +
  geom_point( size = 5.5, aes(shape = Ctype, colour = Ncells)) +
  scale_shape_manual(name="Leachate type", breaks=c("Grasses", "Leaves", "Pine"),
                     labels = c("Grasses", "Leaves", "Pine"), values = c(16,17,15)) +
  scale_colour_gradient(low = "blue", high = 'red', name = "Cell count") +
  xlab(expression('%'~Lability[EEMs]*'')) +
  ylim(c(-10,20)) + coord_cartesian(ylim = c(-0.5,10))+
  scale_x_reverse() +
  geom_rug(sides = 't')+
  #geom_vline(xintercept = 75, linetype = "dashed", size = 1.2) +
  #geom_text(data = doc, aes(x = X.labile, y = 9.8, label = time.mod), size = 3) +
  geom_text_repel(data = doc, aes(x = X.labile, y = 10, label = Time), size = 3, ylim = c(NA, 10.2), direction = "y", point.padding = NA, segment.colour = NA, force = 1.2) +
  ylab(expression(atop("DOM Processing Efficiency", paste("(-"*Delta*"[DOC]/"*10^6~"cells/day)")))) +
  theme(legend.position = c(0.1, 0.65), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.text.y = element_text(size = 22), axis.text.x = element_text(size = 22), 
        axis.title.x = element_text(size = 24), axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 0, l = 0, b = 0)));Feff

tiff("Figure3.tiff", width = 13, height = 10, units = "in", res = 300)
Feff
dev.off()
###############Breakpoint linear regression code#################
full <- function(mainplot = NULL, subplot = NULL, vp = NULL) {
  print(mainplot)
  theme_set(theme_bw(12))
  print(subplot, vp = vp)
  theme_set(theme_bw())
}

full()

library(segmented)

m.seg=segmented(m, seg.Z=~Time, psi=55)
m.int.seg=segmented(m, seg.Z=~Time, psi=55)

##Code for one breakpoint with the SR data only
m.SR=lm(Oxygen ~ Time, data=Resp_R[Resp_R$Carbon=="SR",])
Resp_R[Resp_R$Carbon=="SR",]
m.SR.seg=segmented(m.SR, seg.Z=~Time, psi=55)
summary(m.SR.seg)

#####Don't know if I need this code to pull in the specific dropbox
###directory for each user
get.dropbox.folder <- function() {
  
  if (!require(RCurl)) stop ("You need to install RCurl package.")
  if (Sys.info()["sysname"]!="Windows") stop("Currently, 'get.dropbox.folder' works for Windows only. Sorry.")
  
  db.file <- paste(Sys.getenv('APPDATA'), '\\Dropbox\\host.db', sep='')
  base64coded <- readLines(db.file, warn=FALSE)[2]
  
  base64(base64coded, encode=FALSE)



#####This code was for the initial ggplot but kept crashing

#doc.plot<-ggplot(data=doc, aes(x=Time, y=Per, group=Ctype, colour=Ctype)) +
	#geom_line(size=1.5) + 
	#ylab(expression(paste("DOC remaining (%)"))) +
	#xlab(expression(paste("Time (days)"))) +
	#xlim(0,45) +
	#geom_point(size=3, fill="white") +
	#scale_shape_manual(values=c(22,21,23)) +
	#theme(legend.position="top") +
	#theme(axis.title.x = element_text(face="bold", size=20),
	#axis.text.x = element_text(colour = "black", size=20),
	#axis.title.y = element_text(face="bold", size=20),
	#axis.text.y = element_text(colour = "black", size=20),
	#panel.grid.major = element_blank(),
	#panel.grid.minor = element_blank(),
	#panel.background = element_blank(),
	#panel.border = element_rect(colour = "black", fill=NA),
	#legend.background = element_rect(fill = FALSE),
	#legend.text = element_text(size = 20),
	#legend.text = element_text(size = 20))
	#legend.key = element_rect(colour = "Black", size = 1.0))
#doc.plot
}
``
grid.draw(cbind_gtable_max(F1, F2, F3, F4))

#####for making figures #####
rbind_gtable_max <- function(...){
  
  gtl <- list(...)
  stopifnot(all(sapply(gtl, is.gtable)))
  bind2 <- function (x, y) 
  {
    stopifnot(ncol(x) == ncol(y))
    if (nrow(x) == 0) 
      return(y)
    if (nrow(y) == 0) 
      return(x)
    y$layout$t <- y$layout$t + nrow(x)
    y$layout$b <- y$layout$b + nrow(x)
    x$layout <- rbind(x$layout, y$layout)
    x$heights <- gtable:::insert.unit(x$heights, y$heights)
    x$rownames <- c(x$rownames, y$rownames)
    x$widths <- grid::unit.pmax(x$widths, y$widths)
    x$grobs <- append(x$grobs, y$grobs)
    x
  }
  
  Reduce(bind2, gtl)
}

cbind_gtable_max <- function(...){
  
  gtl <- list(...)
  stopifnot(all(sapply(gtl, is.gtable)))
  bind2 <- function (x, y) 
  {
    stopifnot(nrow(x) == nrow(y))
    if (ncol(x) == 0) 
      return(y)
    if (ncol(y) == 0) 
      return(x)
    y$layout$l <- y$layout$l + ncol(x)
    y$layout$r <- y$layout$r + ncol(x)
    x$layout <- rbind(x$layout, y$layout)
    x$widths <- gtable:::insert.unit(x$widths, y$widths)
    x$colnames <- c(x$colnames, y$colnames)
    x$heights <- grid::unit.pmax(x$heights, y$heights)
    x$grobs <- append(x$grobs, y$grobs)
    x
  }
  Reduce(bind2, gtl)
}
#########
###Fucking around with O2-CO2 stoichiometry###
colnames(resp)[2] = "Days"
o2_resp %>%
  gather(Ctype, O2_consump, Grasses:Pine, factor_key = T) %>%
  group_by(Ctype) %>%
  left_join(resp) %>%
  mutate(O_C = abs(O2_consump)/abs(CO2.accum)*(48/36)) -> resp_o2.df

ggplot(resp_o2.df, aes(x = Days, y = O_C, color = Ctype)) + geom_point(size = 2) + geom_path(size = 1.2) +
  ylab("O2 consumption (uL)/CO2 accum. (uL)")+theme(legend.position = c(0.8,0.2))


o2_resp %>%
  gather(Ctype, O2_consump, Grasses:Pine, factor_key = T) %>%
  group_by(Ctype) %>%
  mutate(dO2 = c(NA, diff(O2_consump))) -> o2.change

co2.change <- resp2 %>% group_by(Ctype,Rep) %>% mutate(dCO2 = c(NA, diff(CO2.accum)))
co2.change %>%
  group_by(Ctype, Time) %>%
  summarize(dCO2 = mean(dCO2, na.rm = T)) -> co2.change
colnames(co2.change)[2] = "Days"
co2.change %>%
  left_join(o2.change, by = c("Ctype","Days")) %>%
  mutate(O_C = abs(dO2)/abs(dCO2)) -> O_C.change

ggplot(O_C.change, aes(x = Days, y = log(O_C), color = Ctype)) + geom_point(size = 2) + geom_path(size = 1.2) +
theme(legend.position = c(0.8,0.2))

