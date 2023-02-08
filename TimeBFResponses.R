# Time of BF responses
# 
library(purrr)
library(tidyr)
library(dplyr)
library(readxl)
library(aod)

library(reshape2)
library(data.table)
library(ggplot2)
library(ggpubr)

library(stats)
library(xlsx)
library(rstatix)

library(egg)

BFData <- read_excel(file.choose()) # use Counting Cells_ResponseTimes

# separated into different kinds of times (cell response, end plasmolysis, end deplasmolysis)
responsetimes <- na.omit(BFData[,c("Time Response (s)", "CPA")])
responsetimes$CPA <- paste0("a.resp", responsetimes$CPA)
plastimes <- na.omit(BFData[,c("End Plasmolysis", "CPA")])
plastimes$CPA <- paste0("b.plas", plastimes$CPA)
deplastimes <- na.omit(BFData[,c("End Deplasmolysis", "CPA")])
deplastimes$CPA <- paste0("c.deplas", deplastimes$CPA)


# make boxplot DF
timeBP <- list(responsetimes, plastimes, deplastimes)
timeBP <- dplyr::bind_rows(timeBP)
timeBP <- na.omit(reshape2::melt(timeBP))
timeBP <- timeBP[,c("CPA","value")]


splitBP <- split(timeBP, timeBP$CPA)

deplasDMSO <- splitBP$c.deplasDMSO
deplasEG <- splitBP$c.deplasEG
deplasGlyc <- splitBP$c.deplasGlyc

plasDMSO <- splitBP$b.plasDMSO
plasEG <- splitBP$b.plasEG
plasGlyc <- splitBP$b.plasGlyc

respDMSO <- splitBP$a.respDMSO
respEG <- splitBP$a.respEG
respGlyc <- splitBP$a.respGlyc

# DMSO
DMSOtimes <- list(deplasDMSO, plasDMSO, respDMSO)
DMSOtimes <- dplyr::bind_rows(DMSOtimes)
DMSOtimes <- na.omit(reshape2::melt(DMSOtimes))
DMSOtimes <- DMSOtimes[,c("CPA","value")]

EGtimes <- list(deplasEG, plasEG, respEG)
EGtimes <- dplyr::bind_rows(EGtimes)
EGtimes <- na.omit(reshape2::melt(EGtimes))
EGtimes <- EGtimes[,c("CPA","value")]

glyctimes <- list(deplasGlyc, plasGlyc, respGlyc)
glyctimes <- dplyr::bind_rows(glyctimes)
glyctimes <- na.omit(reshape2::melt(glyctimes))
glyctimes <- glyctimes[,c("CPA","value")]


# plots
t1 <- theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_line(color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        axis.title=element_text(size=14,face="bold"))

allTimesBP <- ggplot(timeBP, aes(y = value, x = CPA)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(color = "black", aes(fill = CPA)) + t1 +
  theme(legend.position="bottom")

allTimes_BoxplotStats <- boxplot(value ~ CPA, data = timeBP)
allTimes_BPSummary <- allTimes_BoxplotStats$stats
new_row<- t(as.data.frame(allTimes_BoxplotStats$n))
allTimes_BPSummary <- rbind(new_row, allTimes_BPSummary)
colnames(allTimes_BPSummary) <- allTimes_BoxplotStats$names
rownames(allTimes_BPSummary)<-c("n", "Min","First Quartile","Median","Third Quartile","Maximum")

#Save your own file!
mypath <- file.path("D:", "Writing", "Cell response Paper", "Files for Sharing", "Bright Field Data", fsep="/",
                    paste("AllBFresponseTimes.csv"))
write.csv(allTimes_BPSummary, mypath)


DMSOtimesBP <- ggplot(DMSOtimes, aes(y = value, x = CPA)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(color = "black", aes(fill = CPA)) + t1 +
  theme(legend.position="bottom")

EGtimesBP <- ggplot(EGtimes, aes(y = value, x = CPA)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(color = "black", aes(fill = CPA)) + t1 +
  theme(legend.position="bottom")

glyctimesBP <- ggplot(glyctimes, aes(y = value, x = CPA)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(color = "black", aes(fill = CPA)) + t1 +
  theme(legend.position="bottom")

deplastimesBP <- ggplot(deplastimes, aes(y = `End Deplasmolysis`, x = CPA)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(color = "black", aes(fill = CPA)) + t1 +
  theme(legend.position="bottom")

plastimesBP <- ggplot(plastimes, aes(y = `End Plasmolysis`, x = CPA)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(color = "black", aes(fill = CPA)) + t1 +
  theme(legend.position="bottom")

resptimesBP <- ggplot(responsetimes, aes(y = `Time Response (s)`, x = CPA)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(color = "black", aes(fill = CPA)) + t1 +
  theme(legend.position="bottom")

# Statistical Testing ####
# Shapiro Test for normality
shapiro.test(plasDMSO$value) #not normal p<0.0001
shapiro.test(plasEG$value) #not normal p<0.0001
shapiro.test(plasGlyc$value) #not normal p<0.0001

shapiro.test(deplasDMSO$value) #not normal p<0.0001
shapiro.test(deplasEG$value) #normal p=0.08
shapiro.test(deplasGlyc$value) #normal p=0.6

shapiro.test(respDMSO$value) #not normal p<0.0001
shapiro.test(respEG$value) #not normal p<0.0001
shapiro.test(respGlyc$value) #not normal p<0.001

# use non-parametric testing
KWtest_BF <- kruskal.test(value ~ CPA, timeBP) 
KWtest_BF # p<0.0001
# Dunn test on all data ####
dunn.all <- timeBP %>% dunn_test(value ~ CPA, p.adjust.method = "BY")
dunn.all.DMSO <- dunn.all[c(3,6,24),]
dunn.all.EG <- dunn.all[c(11,14,29),]
dunn.all.glyc <- dunn.all[c(18,21,33),]
dunn.all.resp <- dunn.all[c(1,2,9),]
dunn.all.plas <- dunn.all[c(22,23,27),]
dunn.all.deplas <- dunn.all[c(34,35,36),]

mypath <- file.path("D:", "Writing", "Cell response Paper", "Files for Sharing", "Bright Field Data", fsep="/",
                    paste("DunnAllBrightfield.csv"))
write.csv(dunn.all, mypath)
# Dunn test with BY adjustment since these are dependent on one another
dunn.DMSO <- DMSOtimes %>% dunn_test(value ~ CPA, p.adjust.method = "BY") 
dunn.EG <- EGtimes %>% dunn_test(value ~ CPA, p.adjust.method = "BY") 
dunn.glyc <- glyctimes %>% dunn_test(value ~ CPA, p.adjust.method = "BY") 

DMSO_stats <- DMSOtimesBP +
  scale_y_continuous(breaks=seq(0,2000,100)) + 
  labs(y = "Time (s)", fill = "Action:") +
  scale_fill_manual(labels = c("Response", "Plasmolysis", "Deplasmolysis"),
                    values=c("#FFA5AB", "#AD646F","#5B2333")) +  
  stat_pvalue_manual(dunn.all.DMSO, y.position = 710, step.increase = 0.05, size = 5, tip.length = 0.01) +
  expand_limits(y = 1310) +
  expand_limits(y = 0) +
  theme(legend.title = element_text(size=15, 
                                    face="bold")) +
  theme(legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE))

EG_stats <- EGtimesBP +
  scale_y_continuous(breaks=seq(0,2000,100)) + 
  labs(y = "Time (s)") +
  scale_fill_manual(labels = c("Response", "Plasmolysis", "Deplasmolysis"),
                    values=c("#77C3BC", "#388697", "#21525D")) +  
  stat_pvalue_manual(dunn.all.EG, y.position = 540,step.increase = 0.05, size = 5, tip.length = 0.01) +
  expand_limits(y = 1310) +
  expand_limits(y = 0) +
  theme(axis.title.y = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE))

glyc_stats <- glyctimesBP +
  scale_y_continuous(breaks=seq(0,2000,100)) + 
  labs(y = "Time (s)") +
  scale_fill_manual(labels = c("Response", "Plasmolysis", "Deplasmolysis"),
                    values=c("#E29955", "#9C5716", "#844000")) +  
  stat_pvalue_manual(dunn.all.glyc, y.position = 1150,step.increase = 0.05, size = 5, tip.length = 0.01) +
  theme(axis.title.y = element_blank()) +
  expand_limits(y = 1310) +
  expand_limits(y = 0) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE))

deplas_stats <- deplastimesBP +
  scale_y_continuous(breaks=seq(0,2000,100)) + 
  labs(y = "Time (s)") +
  scale_fill_manual(labels = c("DMSO", "Ethyelene Glycol", "Glycerol"),
                    values=c("#5B2333", "#21525D", "#844000")) +  
  stat_pvalue_manual(dunn.all.deplas, y.position = 900,step.increase = 0.05, size = 5, tip.length = 0.01) +
  theme(axis.title.y = element_blank()) +
  expand_limits(y = 1310) +
  expand_limits(y = 0) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE))

plas_stats <- plastimesBP +
  scale_y_continuous(breaks=seq(0,2000,100)) + 
  labs(y = "Time (s)") +
  scale_fill_manual(labels = c("DMSO", "Ethyelene Glycol", "Glycerol"),
                    values=c("#AD646F", "#388697", "#9C5716")) +  
  stat_pvalue_manual(dunn.all.plas, y.position = 1150,step.increase = 0.05, size = 5, tip.length = 0.01) +
  theme(axis.title.y = element_blank()) +
  expand_limits(y = 1310) +
  expand_limits(y = 0) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE))

resp_stats <- resptimesBP +
  scale_y_continuous(breaks=seq(0,2000,100)) + 
  labs(y = "Time (s)", fill = "CPA:") +
  scale_fill_manual(labels = c("DMSO", "Ethyelene Glycol", "Glycerol"),
                    values=c("#FFA5AB", "#77C3BC", "#E29955")) +  
  stat_pvalue_manual(dunn.all.resp, y.position = 900,step.increase = 0.05, size = 5, tip.length = 0.01) +
  expand_limits(y = 1310) +
  expand_limits(y = 0) +
  theme(legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE))


#Final Plots ####
png(file = file.path("D:", "Writing", "Cell response Paper", "Files for Sharing", "Bright Field Data", fsep="/",
                     paste("BF_TimeBoxPlot.png")), units="in", width=10, height=10, res=300)
BF_TimeBoxPlot <- ggarrange(DMSO_stats, EG_stats, glyc_stats, 
                              labels = c("a.", "b.", "c."),
                              ncol = 3, nrow = 1) 
dev.off()

png(file = file.path("D:", "Writing", "Cell response Paper", "Files for Sharing", "Bright Field Data", fsep="/",
                     paste("BF_ResponseType.png")), units="in", width=10, height=10, res=300)
BF_TimeBoxPlot_response <- ggarrange(resp_stats, plas_stats, deplas_stats, 
                            labels = c("a.", "b.", "c."),
                            ncol = 3, nrow = 1) 
dev.off()

# Timing of CARS responses ####
CARSData <- read_excel(file.choose())

plastimes_CARS <- na.omit(CARSData[,c("Time Plas Start", "CPA")])
plastimes_CARS$CPA <- paste0("a.plas", plastimes_CARS$CPA)
names(plastimes_CARS)[1] <- c("values")
deplastimes_CARS <- na.omit(CARSData[,c("Time Deplas", "CPA")])
deplastimes_CARS$CPA <- paste0("b.deplas", deplastimes_CARS$CPA)
names(deplastimes_CARS)[1] <- c("values")

BPplastimes_CARS <- ggplot(plastimes_CARS, aes(y = values, x = CPA)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(color = "black", aes(fill = CPA)) + t1 +
  theme(legend.position="bottom")

BPdeplastimes_CARS <- ggplot(deplastimes_CARS, aes(y = values, x = CPA)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(color = "black", aes(fill = CPA)) + t1 +
  theme(legend.position="bottom")

timeBP_CARS <- list(plastimes_CARS, deplastimes_CARS)
timeBP_CARS <- dplyr::bind_rows(timeBP_CARS)
timeBP_CARS <- na.omit(reshape2::melt(timeBP_CARS))
timeBP_CARS <- timeBP_CARS[,c("CPA","value")]

KWtest_BF <- kruskal.test(value ~ CPA, timeBP_CARS) 
KWtest_BF # p<0.0001

dunn.all.CARS <- timeBP_CARS %>% dunn_test(value ~ CPA, p.adjust.method = "BY")

mypath <- file.path("D:", "Writing", "Cell response Paper", "Files for Sharing", "Bright Field Data", fsep="/",
                    paste("DunnAllCARS.csv"))
write.csv(dunn.all.CARS, mypath)

dunn.all.plas.CARS <- dunn.all.CARS[c(1,2,6),]
dunn.all.deplas.CARS <- dunn.all.CARS[c(13,14,15),]

PlasTimes_Stats <- boxplot(values ~ CPA, data = plastimes_CARS)
PlasTimes_BPSummary <- PlasTimes_Stats$stats
new_row<- t(as.data.frame(PlasTimes_Stats$n))
PlasTimes_BPSummary <- rbind(new_row, PlasTimes_BPSummary)
colnames(PlasTimes_BPSummary) <- PlasTimes_Stats$names
rownames(PlasTimes_BPSummary)<-c("n", "Min","First Quartile","Median","Third Quartile","Maximum")
mypath <- file.path("D:", "Writing", "Cell response Paper", "Files for Sharing", "Bright Field Data", fsep="/",
                    paste("PlasTimesChangePointBF.csv"))
write.csv(PlasTimes_BPSummary, mypath)

BPplastimes_CARS_stats <- BPplastimes_CARS +  t1 +
  stat_pvalue_manual(dunn.all.plas.CARS , y.position = 110, step.increase = 0.05, size = 5, tip.length = 0.01)+
  expand_limits(y = 150) +
  expand_limits(y = -30) +
  scale_y_continuous(breaks=seq(-30,2000,10)) +
  theme(legend.position="right") + 
  labs(y = "Time (s)", fill = "CPA") +
  scale_fill_manual(labels = c("DMSO", "Ethyelene Glycol", "Glycerol"),
                    values=c("#D0001C", "#00C2DC", "#D57124")) +
  theme(legend.title = element_text(size=15, 
                                    face="bold"))  

mypath <- file.path("D:", "Writing", "Cell response Paper", "Files for Sharing", "Bright Field Data", fsep="/",
                    paste("CARS_plastimes.png"))
png(file = mypath, units="in", width=6, height=6, res=300)
BPplastimes_CARS_stats
dev.off()

mypath <- file.path("D:", "Writing", "Cell response Paper", "Files for Sharing", "Bright Field Data", fsep="/",
                    paste("CARS_PlasTimesAndCP.png"))
png(file = mypath, units="in", width=6, height=6, res=300)
comboPlasCP <- ggarrange(BP_stats, BPplastimes_CARS_stats, 
                         labels = c("a.", "b."),
                         ncol = 2, nrow = 1)
dev.off()


BPdeplastimes_CARS_stats <- BPdeplastimes_CARS +  
  stat_pvalue_manual(dunn.all.deplas.CARS , y.position = 1250, step.increase = 0.05, size = 5, tip.length = 0.01) +
  expand_limits(y = 1400) +
  expand_limits(y = 0) +
  scale_y_continuous(breaks=seq(0,2000,100)) +
  theme(legend.text = element_text(size=12)) 

png(file = file.path("D:", "Writing", "Cell response Paper", "Files for Sharing", "Bright Field Data", fsep="/",
                     paste("CARS_ResponseType.png")), units="in", width=10, height=10, res=300)
BF_TimeBoxPlot_CARS <- ggarrange(BPplastimes_CARS_stats, BPdeplastimes_CARS_stats, 
                            labels = c("a.", "b."),
                            ncol = 2, nrow = 1) 
dev.off()
