# Normalize to Background using Change points not sigmoids as t=0
# load AllCPARawData.RData before running this

# Loading library conditions ####
#More info in NormalizingData.R
library(purrr)
library(tidyr)
library(dplyr)
library(readxl)
library(aod)

library(plotly)
library(RColorBrewer)
library(viridisLite)
library(viridis)
library(wesanderson)

library(gtsummary)
library(gridExtra)

library(reshape2)
library(data.table)
library(ggplot2)
library(ggpubr)

library(stats)
library(xlsx)

library(segmented)
library(reshape2)
# Functions ####
# Function to subtract the min and max of the background data from cell data
norm_to_BK <- function(x, cell_data, BK_min, BK_max) {
  (cell_data[[x]] - BK_min[[x]]) / (BK_max[[x]] - BK_min[[x]])
}

# Function to fit segmented linear pieces. Modified from ChangePoints.R
BPsegment <- function(DF, pred) {
  
  plot.new()
  sigfit <- list()
  fitList <- list()
  
  x <- DF$Time # The data frame has to contain Time as a column
  DF <- DF %>% dplyr::select(-contains("Time")) # remove the Time column from the data before looping
  par(mfrow=c(ncol(DF),2))
  for (i in 1:ncol(DF)) {
    y <- DF[i]
    y <- unlist(y)
    #fit simple linear regression model
    fit <-lm(y ~ x)
    #fit piece wise regression model to original model, estimating one or multiple breakpoints
    segmented.fit <- segmented(fit, seg.Z = ~x, psi = pred)
    fitList[[i]] <- segmented.fit
    summary_BP <- summary(segmented.fit)
    Estimates <- as.data.frame(summary_BP$psi)
    sigfit[[i]] <- (Estimates$Est.)
    
    plot(x, y, pch=16, col='steelblue')
    plot(segmented.fit, add=T)
    plot(resid(segmented.fit))
  }

  sigfit <- as.data.frame(do.call(rbind, sigfit))
  colnames(sigfit) <- NULL
  sigfit <- t(sigfit)
  colnames(sigfit) <- colnames(DF)
  sigfit <- as.data.frame((sigfit))

}

# This function works to bump up any data by an average - it adds the average to the data. 
# This assumes the data points are less than 0
bump <- function(x, data, average) {
  data[x] + abs(average[[x]])
}

# Make a list of all the minimums and maximums of the backgrounds ####
BK_Full_DMSO <- DMSOCARS %>% dplyr::select(-contains("Time"))
BK_Full_EG <- EGCARS %>% dplyr::select(-contains("Time"))
BK_Full_Glyc <- GlycCARS %>% dplyr::select(-contains("Time"))

BK_Full_DMSO <- DMSOCARS %>% dplyr::select(contains("BK"))
BK_Full_EG <- EGCARS %>% dplyr::select(contains("BK"))
BK_Full_Glyc <- GlycCARS %>% dplyr::select(contains("BK"))

BK_min_DMSO <- sapply(BK_Full_DMSO, min, na.rm=TRUE)
BK_min_DMSO <- as.data.frame(t(reshape2::melt(BK_min_DMSO, variable.name = "Background", value.name = "Min")))

BK_min_EG <- sapply(BK_Full_EG, min, na.rm=TRUE)
BK_min_EG <- as.data.frame(t(reshape2::melt(BK_min_EG, variable.name = "Background", value.name = "Min")))

BK_min_Glyc <- sapply(BK_Full_Glyc, min, na.rm=TRUE)
BK_min_Glyc <- as.data.frame(t(reshape2::melt(BK_min_Glyc, variable.name = "Background", value.name = "Min")))

BK_max_DMSO <- sapply(BK_Full_DMSO, max, na.rm=TRUE)
BK_max_DMSO <- as.data.frame(t(reshape2::melt(BK_max_DMSO, variable.name = "Background", value.name = "max")))

BK_max_EG <- sapply(BK_Full_EG, max, na.rm=TRUE)
BK_max_EG <- as.data.frame(t(reshape2::melt(BK_max_EG, variable.name = "Background", value.name = "max")))

BK_max_Glyc <- sapply(BK_Full_Glyc, max, na.rm=TRUE)
BK_max_Glyc <- as.data.frame(t(reshape2::melt(BK_max_Glyc, variable.name = "Background", value.name = "max")))


# Pull average min for each experiment ####
BK_min_Aug27E3 <- BK_min_Glyc %>% dplyr::select(contains('Aug27E3'))
BK_min_Aug27E3 <- rowMeans(BK_min_Aug27E3)

BK_min_Aug27E4 <- BK_min_Glyc %>% dplyr::select(contains('Aug27E4'))
BK_min_Aug27E4 <- rowMeans(BK_min_Aug27E4)

BK_min_Oct7E5 <- BK_min_Glyc %>% dplyr::select(contains('Oct7E5'))
BK_min_Oct7E5 <- rowMeans(BK_min_Oct7E5)

BK_min_Oct7E6 <- BK_min_Glyc %>% dplyr::select(contains('Oct7E6'))
BK_min_Oct7E6 <- rowMeans(BK_min_Oct7E6)

BK_min_Aug27E7 <- BK_min_EG %>% dplyr::select(contains('Aug27E7'))
BK_min_Aug27E7 <- rowMeans(BK_min_Aug27E7)

BK_min_Oct8E2 <- BK_min_EG %>% dplyr::select(contains('Oct8E2'))
BK_min_Oct8E2 <- rowMeans(BK_min_Oct8E2)

BK_min_Oct8E3 <- BK_min_EG %>% dplyr::select(contains('Oct8E3'))
BK_min_Oct8E3 <- rowMeans(BK_min_Oct8E3)

BK_min_Oct8E4 <- BK_min_EG %>% dplyr::select(contains('Oct8E4'))
BK_min_Oct8E4 <- rowMeans(BK_min_Oct8E4)

BK_min_July27E4 <- BK_min_DMSO %>% dplyr::select(contains('July27E4'))
BK_min_July27E4 <- rowMeans(BK_min_July27E4)

BK_min_July27E5 <- BK_min_DMSO %>% dplyr::select(contains('July27E5'))
BK_min_July27E5 <- rowMeans(BK_min_July27E5)

BK_min_July27E6 <- BK_min_DMSO %>% dplyr::select(contains('July27E6'))
BK_min_July27E6 <- rowMeans(BK_min_July27E6)

BK_min_Aug10E2 <- BK_min_DMSO %>% dplyr::select(contains('Aug10E2'))
BK_min_Aug10E2 <- rowMeans(BK_min_Aug10E2)

BK_min_Aug10E3 <- BK_min_DMSO %>% dplyr::select(contains('Aug10E3'))
BK_min_Aug10E3 <- rowMeans(BK_min_Aug10E3)

BK_min_Aug10E5 <- BK_min_DMSO %>% dplyr::select(contains('Aug10E5'))
BK_min_Aug10E5 <- rowMeans(BK_min_Aug10E5)

BK_min_Aug26E2 <- BK_min_DMSO %>% dplyr::select(contains('Aug26E2'))
BK_min_Aug26E2 <- rowMeans(BK_min_Aug26E2)

# Pull out average max of each experiment ####
BK_max_Aug27E3 <- BK_max_Glyc %>% dplyr::select(contains('Aug27E3'))
BK_max_Aug27E3 <- rowMeans(BK_max_Aug27E3)

BK_max_Aug27E4 <- BK_max_Glyc %>% dplyr::select(contains('Aug27E4'))
BK_max_Aug27E4 <- rowMeans(BK_max_Aug27E4)

BK_max_Oct7E5 <- BK_max_Glyc %>% dplyr::select(contains('Oct7E5'))
BK_max_Oct7E5 <- rowMeans(BK_max_Oct7E5)

BK_max_Oct7E6 <- BK_max_Glyc %>% dplyr::select(contains('Oct7E6'))
BK_max_Oct7E6 <- rowMeans(BK_max_Oct7E6)

BK_max_Aug27E7 <- BK_max_EG %>% dplyr::select(contains('Aug27E7'))
BK_max_Aug27E7 <- rowMeans(BK_max_Aug27E7)

BK_max_Oct8E2 <- BK_max_EG %>% dplyr::select(contains('Oct8E2'))
BK_max_Oct8E2 <- rowMeans(BK_max_Oct8E2)

BK_max_Oct8E3 <- BK_max_EG %>% dplyr::select(contains('Oct8E3'))
BK_max_Oct8E3 <- rowMeans(BK_max_Oct8E3)

BK_max_Oct8E4 <- BK_max_EG %>% dplyr::select(contains('Oct8E4'))
BK_max_Oct8E4 <- rowMeans(BK_max_Oct8E4)

BK_max_July27E4 <- BK_max_DMSO %>% dplyr::select(contains('July27E4'))
BK_max_July27E4 <- rowMeans(BK_max_July27E4)

BK_max_July27E5 <- BK_max_DMSO %>% dplyr::select(contains('July27E5'))
BK_max_July27E5 <- rowMeans(BK_max_July27E5)

BK_max_July27E6 <- BK_max_DMSO %>% dplyr::select(contains('July27E6'))
BK_max_July27E6 <- rowMeans(BK_max_July27E6)

BK_max_Aug10E2 <- BK_max_DMSO %>% dplyr::select(contains('Aug10E2'))
BK_max_Aug10E2 <- rowMeans(BK_max_Aug10E2)

BK_max_Aug10E3 <- BK_max_DMSO %>% dplyr::select(contains('Aug10E3'))
BK_max_Aug10E3 <- rowMeans(BK_max_Aug10E3)

BK_max_Aug10E5 <- BK_max_DMSO %>% dplyr::select(contains('Aug10E5'))
BK_max_Aug10E5 <- rowMeans(BK_max_Aug10E5)

BK_max_Aug26E2 <- BK_max_DMSO %>% dplyr::select(contains('Aug26E2'))
BK_max_Aug26E2 <- rowMeans(BK_max_Aug26E2)

# Min and max lists and full sheets ####
BK_min_list <- list(BK_min_Aug27E3, BK_min_Oct7E5, BK_min_Oct7E6, BK_min_Aug27E4, 
                    BK_min_Aug27E7, BK_min_Oct8E2, BK_min_Oct8E3, BK_min_Oct8E4,
                    BK_min_July27E4, BK_min_July27E5, BK_min_July27E6, BK_min_Aug10E2, BK_min_Aug10E3, 
                    BK_min_Aug10E5, BK_min_Aug26E2)
names(BK_min_list) <- c("Aug27E3", "Oct7E5", "Oct7E6", "Aug27E4", 
                        "Aug27E7", "Oct8E2", "Oct8E3", "Oct8E4",
                        "July27E4", "July27E5", "July27E6", "Aug10E2", "Aug10E3", 
                        "Aug10E5", "Aug26E2")
BK_min_list <- BK_min_list[order(names(BK_min_list))]

BK_max_list <- list(BK_max_Aug27E3, BK_max_Oct7E5, BK_max_Oct7E6, BK_max_Aug27E4, 
                    BK_max_Aug27E7, BK_max_Oct8E2, BK_max_Oct8E3, BK_max_Oct8E4,
                    BK_max_July27E4, BK_max_July27E5, BK_max_July27E6, BK_max_Aug10E2, BK_max_Aug10E3, 
                    BK_max_Aug10E5, BK_max_Aug26E2)
names(BK_max_list) <- c("Aug27E3", "Oct7E5", "Oct7E6", "Aug27E4", 
                        "Aug27E7", "Oct8E2", "Oct8E3", "Oct8E4",
                        "July27E4", "July27E5", "July27E6", "Aug10E2", "Aug10E3", 
                        "Aug10E5", "Aug26E2")
BK_max_list <- BK_max_list[order(names(BK_max_list))]

sheets_all <- c(sheetsDMSO, sheetsEG, sheetsglyc)
sheets_all <- sheets_all[order(names(sheets_all))]
sheets_all <- lapply(sheets_all, function(x) x[!names(x)=="Time"])
sheets_all <- lapply(sheets_all,function(x) replace(x,x==0,NA))


# Applying function to normalize data for background ####

norm_to_BK_cells <- lapply(names(sheets_all), norm_to_BK, cell_data=sheets_all, 
                           BK_min=BK_min_list, BK_max=BK_max_list)

names(norm_to_BK_cells) <- names(sheets_all)


# Add original Times back into norm_to_BK_cells ####
# Pull out times for each CPA experiment
TimeDMSO <- lapply(sheetsDMSO, function(x) x[names(x)=="Time"])
TimeEG <- lapply(sheetsEG, function(x) x[names(x)=="Time"])
Timeglyc <- lapply(sheetsglyc, function(x) x[names(x)=="Time"])

# Make a list of times for each experiment
sheetsTime <- c(TimeDMSO, TimeEG, Timeglyc)
sheetsTime <- sheetsTime[order(names(sheetsTime))]

# Combine Times list with normalized experiment list
normalized_to_BK <- mapply(c, norm_to_BK_cells, sheetsTime, SIMPLIFY = FALSE)

# Make a new DF from each experiment in list normalized_to_BK ####
list2env(normalized_to_BK,envir=.GlobalEnv)

# _int indicated that this is the first iteration of these dataframes, not the final version
Aug10E2_int <- dplyr::bind_rows(Aug10E2)
Aug10E3_int <- dplyr::bind_rows(Aug10E3)
Aug10E5_int <- dplyr::bind_rows(Aug10E5)
Aug26E2_int <- dplyr::bind_rows(Aug26E2)
Aug27E3_int <- dplyr::bind_rows(Aug27E3)
Aug27E4_int <- dplyr::bind_rows(Aug27E4)
Aug27E7_int <- dplyr::bind_rows(Aug27E7)
July27E4_int <- dplyr::bind_rows(July27E4)
July27E5_int <- dplyr::bind_rows(July27E5)
July27E6_int <- dplyr::bind_rows(July27E6)
Oct7E5_int <- dplyr::bind_rows(Oct7E5)
Oct7E6_int <- dplyr::bind_rows(Oct7E6)
Oct8E2_int <- dplyr::bind_rows(Oct8E2)
Oct8E3_int <- dplyr::bind_rows(Oct8E3)
Oct8E4_int <- dplyr::bind_rows(Oct8E4)

# Find Change Points of each background ####
# Pull out time and BK for each experiment
Aug10E2_BK <- Aug10E2_int %>% dplyr::select(contains("Time"), contains("BK"))
Aug10E3_BK <- Aug10E3_int %>% dplyr::select(contains("Time"), contains("BK"))
Aug10E5_BK <- Aug10E5_int %>% dplyr::select(contains("Time"), contains("BK"))
Aug26E2_BK <- Aug26E2_int %>% dplyr::select(contains("Time"), contains("BK"))
Aug27E3_BK <- Aug27E3_int %>% dplyr::select(contains("Time"), contains("BK"))
Aug27E4_BK <- Aug27E4_int %>% dplyr::select(contains("Time"), contains("BK"))
Aug27E7_BK <- Aug27E7_int %>% dplyr::select(contains("Time"), contains("BK"))
July27E4_BK <- July27E4_int %>% dplyr::select(contains("Time"), contains("BK"))
July27E5_BK <- July27E5_int %>% dplyr::select(contains("Time"), contains("BK"))
July27E6_BK <- July27E6_int %>% dplyr::select(contains("Time"), contains("BK"))
Oct7E5_BK <- Oct7E5_int %>% dplyr::select(contains("Time"), contains("BK"))
Oct7E6_BK <- Oct7E6_int %>% dplyr::select(contains("Time"), contains("BK"))
Oct8E2_BK <- Oct8E2_int %>% dplyr::select(contains("Time"), contains("BK"))
Oct8E3_BK <- Oct8E3_int %>% dplyr::select(contains("Time"), contains("BK"))
Oct8E4_BK <- Oct8E4_int %>% dplyr::select(contains("Time"), contains("BK"))

Aug10E2_BK <- subset(Aug10E2_BK, Time<400)
#I skrewed up the names of the background ROIs in the original excel file.
colnames(Aug10E3_BK)[2] <- "Aug10E3CARSBK1"
colnames(Aug10E3_BK)[3] <- "Aug10E3CARSBK2"
colnames(Aug10E3_BK)[4] <- "Aug10E3CARSBK3"
Aug10E3_BK <- subset(Aug10E3_BK, Time<300)
Aug10E3_BK <- subset(Aug10E3_BK, Time>75)
Aug10E5_BK <- subset(Aug10E5_BK, Time<400)
Aug26E2_BK <- subset(Aug26E2_BK, Time<300)
Aug27E3_BK <- subset(Aug27E3_BK, Time<600)
Aug27E4_BK <- subset(Aug27E4_BK, Time<600)
Aug27E4_BK <- subset(Aug27E4_BK, Time>200)
Aug27E7_BK <- subset(Aug27E7_BK, Time<600)
July27E4_BK <- subset(July27E4_BK, Time<400)
July27E4_BK <- subset(July27E4_BK, Time>50)
July27E5_BK <- subset(July27E5_BK, Time<500)
July27E5_BK <- subset(July27E5_BK, Time>50)
July27E6_BK <- subset(July27E6_BK, Time<600)
Oct7E5_BK <- subset(Oct7E5_BK, Time<800)
Oct7E5_BK <- subset(Oct7E5_BK, Time>100)
Oct7E6_BK <- subset(Oct7E6_BK, Time<500)
Oct7E6_BK <- subset(Oct7E6_BK, Time>200)
Oct8E2_BK <- subset(Oct8E2_BK, Time<600)
Oct8E2_BK <- subset(Oct8E2_BK, Time>100)
Oct8E3_BK <- subset(Oct8E3_BK, Time<550)
Oct8E3_BK <- subset(Oct8E3_BK, Time>150)
Oct8E4_BK <- subset(Oct8E4_BK, Time<450)
Oct8E4_BK <- subset(Oct8E4_BK, Time>200)

# Run change point function for every BK
Aug10E2_cpBK <- BPsegment(Aug10E2_BK, c(120, 150))
Aug10E3_cpBK <- BPsegment(Aug10E3_BK, c(175, 250)) 
Aug10E5_cpBK <- BPsegment(Aug10E5_BK, c(100, 175, 200)) 
Aug26E2_cpBK <- BPsegment(Aug26E2_BK, c(175, 200, 250)) 
Aug27E3_cpBK <- BPsegment(Aug27E3_BK, c(375, 450)) 
Aug27E4_cpBK <- BPsegment(Aug27E4_BK, c(350, 400, 450)) 
Aug27E7_cpBK <- BPsegment(Aug27E7_BK, c(400, 450, 500)) 
July27E4_cpBK <- BPsegment(July27E4_BK, c(175, 275, 300)) 
July27E5_cpBK <- BPsegment(July27E5_BK, c(225, 300, 400)) 
July27E6_cpBK <- BPsegment(July27E6_BK, c(250,350, 450))
Oct7E5_cpBK <- BPsegment(Oct7E5_BK, c(325, 400, 600)) 
Oct7E6_cpBK <- BPsegment(Oct7E6_BK, c(275, 350, 400) ) 
Oct8E2_cpBK <- BPsegment(Oct8E2_BK, c(350, 450, 500)) 
Oct8E3_cpBK <- BPsegment(Oct8E3_BK, c(350, 450, 500)) 
Oct8E4_cpBK <- BPsegment(Oct8E4_BK, c(325, 400))

# Find average of the first row (ie the first change point)
Aug10E2_cpBK_AV <- as.list(rowMeans(Aug10E2_cpBK[1,]))
Aug10E3_cpBK_AV <- as.list(rowMeans(Aug10E3_cpBK[1,])) 
Aug10E5_cpBK_AV <- as.list(rowMeans(Aug10E5_cpBK[1,])) 
Aug26E2_cpBK_AV <- as.list(rowMeans(Aug26E2_cpBK[1,])) 
Aug27E3_cpBK_AV <- as.list(rowMeans(Aug27E3_cpBK[1,])) 
Aug27E4_cpBK_AV <- as.list(rowMeans(Aug27E4_cpBK[1,])) 
Aug27E7_cpBK_AV <- as.list(rowMeans(Aug27E7_cpBK[1,])) 
July27E4_cpBK_AV <- as.list(rowMeans(July27E4_cpBK[1,])) 
July27E5_cpBK_AV <- as.list(rowMeans(July27E5_cpBK[1,])) 
July27E6_cpBK_AV <- as.list(rowMeans(July27E6_cpBK[1,]))
Oct7E5_cpBK_AV <- as.list(rowMeans(Oct7E5_cpBK[1,])) 
Oct7E6_cpBK_AV <- as.list(rowMeans(Oct7E6_cpBK[1,])) 
Oct8E2_cpBK_AV <- as.list(rowMeans(Oct8E2_cpBK[1,])) 
Oct8E3_cpBK_AV <- as.list(rowMeans(Oct8E3_cpBK[1,])) 
Oct8E4_cpBK_AV <- as.list(rowMeans(Oct8E4_cpBK[1,])) 

# Make average first cp into a list
fcplist_allexp <- list(Aug10E2_cpBK_AV, Aug10E3_cpBK_AV,
                        Aug10E5_cpBK_AV, Aug26E2_cpBK_AV,
                        Aug27E3_cpBK_AV, Aug27E4_cpBK_AV,
                        Aug27E7_cpBK_AV, July27E4_cpBK_AV,
                        July27E5_cpBK_AV, July27E6_cpBK_AV,
                        Oct7E5_cpBK_AV, Oct7E6_cpBK_AV,
                        Oct8E2_cpBK_AV, Oct8E3_cpBK_AV, 
                        Oct8E4_cpBK_AV)
names(fcplist_allexp) <- names(normalized_to_BK)
fcptable <- dplyr::bind_rows(fcplist_allexp)

# Move data based off of first change point. In this case xmid is the first change point ####
# f is loaded with AllCPARawData
newtimes <- lapply(names(sheetsTime), f, CPA = sheetsTime, xmid = fcplist_allexp)
names(newtimes) <- names(normalized_to_BK)

# add newtimes back into the normalized data and put into global environment ####
finalDF_norm_to_BK <- mapply(c, norm_to_BK_cells, newtimes, SIMPLIFY = FALSE)

list2env(finalDF_norm_to_BK,envir=.GlobalEnv)

# These replace the last set of global enviroment things from normalizing.
# _sft shows these are the shifted graphs
Aug10E2_sft <- dplyr::bind_rows(Aug10E2)
Aug10E3_sft <- dplyr::bind_rows(Aug10E3)
Aug10E5_sft <- dplyr::bind_rows(Aug10E5)
Aug26E2_sft <- dplyr::bind_rows(Aug26E2)
Aug27E3_sft <- dplyr::bind_rows(Aug27E3)
Aug27E4_sft <- dplyr::bind_rows(Aug27E4)
Aug27E7_sft <- dplyr::bind_rows(Aug27E7)
July27E4_sft <- dplyr::bind_rows(July27E4)
July27E5_sft <- dplyr::bind_rows(July27E5)
July27E6_sft <- dplyr::bind_rows(July27E6)
Oct7E5_sft <- dplyr::bind_rows(Oct7E5)
Oct7E6_sft <- dplyr::bind_rows(Oct7E6)
Oct8E2_sft <- dplyr::bind_rows(Oct8E2)
Oct8E3_sft <- dplyr::bind_rows(Oct8E3)
Oct8E4_sft <- dplyr::bind_rows(Oct8E4)

# Finding average of points before 0 and taking the average of them ####
# All the points less than -50, this should be far enough before t=0 that it won't include the increase
Aug10E2_sft_pts <- subset(Aug10E2_sft, Time<(-50)) %>% dplyr::select(-contains("Time"), -contains("BK"))
Aug10E3_sft_pts <- subset(Aug10E3_sft, Time<(-50)) %>% dplyr::select(-contains("Time"), -contains("BK"))
Aug10E5_sft_pts <- subset(Aug10E5_sft, Time<(-50)) %>% dplyr::select(-contains("Time"), -contains("BK"))
Aug26E2_sft_pts <- subset(Aug26E2_sft, Time<(-50)) %>% dplyr::select(-contains("Time"), -contains("BK"))
Aug27E3_sft_pts <- subset(Aug27E3_sft, Time<(-50)) %>% dplyr::select(-contains("Time"), -contains("BK"))
Aug27E4_sft_pts <- subset(Aug27E4_sft, Time<(-50)) %>% dplyr::select(-contains("Time"), -contains("BK"))
Aug27E7_sft_pts <- subset(Aug27E7_sft, Time<(-50)) %>% dplyr::select(-contains("Time"), -contains("BK"))
July27E4_sft_pts <- subset(July27E4_sft, Time<(-50)) %>% dplyr::select(-contains("Time"), -contains("BK"))
July27E5_sft_pts <- subset(July27E5_sft, Time<(-50)) %>% dplyr::select(-contains("Time"), -contains("BK"))
July27E6_sft_pts <- subset(July27E6_sft, Time<(-50)) %>% dplyr::select(-contains("Time"), -contains("BK"))
Oct7E5_sft_pts <- subset(Oct7E5_sft, Time<(-50)) %>% dplyr::select(-contains("Time"), -contains("BK"))
Oct7E6_sft_pts <- subset(Oct7E6_sft, Time<(-50)) %>% dplyr::select(-contains("Time"), -contains("BK"))
Oct8E2_sft_pts <- subset(Oct8E2_sft, Time<(-50)) %>% dplyr::select(-contains("Time"), -contains("BK"))
Oct8E3_sft_pts <- subset(Oct8E3_sft, Time<(-50)) %>% dplyr::select(-contains("Time"), -contains("BK"))
Oct8E4_sft_pts <- subset(Oct8E4_sft, Time<(-50)) %>% dplyr::select(-contains("Time"), -contains("BK"))

# Find averages of each column
Aug10E2_sft_avs <- sapply(Aug10E2_sft_pts, mean, na.rm =TRUE)
Aug10E3_sft_avs <- sapply(Aug10E3_sft_pts, mean, na.rm =TRUE)
Aug10E5_sft_avs <- sapply(Aug10E5_sft_pts, mean, na.rm =TRUE)
Aug26E2_sft_avs <- sapply(Aug26E2_sft_pts, mean, na.rm =TRUE)
Aug27E3_sft_avs <- sapply(Aug27E3_sft_pts, mean, na.rm =TRUE)
Aug27E4_sft_avs <- sapply(Aug27E4_sft_pts, mean, na.rm =TRUE)
Aug27E7_sft_avs <- sapply(Aug27E7_sft_pts, mean, na.rm =TRUE)
July27E4_sft_avs <- sapply(July27E4_sft_pts, mean, na.rm =TRUE)
July27E5_sft_avs <- sapply(July27E5_sft_pts, mean, na.rm =TRUE)
July27E6_sft_avs <- sapply(July27E6_sft_pts, mean, na.rm =TRUE)
Oct7E5_sft_avs <- sapply(Oct7E5_sft_pts, mean, na.rm =TRUE)
Oct7E6_sft_avs <- sapply(Oct7E6_sft_pts, mean, na.rm =TRUE)
Oct8E2_sft_avs <- sapply(Oct8E2_sft_pts, mean, na.rm =TRUE)
Oct8E3_sft_avs <- sapply(Oct8E3_sft_pts, mean, na.rm =TRUE)
Oct8E4_sft_avs <- sapply(Oct8E4_sft_pts, mean, na.rm =TRUE)

cell_sft_avs <- list(Aug10E2_sft_avs, Aug10E3_sft_avs, Aug10E5_sft_avs, 
                     Aug26E2_sft_avs, Aug27E3_sft_avs, Aug27E4_sft_avs, 
                     Aug27E7_sft_avs, July27E4_sft_avs, July27E5_sft_avs, 
                     July27E6_sft_avs, Oct7E5_sft_avs, Oct7E6_sft_avs, 
                     Oct8E2_sft_avs, Oct8E3_sft_avs, Oct8E4_sft_avs)

names(cell_sft_avs) <- names(newtimes)

# Only the cells that are normalized to background
Aug10E2_sft_cells <- Aug10E2_sft %>% dplyr::select(-contains("Time"), -contains("BK"))
Aug10E3_sft_cells <- Aug10E3_sft %>% dplyr::select(-contains("Time"), -contains("BK"))
Aug10E5_sft_cells <- Aug10E5_sft %>% dplyr::select(-contains("Time"), -contains("BK"))
Aug26E2_sft_cells <- Aug26E2_sft %>% dplyr::select(-contains("Time"), -contains("BK"))
Aug27E3_sft_cells <- Aug27E3_sft %>% dplyr::select(-contains("Time"), -contains("BK"))
Aug27E4_sft_cells <- Aug27E4_sft %>% dplyr::select(-contains("Time"), -contains("BK"))
Aug27E7_sft_cells <- Aug27E7_sft %>% dplyr::select(-contains("Time"), -contains("BK"))
July27E4_sft_cells <- July27E4_sft %>% dplyr::select(-contains("Time"), -contains("BK"))
July27E5_sft_cells <- July27E5_sft %>% dplyr::select(-contains("Time"), -contains("BK"))
July27E6_sft_cells <- July27E6_sft %>% dplyr::select(-contains("Time"), -contains("BK"))
Oct7E5_sft_cells <- Oct7E5_sft %>% dplyr::select(-contains("Time"), -contains("BK"))
Oct7E6_sft_cells <- Oct7E6_sft %>% dplyr::select(-contains("Time"), -contains("BK"))
Oct8E2_sft_cells <- Oct8E2_sft %>% dplyr::select(-contains("Time"), -contains("BK"))
Oct8E3_sft_cells <- Oct8E3_sft %>% dplyr::select(-contains("Time"), -contains("BK"))
Oct8E4_sft_cells <- Oct8E4_sft %>% dplyr::select(-contains("Time"), -contains("BK"))

# This list only contains the cell ROIs, not time or backgrounds
cell_sft_cells <- list(Aug10E2_sft_cells, Aug10E3_sft_cells, Aug10E5_sft_cells, 
                       Aug26E2_sft_cells, Aug27E3_sft_cells, Aug27E4_sft_cells, 
                       Aug27E7_sft_cells, July27E4_sft_cells, July27E5_sft_cells, 
                       July27E6_sft_cells, Oct7E5_sft_cells, Oct7E6_sft_cells, 
                       Oct8E2_sft_cells, Oct8E3_sft_cells, Oct8E4_sft_cells)

names(cell_sft_cells) <- names(newtimes)

# Bump up all the cell data using the bump function ####
# Bumped experiments only containing cells
Aug10E2_sft_bump <- as.data.frame(lapply(names(Aug10E2_sft_cells), bump, data = Aug10E2_sft_cells, average = Aug10E2_sft_avs))
Aug10E3_sft_bump <- as.data.frame(lapply(names(Aug10E3_sft_cells), bump, data = Aug10E3_sft_cells, average = Aug10E3_sft_avs))
Aug10E5_sft_bump <- as.data.frame(lapply(names(Aug10E5_sft_cells), bump, data = Aug10E5_sft_cells, average = Aug10E5_sft_avs))
Aug26E2_sft_bump <- as.data.frame(lapply(names(Aug26E2_sft_cells), bump, data = Aug26E2_sft_cells, average = Aug26E2_sft_avs))
Aug27E3_sft_bump <- as.data.frame(lapply(names(Aug27E3_sft_cells), bump, data = Aug27E3_sft_cells, average = Aug27E3_sft_avs))
Aug27E4_sft_bump <- as.data.frame(lapply(names(Aug27E4_sft_cells), bump, data = Aug27E4_sft_cells, average = Aug27E4_sft_avs))
Aug27E7_sft_bump <- as.data.frame(lapply(names(Aug27E7_sft_cells), bump, data = Aug27E7_sft_cells, average = Aug27E7_sft_avs))
July27E4_sft_bump <- as.data.frame(lapply(names(July27E4_sft_cells), bump, data = July27E4_sft_cells, average = July27E4_sft_avs))
July27E5_sft_bump <- as.data.frame(lapply(names(July27E5_sft_cells), bump, data = July27E5_sft_cells, average = July27E5_sft_avs))
July27E6_sft_bump <- as.data.frame(lapply(names(July27E6_sft_cells), bump, data = July27E6_sft_cells, average = July27E6_sft_avs))
Oct7E5_sft_bump <- as.data.frame(lapply(names(Oct7E5_sft_cells), bump, data = Oct7E5_sft_cells, average = Oct7E5_sft_avs))
Oct7E6_sft_bump <- as.data.frame(lapply(names(Oct7E6_sft_cells), bump, data = Oct7E6_sft_cells, average = Oct7E6_sft_avs))
Oct8E2_sft_bump <- as.data.frame(lapply(names(Oct8E2_sft_cells), bump, data = Oct8E2_sft_cells, average = Oct8E2_sft_avs))
Oct8E3_sft_bump <- as.data.frame(lapply(names(Oct8E3_sft_cells), bump, data = Oct8E3_sft_cells, average = Oct8E3_sft_avs))
Oct8E4_sft_bump <- as.data.frame(lapply(names(Oct8E4_sft_cells), bump, data = Oct8E4_sft_cells, average = Oct8E4_sft_avs))

# Bumped list containing ONLY cells
bumpedDF <- list(Aug10E2_sft_bump, Aug10E3_sft_bump, Aug10E5_sft_bump,
                 Aug26E2_sft_bump, Aug27E3_sft_bump, Aug27E4_sft_bump,
                 Aug27E7_sft_bump, July27E4_sft_bump, July27E5_sft_bump,
                 July27E6_sft_bump, Oct7E5_sft_bump, Oct7E6_sft_bump, 
                 Oct8E2_sft_bump, Oct8E3_sft_bump, Oct8E4_sft_bump)
names(bumpedDF) <- names(newtimes)

# Only the background and time
Aug10E2_sft_BK_T <- Aug10E2_sft %>% dplyr::select(contains("Time"), contains("BK"))
Aug10E3_sft_BK_T <- Aug10E3_sft %>% dplyr::select(contains("Time"), contains("BK"))
Aug10E5_sft_BK_T <- Aug10E5_sft %>% dplyr::select(contains("Time"), contains("BK"))
Aug26E2_sft_BK_T <- Aug26E2_sft %>% dplyr::select(contains("Time"), contains("BK"))
Aug27E3_sft_BK_T <- Aug27E3_sft %>% dplyr::select(contains("Time"), contains("BK"))
Aug27E4_sft_BK_T <- Aug27E4_sft %>% dplyr::select(contains("Time"), contains("BK"))
Aug27E7_sft_BK_T <- Aug27E7_sft %>% dplyr::select(contains("Time"), contains("BK"))
July27E4_sft_BK_T <- July27E4_sft %>% dplyr::select(contains("Time"), contains("BK"))
July27E5_sft_BK_T <- July27E5_sft %>% dplyr::select(contains("Time"), contains("BK"))
July27E6_sft_BK_T <- July27E6_sft %>% dplyr::select(contains("Time"), contains("BK"))
Oct7E5_sft_BK_T <- Oct7E5_sft %>% dplyr::select(contains("Time"), contains("BK"))
Oct7E6_sft_BK_T <- Oct7E6_sft %>% dplyr::select(contains("Time"), contains("BK"))
Oct8E2_sft_BK_T <- Oct8E2_sft %>% dplyr::select(contains("Time"), contains("BK"))
Oct8E3_sft_BK_T <- Oct8E3_sft %>% dplyr::select(contains("Time"), contains("BK"))
Oct8E4_sft_BK_T <- Oct8E4_sft %>% dplyr::select(contains("Time"), contains("BK"))

# Bumped list containing ONLY background and time
BK_T_DF <- list(Aug10E2_sft_BK_T, Aug10E3_sft_BK_T, Aug10E5_sft_BK_T,
                Aug26E2_sft_BK_T, Aug27E3_sft_BK_T, Aug27E4_sft_BK_T,
                Aug27E7_sft_BK_T, July27E4_sft_BK_T, July27E5_sft_BK_T,
                July27E6_sft_BK_T, Oct7E5_sft_BK_T, Oct7E6_sft_BK_T, 
                Oct8E2_sft_BK_T, Oct8E3_sft_BK_T, Oct8E4_sft_BK_T)
names(BK_T_DF) <- names(newtimes)

# Put the background and times back into the bumped data set and add to global envir ####
FINALDF_norm_bump <- mapply(c, bumpedDF, BK_T_DF, SIMPLIFY = FALSE)
list2env(FINALDF_norm_bump,envir=.GlobalEnv)

# These should be the final products!!!!!!!!!!
Aug10E2 <- dplyr::bind_rows(Aug10E2)
Aug10E3 <- dplyr::bind_rows(Aug10E3)
Aug10E5 <- dplyr::bind_rows(Aug10E5)
Aug26E2 <- dplyr::bind_rows(Aug26E2)
Aug27E3 <- dplyr::bind_rows(Aug27E3)
Aug27E4 <- dplyr::bind_rows(Aug27E4)
Aug27E7 <- dplyr::bind_rows(Aug27E7)
July27E4 <- dplyr::bind_rows(July27E4)
July27E5 <- dplyr::bind_rows(July27E5)
July27E6 <- dplyr::bind_rows(July27E6)
Oct7E5 <- dplyr::bind_rows(Oct7E5)
Oct7E6 <- dplyr::bind_rows(Oct7E6)
Oct8E2 <- dplyr::bind_rows(Oct8E2)
Oct8E3 <- dplyr::bind_rows(Oct8E3)
Oct8E4 <- dplyr::bind_rows(Oct8E4)

# Melt the data to plot the experiments ####

meltAug10E2 <- reshape2::melt(Aug10E2, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
meltAug10E3 <- reshape2::melt(Aug10E3, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
meltAug10E5 <- reshape2::melt(Aug10E5, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
meltAug26E2 <- reshape2::melt(Aug26E2, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
meltAug27E3 <- reshape2::melt(Aug27E3, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
meltAug27E4 <- reshape2::melt(Aug27E4, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
meltAug27E7 <- reshape2::melt(Aug27E7, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
meltJuly27E4 <- reshape2::melt(July27E4, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
meltJuly27E5 <- reshape2::melt(July27E5, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
meltJuly27E6 <- reshape2::melt(July27E6, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
meltOct7E5 <- reshape2::melt(Oct7E5, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
meltOct7E6 <- reshape2::melt(Oct7E6, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
meltOct8E2 <- reshape2::melt(Oct8E2, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
meltOct8E3 <- reshape2::melt(Oct8E3, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
meltOct8E4 <- reshape2::melt(Oct8E4, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")

# Plot each experiment ####
colgry <- c("#000000", "#969696")
mypal <- colorRampPalette(rev(colgry))
colgb <- c("#2b8cbe", "#4eb3d3", "#238b45", "#006d2c")
mypal2 <- colorRampPalette(rev(colgb))
brewpal <- colorRampPalette(brewer.pal(12, "Paired"))
rainpal <- colorRampPalette(hsv(seq(0,1 - 1/5,length.out = 100), 0.60 , 0.85))
rain2pal <- c("#001219", "#005f73", "#0a9396", 
              "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702", "#bb3e03", "#ae2012", "#9b2226","#f0795c")
rain3pal <- c("#591f0d", "#8d3316", "#A8570F", "#c57d07", "#f8ae35", "#25562e", "#378145", "#0a4d65", "#4ca893", 
              "#9b5bb0", "#722c82")
rain4pal <- c("#96331c", "#e4500c", "#e9730c", "#e6a00a", "#dcc40f", "#7fd05d", 
              "#20d4a7", "#68aae3", "#254AB5", "#3B248E", "#A36FB5")

plotAug10E2 <- ggplot(data = meltAug10E2, aes(Time, Intensity, col = ROI)) + geom_point() + 
  geom_line(aes(Time, Intensity, col = ROI), meltAug10E2)  +
  scale_colour_manual(values = c(rain4pal[1:4], mypal(2))) +
  theme_classic() + ggtitle("Aug10 E2") + xlab("Time (s)") + ylab("Intensity (a.u.)") +
  theme(plot.title = element_text(size=10, face = "bold"), 
        panel.grid.major.x=element_line(), legend.position =  "right") +
  xlim(-200,600) + 
  guides(color = guide_legend(ncol = 2, bycol = TRUE)) +
  #geom_line(aes(faketime, Aug10E2_sfBK_sigmoid), Aug10E2_sfBK_sigmoid, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -25.9, linetype="dashed", color = "red") +
  geom_vline(xintercept = -4.1, linetype="dashed", color = "red") +
  geom_vline(xintercept = 34.1, linetype="dashed", color = "red") +
  geom_vline(xintercept = 94.1, linetype="dashed", color = "red") +
  geom_vline(xintercept = 244.1, linetype="dashed", color = "red") +
  scale_x_continuous(breaks = seq(-200, 600, 20), labels = label_at(200), limits = c(-200, 600))

plotAug10E3 <- ggplot(data = meltAug10E3, aes(Time, Intensity, col = ROI)) + geom_point() + 
  geom_line(aes(Time, Intensity, col = ROI), meltAug10E3) +
  scale_colour_manual(values = c(rain4pal[1:11], mypal(3))) +
  theme_classic() + ggtitle("Aug10 E3") + xlab("Time (s)") + ylab("Intensity (a.u.)") +
  theme(plot.title = element_text(size=10, face = "bold"), 
        panel.grid.major.x=element_line(), legend.position =  "right") +
  xlim(-200,600) + 
  guides(color = guide_legend(ncol = 2, bycol = TRUE)) +
  #geom_line(aes(faketime, Aug10E2_sfBK_sigmoid), Aug10E2_sfBK_sigmoid, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -19.8, linetype="dashed", color = "red") +
  geom_vline(xintercept = 0.2, linetype="dashed", color = "red") +
  geom_vline(xintercept = 20.2, linetype="dashed", color = "red") +
  geom_vline(xintercept = 60.2, linetype="dashed", color = "red") +
  geom_vline(xintercept = 160.2, linetype="dashed", color = "red") +
  scale_x_continuous(breaks = seq(-200, 600, 20), labels = label_at(200), limits = c(-200, 600))

plotAug10E5 <- ggplot(data = meltAug10E5, aes(Time, Intensity, col = ROI)) + geom_point() + 
  geom_line(aes(Time, Intensity, col = ROI), meltAug10E5) +
  scale_colour_manual(values = c(rain4pal[1:2], mypal(3))) +
  theme_classic() + ggtitle("Aug10 E5") + xlab("Time (s)") + ylab("Intensity (a.u.)") +
  theme(plot.title = element_text(size=10, face = "bold"), 
        panel.grid.major.x=element_line(), legend.position =  "right") +
  xlim(-200,600) + 
  guides(color = guide_legend(ncol = 2, bycol = TRUE)) +
  #geom_line(aes(faketime, Aug10E5_sfBK_sigmoid), Aug10E5_sfBK_sigmoid, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -19.7, linetype="dashed", color = "red") +
  geom_vline(xintercept = 0.3, linetype="dashed", color = "red") +
  geom_vline(xintercept = 20.3, linetype="dashed", color = "red") +
  geom_vline(xintercept = 60.3, linetype="dashed", color = "red") +
  geom_vline(xintercept = 160.3, linetype="dashed", color = "red") +
  scale_x_continuous(breaks = seq(-200, 600, 20), labels = label_at(200), limits = c(-200, 600))

plotAug26E2 <- ggplot(data = meltAug26E2, aes(Time, Intensity, col = ROI)) + geom_point() + 
  geom_line(aes(Time, Intensity, col = ROI), meltAug26E2) +
  scale_colour_manual(values = c(rain4pal[1:7], mypal(3))) +
  theme_classic() + ggtitle("Aug26 E2") + xlab("Time (s)") + ylab("Intensity (a.u.)") +
  theme(plot.title = element_text(size=10, face = "bold"), 
        panel.grid.major.x=element_line(), legend.position =  "right") +
  xlim(-200,600) + 
  guides(color = guide_legend(ncol = 2, bycol = TRUE)) +
  #geom_line(aes(faketime, Aug26E2_sfBK_sigmoid), Aug26E2_sfBK_sigmoid, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -18.6, linetype="dashed", color = "red") +
  geom_vline(xintercept = 1.4, linetype="dashed", color = "red") +
  geom_vline(xintercept = 21.4, linetype="dashed", color = "red") +
  geom_vline(xintercept = 61.4, linetype="dashed", color = "red") +
  geom_vline(xintercept = 161.4, linetype="dashed", color = "red") +
  scale_x_continuous(breaks = seq(-200, 600, 20), labels = label_at(200), limits = c(-200, 600))

plotAug27E3 <- ggplot(data = meltAug27E3, aes(Time, Intensity, col = ROI)) + geom_point() + 
  geom_line(aes(Time, Intensity, col = ROI), meltAug27E3) +
  scale_colour_manual(values = c(rain4pal[1:8], mypal(3))) +
  theme_classic() + ggtitle("Aug27 E3") + xlab("Time (s)") + ylab("Intensity (a.u.)") +
  theme(plot.title = element_text(size=10, face = "bold"), 
        panel.grid.major.x=element_line(), legend.position =  "right") +
  xlim(-200,600) + 
  guides(color = guide_legend(ncol = 2, bycol = TRUE)) +
  #geom_line(aes(faketime, Aug27E3_sfBK_sigmoid), Aug27E3_sfBK_sigmoid, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -14.5, linetype="dashed", color = "red") +
  geom_vline(xintercept = 5.5, linetype="dashed", color = "red") +
  geom_vline(xintercept = 25.5, linetype="dashed", color = "red") +
  geom_vline(xintercept = 65.5, linetype="dashed", color = "red") +
  geom_vline(xintercept = 165.5, linetype="dashed", color = "red") +
  scale_x_continuous(breaks = seq(-200, 600, 20), labels = label_at(200), limits = c(-200, 600))

plotAug27E4 <- ggplot(data = meltAug27E4, aes(Time, Intensity, col = ROI)) + geom_point() + 
  geom_line(aes(Time, Intensity, col = ROI), meltAug27E4) +
  scale_colour_manual(values = c(rain4pal[1:7], mypal(3))) +
  theme_classic() + ggtitle("Aug27 E4") + xlab("Time (s)") + ylab("Intensity (a.u.)") +
  theme(plot.title = element_text(size=10, face = "bold"), 
        panel.grid.major.x=element_line(), legend.position =  "right") +
  xlim(-200,600) + 
  guides(color = guide_legend(ncol = 2, bycol = TRUE)) +
  #geom_line(aes(faketime, Aug27E4_sfBK_sigmoid), Aug27E4_sfBK_sigmoid, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -20.2, linetype="dashed", color = "red") +
  geom_vline(xintercept = -0.2, linetype="dashed", color = "red") +
  geom_vline(xintercept = 19.8, linetype="dashed", color = "red") +
  geom_vline(xintercept = 59.8, linetype="dashed", color = "red") +
  geom_vline(xintercept = 159.8, linetype="dashed", color = "red") +
  scale_x_continuous(breaks = seq(-200, 600, 20), labels = label_at(200), limits = c(-200, 600))

plotAug27E7 <- ggplot(data = meltAug27E7, aes(Time, Intensity, col = ROI)) + geom_point() + 
  geom_line(aes(Time, Intensity, col = ROI), meltAug27E7) +
  scale_colour_manual(values = c(rain4pal[1:10], mypal(3))) +
  theme_classic() + ggtitle("Aug27 E7") + xlab("Time (s)") + ylab("Intensity (a.u.)") +
  theme(plot.title = element_text(size=10, face = "bold"), 
        panel.grid.major.x=element_line(), legend.position =  "right") +
  xlim(-200,600) + 
  guides(color = guide_legend(ncol = 2, bycol = TRUE)) +
  #geom_line(aes(faketime, Aug27E7_sfBK_sigmoid), Aug27E7_sfBK_sigmoid, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -20.4, linetype="dashed", color = "red") +
  geom_vline(xintercept = -0.4, linetype="dashed", color = "red") +
  geom_vline(xintercept = 19.6, linetype="dashed", color = "red") +
  geom_vline(xintercept = 59.6, linetype="dashed", color = "red") +
  geom_vline(xintercept = 159.6, linetype="dashed", color = "red") +
  scale_x_continuous(breaks = seq(-200, 600, 20), labels = label_at(200), limits = c(-200, 600))

plotJuly27E4 <- ggplot(data = meltJuly27E4, aes(Time, Intensity, col = ROI)) + geom_point() + 
  geom_line(aes(Time, Intensity, col = ROI), meltJuly27E4) +
  scale_colour_manual(values = c(rain4pal[1:5], mypal(3))) +
  theme_classic() + ggtitle("July27 E4") + xlab("Time (s)") + ylab("Intensity (a.u.)") +
  theme(plot.title = element_text(size=10, face = "bold"), 
        panel.grid.major.x=element_line(), legend.position =  "right") +
  xlim(-200,600) + 
  guides(color = guide_legend(ncol = 2, bycol = TRUE)) +
  #geom_line(aes(faketime, July27E4_sfBK_sigmoid), July27E4_sfBK_sigmoid, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -18.5, linetype="dashed", color = "red") +
  geom_vline(xintercept = 1.5, linetype="dashed", color = "red") +
  geom_vline(xintercept = 21.5, linetype="dashed", color = "red") +
  geom_vline(xintercept = 61.5, linetype="dashed", color = "red") +
  geom_vline(xintercept = 161.5, linetype="dashed", color = "red") +
  scale_x_continuous(breaks = seq(-200, 600, 20), labels = label_at(200), limits = c(-200, 600))

plotJuly27E5 <- ggplot(data = meltJuly27E5, aes(Time, Intensity, col = ROI)) + geom_point() + 
  geom_line(aes(Time, Intensity, col = ROI), meltJuly27E5) +
  scale_colour_manual(values = c(rain4pal[1:8], mypal(3))) +
  theme_classic() + ggtitle("July27 E5") + xlab("Time (s)") + ylab("Intensity (a.u.)") +
  theme(plot.title = element_text(size=10, face = "bold"), 
        panel.grid.major.x=element_line(), legend.position =  "right") +
  xlim(-200,600) + 
  guides(color = guide_legend(ncol = 2, bycol = TRUE)) +
  #geom_line(aes(faketime, July27E5_sfBK_sigmoid), July27E5_sfBK_sigmoid, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -20.4, linetype="dashed", color = "red") +
  geom_vline(xintercept = -0.4, linetype="dashed", color = "red") +
  geom_vline(xintercept = 19.6, linetype="dashed", color = "red") +
  geom_vline(xintercept = 59.6, linetype="dashed", color = "red") +
  geom_vline(xintercept = 159.6, linetype="dashed", color = "red") +
  scale_x_continuous(breaks = seq(-200, 600, 20), labels = label_at(200), limits = c(-200, 600))

plotJuly27E6 <- ggplot(data = meltJuly27E6, aes(Time, Intensity, col = ROI)) + geom_point() + 
  geom_line(aes(Time, Intensity, col = ROI), meltJuly27E6) +
  scale_colour_manual(values = c(rain4pal[1:10], mypal(3))) +
  theme_classic() + ggtitle("July27 E6") + xlab("Time (s)") + ylab("Intensity (a.u.)") +
  theme(plot.title = element_text(size=10, face = "bold"), 
        panel.grid.major.x=element_line(), legend.position =  "right") +
  xlim(-200,600) + 
  guides(color = guide_legend(ncol = 2, bycol = TRUE)) +
  #geom_line(aes(faketime, July27E6_sfBK_sigmoid), July27E6_sfBK_sigmoid, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -31.2, linetype="dashed", color = "red") +
  geom_vline(xintercept = -1.2, linetype="dashed", color = "red") +
  geom_vline(xintercept = 28.8, linetype="dashed", color = "red") +
  geom_vline(xintercept = 88.8, linetype="dashed", color = "red") +
  geom_vline(xintercept = 238.8, linetype="dashed", color = "red") +
  scale_x_continuous(breaks = seq(-200, 600, 20), labels = label_at(200), limits = c(-200, 600))

plotOct7E5 <- ggplot(data = meltOct7E5, aes(Time, Intensity, col = ROI)) + geom_point() + 
  geom_line(aes(Time, Intensity, col = ROI), meltOct7E5) +
  scale_colour_manual(values = c(rain4pal[1:5], mypal(3))) +
  theme_classic() + ggtitle("Oct7 E5") + xlab("Time (s)") + ylab("Intensity (a.u.)") +
  theme(plot.title = element_text(size=10, face = "bold"), 
        panel.grid.major.x=element_line(), legend.position =  "right") +
  xlim(-200,600) + 
  guides(color = guide_legend(ncol = 2, bycol = TRUE)) +
  #geom_line(aes(faketime, Oct7E5_sfBK_sigmoid), Oct7E5_sfBK_sigmoid, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -22.7, linetype="dashed", color = "red") +
  geom_vline(xintercept = -3.1, linetype="dashed", color = "red") +
  geom_vline(xintercept = 16.5, linetype="dashed", color = "red") +
  geom_vline(xintercept = 55.6, linetype="dashed", color = "red") +
  geom_vline(xintercept = 153.5, linetype="dashed", color = "red") +
  scale_x_continuous(breaks = seq(-200, 600, 20), labels = label_at(200), limits = c(-200, 600))

plotOct7E6 <- ggplot(data = meltOct7E6, aes(Time, Intensity, col = ROI)) + geom_point() + 
  geom_line(aes(Time, Intensity, col = ROI), meltOct7E6) +
  scale_colour_manual(values = c(rain4pal[1:7], mypal(3))) +
  theme_classic() + ggtitle("Oct7 E6") + xlab("Time (s)") + ylab("Intensity (a.u.)") +
  theme(plot.title = element_text(size=10, face = "bold"), 
        panel.grid.major.x=element_line(), legend.position =  "right") +
  xlim(-200,600) + 
  guides(color = guide_legend(ncol = 2, bycol = TRUE)) +
  #geom_line(aes(faketime, Oct7E6_sfBK_sigmoid), Oct7E6_sfBK_sigmoid, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -18.8, linetype="dashed", color = "red") +
  geom_vline(xintercept = 0.8, linetype="dashed", color = "red") +
  geom_vline(xintercept = 20.4, linetype="dashed", color = "red") +
  geom_vline(xintercept = 59.5, linetype="dashed", color = "red") +
  geom_vline(xintercept = 157.4, linetype="dashed", color = "red") +
  scale_x_continuous(breaks = seq(-200, 600, 20), labels = label_at(200), limits = c(-200, 600))

plotOct8E2 <- ggplot(data = meltOct8E2, aes(Time, Intensity, col = ROI)) + geom_point() + 
  geom_line(aes(Time, Intensity, col = ROI), meltOct8E2) +
  scale_colour_manual(values = c(rain4pal[1:10], mypal(3))) +
  theme_classic() + ggtitle("Oct8 E2") + xlab("Time (s)") + ylab("Intensity (a.u.)") +
  theme(plot.title = element_text(size=10, face = "bold"), 
        panel.grid.major.x=element_line(), legend.position =  "right") +
  xlim(-200,600) + 
  guides(color = guide_legend(ncol = 2, bycol = TRUE)) +
  #geom_line(aes(faketime, Oct8E2_sfBK_sigmoid), Oct8E2_sfBK_sigmoid, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -12.9, linetype="dashed", color = "red") +
  geom_vline(xintercept = -3.3, linetype="dashed", color = "red") +
  geom_vline(xintercept = 22.0, linetype="dashed", color = "red") +
  geom_vline(xintercept = 61.2, linetype="dashed", color = "red") +
  geom_vline(xintercept = 159.1, linetype="dashed", color = "red") +
  scale_x_continuous(breaks = seq(-200, 600, 20), labels = label_at(200), limits = c(-200, 600))

plotOct8E3 <- ggplot(data = meltOct8E3, aes(Time, Intensity, col = ROI)) + geom_point() + 
  geom_line(aes(Time, Intensity, col = ROI), meltOct8E3) +
  scale_colour_manual(values = c(rain4pal[1:8], mypal(3))) +
  theme_classic() + ggtitle("Oct8 E3") + xlab("Time (s)") + ylab("Intensity (a.u.)") +
  theme(plot.title = element_text(size=10, face = "bold"), 
        panel.grid.major.x=element_line(), legend.position =  "right") +
  xlim(-200,600) + 
  guides(color = guide_legend(ncol = 2, bycol = TRUE)) +
  #geom_line(aes(faketime, Oct8E3_sfBK_sigmoid), Oct8E3_sfBK_sigmoid, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -17.1, linetype="dashed", color = "red") +
  geom_vline(xintercept = 2.4, linetype="dashed", color = "red") +
  geom_vline(xintercept = 22.0, linetype="dashed", color = "red") +
  geom_vline(xintercept = 61.2, linetype="dashed", color = "red") +
  geom_vline(xintercept = 159.1, linetype="dashed", color = "red") +
  scale_x_continuous(breaks = seq(-200, 600, 20), labels = label_at(200), limits = c(-200, 600))

plotOct8E4 <- ggplot(data = meltOct8E4, aes(Time, Intensity, col = ROI)) + geom_point() + 
  geom_line(aes(Time, Intensity, col = ROI), meltOct8E4) +
  scale_colour_manual(values = c(rain4pal[1:9], mypal(3))) +
  theme_classic() + ggtitle("Oct8 E4") + xlab("Time (s)") + ylab("Intensity (a.u.)") +
  theme(plot.title = element_text(size=10, face = "bold"), 
        panel.grid.major.x=element_line(), legend.position =  "right") +
  xlim(-200,600) + 
  guides(color = guide_legend(ncol = 2, bycol = TRUE)) +
  #geom_line(aes(faketime, Oct8E4_sfBK_sigmoid), Oct8E4_sfBK_sigmoid, color = "black", linetype = "dashed") +
  geom_vline(xintercept = -15.1, linetype="dashed", color = "red") +
  geom_vline(xintercept = 4.4, linetype="dashed", color = "red") +
  geom_vline(xintercept = 24.0, linetype="dashed", color = "red") +
  geom_vline(xintercept = 63.2, linetype="dashed", color = "red") +
  geom_vline(xintercept = 161.1, linetype="dashed", color = "red") +
  scale_x_continuous(breaks = seq(-200, 600, 20), labels = label_at(200), limits = c(-200, 600))

# Saving every experiment ####
library(egg) #this one does not have the same ggarrange function as ggpubr

 

# Putting plots on the same page ####
# I changed these every time I changed the plots and saved them as new files.
# DMSO experiments are: "July27E4" "July27E5" "July27E6" "Aug10E2"  "Aug10E3"  "Aug10E5"  "Aug26E2" 
# EG experiments are: "Oct8E2"  "Oct8E3"  "Oct8E4"  "Aug27E7"
# Glycerol experiments are: "Oct7E5"  "Oct7E6"  "Aug27E3" "Aug27E4"
library(egg) #this one does not have the same ggarrange function as ggpubr

png(file = file.path("C:", "Users", "fionn", "Dropbox", "Levinger Group", "Levinger Group members", "Fionna", 
                     "Research", "Writing", "Cell response Paper", "Figures", fsep="/",
                     paste("CPDMSO.png")), units="in", width=10, height=17.5, res=300)
allBKnormDMSOexp <- ggarrange(plotJuly27E4, plotJuly27E5, plotJuly27E6, plotAug10E2, 
                              plotAug10E3, plotAug10E5, plotAug26E2, 
                              labels = c("a.", "b.", "c.", "d.", "e.", "f.", "g."),
                              ncol = 1, nrow = 7) 

annotate_figure(allBKnormDMSOexp, top = text_grob("Dimethyl Sulfoxide Cell Intensity (Background Normalized)",
                                                  color = "black", face = "bold", size = 14))
dev.off()

png(file = file.path("C:", "Users", "fionn", "Dropbox", "Levinger Group", "Levinger Group members", "Fionna", 
                     "Research", "Writing", "Cell response Paper", "Figures", fsep="/",
                     paste("CPEG.png")), units="in", width=10, height=10, res=300)
allBKnormEGexp <- ggarrange(plotAug27E7, plotOct8E2, plotOct8E3, plotOct8E4, 
                            labels = c("a.", "b.", "c.", "d."),
                            ncol = 1, nrow = 4) 

annotate_figure(allBKnormEGexp, top = text_grob("Ethylene Glycol Cell Intensity (Background Normalized)",
                                                color = "black", face = "bold", size = 14))
dev.off()


png(file = file.path("C:", "Users", "fionn", "Dropbox", "Levinger Group", "Levinger Group members", "Fionna", 
                     "Research", "Writing", "Cell response Paper", "Figures", fsep="/",
                     paste("CPGlyc.png")), units="in", width=10, height=10, res=300)
allBKnormGLYCexp <- ggarrange(plotAug27E3, plotAug27E4, plotOct7E5, plotOct7E6, 
                              labels = c("a.", "b.", "c.", "d."),
                              ncol = 1, nrow = 4) 

annotate_figure(allBKnormGLYCexp, top = text_grob("Glycerol Cell Intensity (Background Normalized)",
                                                  color = "black", face = "bold", size = 14))
dev.off()

