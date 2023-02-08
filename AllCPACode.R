#For this code to work, load AllCPARawData.RData

#For each CPA, this code generates the background graphs then the Cell graphs

### Functions and loading conditions ####
#More info in NormalizingData.R
library(purrr)
library(tidyr)
library(dplyr)
library(readxl)
library(aod)

library(plotly)
library(RColorBrewer)

library(gtsummary)
library(gridExtra)

library(reshape2)
library(data.table)
library(ggplot2)
library(ggpubr)

library(stats)
library(xlsx)

#importing data - possible to make this more user-friendly using "readline"
read_allsheets_excel <- function(filename, tibble = TRUE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#choose your data
importdata <- function (AF, CARS) { 
  df <- read_allsheets_excel(file.choose())
  null_AF_df <- df$Autofluorescence
  null_AF_df[, 2:ncol(null_AF_df)][null_AF_df[, 2:ncol(null_AF_df)] == 0] <- NA
  AF <<- null_AF_df
  
  null_CARS_df <- df$CARS
  null_CARS_df[, 2:ncol(null_CARS_df)][null_CARS_df[, 2:ncol(null_CARS_df)] == 0] <- NA
  CARS <<- null_CARS_df
}

#normalize data - define Min-Max normalization function: This is the function I usually use in Excel
#this definition takes into account the NA present in the data set.
normalize <- function(x) {
  (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))
}

#Function that subracts xmid from the Time column of whatever CPA data you want to work with.
f <- function(x, CPA, xmid) {
  (CPA[[x]]["Time"] - xmid[[x]])
}


### DMSO ####

# Subset All Data 
{#
  new_DMSOAF <- subset(DMSOAF, Time<500)
  new_DMSOCARS <- subset(DMSOCARS, Time<500)
  
  ggplot(data = new_DMSOCARS, aes(Time, July27E4CARSC1)) + geom_point()
  
  #use lapply to apply it to the columns in the data set - 
  #you don't want to include time column, we will add that back in later
  
  norm_DMSOAF <- as.data.frame(lapply(new_DMSOAF[2:ncol(DMSOAF)], normalize))
  norm_DMSOCARS <- as.data.frame(lapply(new_DMSOCARS[2:ncol(DMSOCARS)], normalize))
  
  #add the time columns back into the data, and reorder them so they are first
  
  norm_DMSOAF$Time <- new_DMSOAF$Time
  norm_DMSOAF <- norm_DMSOAF %>% select(Time, everything())
  
  norm_DMSOCARS$Time <- new_DMSOCARS$Time
  norm_DMSOCARS <- norm_DMSOCARS %>% select(Time, everything())
  
  #Pull out background data
  
  BK_DMSOAF <- norm_DMSOAF %>% select(Time, contains('BK'))
  BK_DMSOCARS <- norm_DMSOCARS %>% select(Time, contains('BK'))
  
  BK_melt_DMSO <- reshape2::melt(BK_DMSOCARS, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
  
}
#Full Time Cell Data
{#
  fulltime_norm_DMSOCARS <- as.data.frame(lapply(DMSOCARS[2:ncol(DMSOCARS)], normalize)) 
  
  fulltime_norm_DMSOCARS$Time <- DMSOCARS$Time
  fulltime_norm_DMSOCARS <- fulltime_norm_DMSOCARS %>% select(Time, everything())
  
  #cells only
  fulltime_cells_DMSOCARS <- fulltime_norm_DMSOCARS %>% select(Time, everything(), -contains('BK'))
  
  ftDMSO <- reshape2::melt(fulltime_cells_DMSOCARS, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
  
  #plot of the background data
  ggplot(ftDMSO, aes(Time,Intensity, col=ROI)) + 
    geom_point()
}

#fitting sigmoids
{#Use sigmoid data previously made in AllCPARawData.RData
  #Look at all the plots
  
  DMSOyfit <- list()
  x <- BK_DMSOCARS$Time
  sfDMSO <- as.data.frame(sfDMSO)
  
  for (i in 1:ncol(sfDMSO)) {
    
    DMSOyf <- (sfDMSO[1,i])/(1+exp((sfDMSO[2,i]-x)/sfDMSO[3,i]))
    
    DMSOyfit[[i]] <- DMSOyf
  }
  
  names(DMSOyfit) <- colnames(BK_DMSOCARS[2:ncol(BK_DMSOCARS)])
  DMSOyfit <- data.frame(t(sapply(DMSOyfit,c)))
  DMSOyfit <- as.data.frame(t(DMSOyfit))
  rownames(DMSOyfit) <- NULL
  DMSOyfit$Time <- BK_DMSOCARS$Time
  DMSOyfit <- DMSOyfit %>% select(Time, everything())
  
  # merge by Time the BK ROIs and the fitted lines together into one data frame (this is not necessary)
  mergeROIfit <- merge(DMSOyfit, BK_DMSOCARS, by = "Time")
  
  
  
  meltfit <- reshape2::melt(DMSOyfit, id.vars = "Time", value.name = "Intensity", variable.name = "FitROI")
  
  #plot the fitted lines with the actual data
  DMSOallplots <- ggplot(meltfit, aes(Time,Intensity, col=FitROI)) + 
    geom_line() + geom_point(aes(Time, Intensity, col = ROI), BK_melt_DMSO) +
    theme_classic() + xlab("Time (s)") + ylab("Intensity (a.u.)") +
    theme(panel.grid.major.x=element_line(), legend.position =  "right") +
    guides(color = guide_legend(nrow = 5, byrow = TRUE)) + 
    ggtitle("Unshifted")

  print(DMSOallplots)
}

#Moving all the background data so that the x-mid is at 0 sec
{#
  #Making a new list that contains each ROI with it's own Time column
  DMSOsinglesets <- list()
  
  for (i in 2:ncol(BK_DMSOCARS)) {
    Time <- BK_DMSOCARS$Time
    DMSOy <- BK_DMSOCARS[i]
    DMSOsinglesets[[i]] <- data.frame(Time,DMSOy)
  }
  
  DMSOsinglesets <- compact(DMSOsinglesets)
  names(DMSOsinglesets) <- colnames(BK_DMSOCARS[2:ncol(BK_DMSOCARS)])
  
  #turn the xmid values into a list
  DMSOxmids <- sfDMSO[2,]
  DMSOxmids <- as.list(DMSOxmids)
  
  # use f function to take the "Time" column in the DMSOsinglesets list and subtracts the associated xmid
  DMSOnewtimelist <- lapply(names(DMSOsinglesets), f, some = DMSOsinglesets, other = DMSOxmids)
  names(DMSOnewtimelist) <- colnames(BK_DMSOCARS[2:ncol(BK_DMSOCARS)])
  
  DMSOROIonlylist <- lapply(DMSOsinglesets, function(x) x[!(names(x) %in% c("Time"))])
  
  DMSOcombolist <- mapply(c, DMSOROIonlylist, DMSOnewtimelist, SIMPLIFY = FALSE)
  DMSOcombodf <- compact(DMSOcombolist)
  DMSOcombodf <- dplyr::bind_rows(DMSOcombodf)
  
  # Move Time column to front
  DMSOcombodf <- DMSOcombodf %>%
    select(Time, everything())
  
  meltDMSOcombodf <- reshape2::melt(DMSOcombodf, id.vars = "Time", value.name = "Intensity", variable.name = "ShiftROI")
  
  #Plotting shifted data
  shiftplots <- ggplot(meltDMSOcombodf, aes(Time,Intensity, col=ShiftROI)) + 
    geom_point() +
    theme_classic() + xlab("Time (s)") + ylab("Intensity (a.u.)")  +
    theme(panel.grid.major.x=element_line()) +
    guides(color = guide_legend(ncol = 2, bycol = TRUE))
  
  print(shiftplots)
}

##Shifted fits
{#replaced calculated xmid values with 0
  #create new times data frame to generate lines for sigmoid model
  DMSOnewtimes <- compact(DMSOnewtimelist)
  DMSOnewtimes <- dplyr::bind_rows(DMSOnewtimes)
  
  DMSOyshiftfit <- list()
  x <- DMSOnewtimes$Time
  
  for (i in 1:ncol(sfDMSO)) {
    
    DMSOyshiftf <- (sfDMSO[1,i])/(1+exp((0-x)/sfDMSO[3,i]))
    
    DMSOyshiftfit[[i]] <- DMSOyshiftf
  }
  
  names(DMSOyshiftfit) <- colnames(BK_DMSOCARS[2:ncol(BK_DMSOCARS)])
  DMSOyshiftfit <- data.frame(t(sapply(DMSOyshiftfit,c)))
  DMSOyshiftfit <- as.data.frame(t(DMSOyshiftfit))
  rownames(DMSOyshiftfit) <- NULL
  #colnames(DMSOyfit) <- paste("Fit", colnames(DMSOyfit), sep = ".")
  DMSOyshiftfit$Time <- DMSOnewtimes$Time
  DMSOyshiftfit <- DMSOyshiftfit %>% select(Time, everything())
  
  
  DMSOmeltshiftfit <- reshape2::melt(DMSOyshiftfit, id.vars = "Time", value.name = "Intensity", variable.name = "ShiftFitROI")
  DMSOshiftfitplot <- ggplot(DMSOmeltshiftfit, aes(Time, Intensity, col = ShiftFitROI)) + geom_line()
  
  #plot the fitted lines with the actual data
  DMSOallshiftplots <- ggplot(DMSOmeltshiftfit, aes(Time,Intensity, col=ShiftFitROI)) + 
    geom_line() + geom_point(aes(Time, Intensity, col = ShiftROI), meltDMSOcombodf) +
    theme_classic() + xlab("Time (s)") + ylab("Intensity (a.u.)")  +
    theme(panel.grid.major.x=element_line(), legend.position =  "right") +
    guides(color = guide_legend(nrow = 5, byrow = TRUE)) + 
    ggtitle("Shifted")

  print(DMSOallshiftplots)
}

#put both plots on the same plot
library(ggpubr)
#you have to uninstall the egg package for this to work.
{

  #png(file = "normalizedbackgroundDMSO.png", units="in", width=7, height=5, res=300)
  DMSOcomboplot <- ggarrange(DMSOallplots, DMSOallshiftplots, common.legend = TRUE, legend="bottom") 
  
  print(DMSOcomboplot) 
  
  annotate_figure(comboplot, top = text_grob("Dimethyl Sulfoxide Background Intensity",
                                             color = "black", face = "bold", size = 14))
  #dev.off()
}



### Ethylene Glycol ####

# Subset All Data 
{#
  new_EGAF <- subset(EGAF, Time<800)
  new_EGAF <- subset(new_EGAF, Time>200)
  
  new_EGCARS <- subset(EGCARS, Time<800)
  new_EGCARS <- subset(new_EGCARS, Time>200)
  
  ggplot(data = new_EGCARS, aes(Time, Oct8E2CARSBK1)) + geom_point()
  
  #use lapply to apply it to the columns in the data set - 
  #you don't want to include time column, we will add that back in later
  
  norm_EGAF <- as.data.frame(lapply(new_EGAF[2:ncol(EGAF)], normalize))
  norm_EGCARS <- as.data.frame(lapply(new_EGCARS[2:ncol(EGCARS)], normalize))
  
  #add the time columns back into the data, and reorder them so they are first
  
  norm_EGAF$Time <- new_EGAF$Time
  norm_EGAF <- norm_EGAF %>% select(Time, everything())
  
  norm_EGCARS$Time <- new_EGCARS$Time
  norm_EGCARS <- norm_EGCARS %>% select(Time, everything())
  
  #Pull out background data
  
  BK_EGAF <- norm_EGAF %>% select(Time, contains('BK'))
  BK_EGCARS <- norm_EGCARS %>% select(Time, contains('BK'))
  
  BK_melt_EG <- reshape2::melt(BK_EGCARS, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
  
}
#Full Time Cell Data
{#
  fulltime_norm_EGCARS <- as.data.frame(lapply(EGCARS[2:ncol(EGCARS)], normalize)) 
  
  fulltime_norm_EGCARS$Time <- EGCARS$Time
  fulltime_norm_EGCARS <- fulltime_norm_EGCARS %>% select(Time, everything())
  
  #cells only
  fulltime_cells_EGCARS <- fulltime_norm_EGCARS %>% select(Time, everything(), -contains('BK'))
  
  ftEG <- reshape2::melt(fulltime_cells_EGCARS, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
  
  #plot of the background data
  ggplot(ftEG, aes(Time,Intensity, col=ROI)) + 
    geom_point()
}

#fitting sigmoids
{#Use sigmoid data previously made in AllCPARawData.RData
  #Look at all the plots
  
  EGyfit <- list()
  x <- BK_EGCARS$Time
  sfEG <- as.data.frame(sfEG)
  
  for (i in 1:ncol(sfEG)) {
    
    EGyf <- (sfEG[1,i])/(1+exp((sfEG[2,i]-x)/sfEG[3,i]))
    
    EGyfit[[i]] <- EGyf
  }
  
  names(EGyfit) <- colnames(BK_EGCARS[2:ncol(BK_EGCARS)])
  EGyfit <- data.frame(t(sapply(EGyfit,c)))
  EGyfit <- as.data.frame(t(EGyfit))
  rownames(EGyfit) <- NULL
  EGyfit$Time <- BK_EGCARS$Time
  EGyfit <- EGyfit %>% select(Time, everything())
  
  # merge by Time the BK ROIs and the fitted lines together into one data frame (this is not necessary)
  mergeROIfit <- merge(EGyfit, BK_EGCARS, by = "Time")
  
  
  
  meltfit <- reshape2::melt(EGyfit, id.vars = "Time", value.name = "Intensity", variable.name = "FitROI")
  
  #plot the fitted lines with the actual data
  EGallplots <- ggplot(meltfit, aes(Time,Intensity, col=FitROI)) + 
    geom_line() + geom_point(aes(Time, Intensity, col = ROI), BK_melt_EG) +
    theme_classic() + xlab("Time (s)") + ylab("Intensity (a.u.)")  +
    theme(panel.grid.major.x=element_line(), legend.position =  "right") +
    guides(color = guide_legend(nrow = 5, byrow = TRUE)) + 
    ggtitle("Unshifted")

  
  print(EGallplots)
}

#Moving all the background data so that the x-mid is at 0 sec
{#
  #Making a new list that contains each ROI with it's own Time column
  EGsinglesets <- list()
  
  for (i in 2:ncol(BK_EGCARS)) {
    Time <- BK_EGCARS$Time
    EGy <- BK_EGCARS[i]
    EGsinglesets[[i]] <- data.frame(Time,EGy)
  }
  
  EGsinglesets <- compact(EGsinglesets)
  names(EGsinglesets) <- colnames(BK_EGCARS[2:ncol(BK_EGCARS)])
  
  #turn the xmid values into a list
  EGxmids <- sfEG[2,]
  EGxmids <- as.list(EGxmids)
  
  # use f function to take the "Time" column in the EGsinglesets list and subtracts the associated xmid 
  EGnewtimelist <- lapply(names(EGsinglesets), f, some = EGsinglesets, other = EGxmids)
  names(EGnewtimelist) <- colnames(BK_EGCARS[2:ncol(BK_EGCARS)])
  
  EGROIonlylist <- lapply(EGsinglesets, function(x) x[!(names(x) %in% c("Time"))])
  
  EGcombolist <- mapply(c, EGROIonlylist, EGnewtimelist, SIMPLIFY = FALSE)
  EGcombodf <- compact(EGcombolist)
  EGcombodf <- dplyr::bind_rows(EGcombodf)
  
  # Move Time column to front
  EGcombodf <- EGcombodf %>%
    select(Time, everything())
  
  meltEGcombodf <- reshape2::melt(EGcombodf, id.vars = "Time", value.name = "Intensity", variable.name = "ShiftROI")
  
  #Plotting shifted data
  shiftplots <- ggplot(meltEGcombodf, aes(Time,Intensity, col=ShiftROI)) + geom_point()
  print(shiftplots)
}

##Shifted fits
{#replaced calculated xmid values with 0
  #create new times data frame to generate lines for sigmoid model
  EGnewtimes <- compact(EGnewtimelist)
  EGnewtimes <- dplyr::bind_rows(EGnewtimes)
  
  EGyshiftfit <- list()
  x <- EGnewtimes$Time
  
  for (i in 1:ncol(sfEG)) {
    
    EGyshiftf <- (sfEG[1,i])/(1+exp((0-x)/sfEG[3,i]))
    
    EGyshiftfit[[i]] <- EGyshiftf
  }
  
  names(EGyshiftfit) <- colnames(BK_EGCARS[2:ncol(BK_EGCARS)])
  EGyshiftfit <- data.frame(t(sapply(EGyshiftfit,c)))
  EGyshiftfit <- as.data.frame(t(EGyshiftfit))
  rownames(EGyshiftfit) <- NULL
  #colnames(EGyfit) <- paste("Fit", colnames(EGyfit), sep = ".")
  EGyshiftfit$Time <- EGnewtimes$Time
  EGyshiftfit <- EGyshiftfit %>% select(Time, everything())
  
  
  EGmeltshiftfit <- reshape2::melt(EGyshiftfit, id.vars = "Time", value.name = "Intensity", variable.name = "ShiftFitROI")
  EGshiftfitplot <- ggplot(EGmeltshiftfit, aes(Time, Intensity, col = ShiftFitROI)) + geom_line()
  
  #plot the fitted lines with the actual data
  EGallshiftplots <- ggplot(EGmeltshiftfit, aes(Time,Intensity, col=ShiftFitROI)) + 
    geom_line() + geom_point(aes(Time, Intensity, col = ShiftROI), meltEGcombodf) +
    theme_classic() + xlab("Time (s)") + ylab("Intensity (a.u.)")  +
    theme(panel.grid.major.x=element_line(), legend.position =  "right") +
    guides(color = guide_legend(nrow = 5, byrow = TRUE)) + 
    ggtitle("Shifted")

  
  print(EGallshiftplots)
}

#put both plots on the same plot
{
  
  #png(file = "normalizedbackgroundEG.png", units="in", width=7, height=5, res=300)
 
  EGcomboplot <- ggarrange(EGallplots, EGallshiftplots, common.legend = TRUE, legend="bottom")
  print(EGcomboplot)
  
  annotate_figure(comboplot, top = text_grob("Ethylene Glycol Background Intensity",
                                             color = "black", face = "bold", size = 14))
  #dev.off()

}



### Glycerol ####

# Subset All Data 
{#
new_GlycAF <- subset(GlycAF, Time>200)
new_GlycAF <- subset(GlycAF, Time<800)

new_GlycCARS <- subset(GlycCARS, Time>200)
new_GlycCARS <- subset(new_GlycCARS, Time<800)
  
  ggplot(data = new_GlycCARS, aes(Time, Oct7E5CARSBK1)) + geom_point()
  
  #use lapply to apply it to the columns in the data set - 
  #you don't want to include time column, we will add that back in later
  
  norm_GlycAF <- as.data.frame(lapply(new_GlycAF[2:ncol(GlycAF)], normalize))
  norm_GlycCARS <- as.data.frame(lapply(new_GlycCARS[2:ncol(GlycCARS)], normalize))
  
  #add the time columns back into the data, and reorder them so they are first
  
  norm_GlycAF$Time <- new_GlycAF$Time
  norm_GlycAF <- norm_GlycAF %>% select(Time, everything())
  
  norm_GlycCARS$Time <- new_GlycCARS$Time
  norm_GlycCARS <- norm_GlycCARS %>% select(Time, everything())
  
  #Pull out background data
  
  BK_GlycAF <- norm_GlycAF %>% select(Time, contains('BK'))
  BK_GlycCARS <- norm_GlycCARS %>% select(Time, contains('BK'))
  
  BK_melt_glyc <- reshape2::melt(BK_GlycCARS, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
  
}
#Full Time Cell Data
{#
  fulltime_norm_GlycCARS <- as.data.frame(lapply(GlycCARS[2:ncol(GlycCARS)], normalize)) 
  
  fulltime_norm_GlycCARS$Time <- GlycCARS$Time
  fulltime_norm_GlycCARS <- fulltime_norm_GlycCARS %>% select(Time, everything())
  
  #cells only
  fulltime_cells_GlycCARS <- fulltime_norm_GlycCARS %>% select(Time, everything(), -contains('BK'))
  
  ftglyc <- reshape2::melt(fulltime_cells_GlycCARS, id.vars = "Time", value.name = "Intensity", variable.name = "ROI")
  
  #plot of the background data
  ggplot(ftglyc, aes(Time,Intensity, col=ROI)) + 
    geom_point()
}

#fitting sigmoids
{#Use sigmoid data previously made in AllCPARawData.RData
  #Look at all the plots
  
  glycyfit <- list()
  x <- BK_GlycCARS$Time
  sfGlyc <- as.data.frame(sfGlyc)
  
  for (i in 1:ncol(sfGlyc)) {
    
    glycyf <- (sfGlyc[1,i])/(1+exp((sfGlyc[2,i]-x)/sfGlyc[3,i]))
    
    glycyfit[[i]] <- glycyf
  }
  
  names(glycyfit) <- colnames(BK_GlycCARS[2:ncol(BK_GlycCARS)])
  glycyfit <- data.frame(t(sapply(glycyfit,c)))
  glycyfit <- as.data.frame(t(glycyfit))
  rownames(glycyfit) <- NULL
  glycyfit$Time <- BK_GlycCARS$Time
  glycyfit <- glycyfit %>% select(Time, everything())
  
  # merge by Time the BK ROIs and the fitted lines together into one data frame (this is not necessary)
  mergeROIfit <- merge(glycyfit, BK_GlycCARS, by = "Time")
  
  
  
  meltfit <- reshape2::melt(glycyfit, id.vars = "Time", value.name = "Intensity", variable.name = "FitROI")
  
  #plot the fitted lines with the actual data
  glycallplots <- ggplot(meltfit, aes(Time,Intensity, col=FitROI)) + 
    geom_line() + geom_point(aes(Time, Intensity, col = ROI), BK_melt_glyc)+
    theme_classic() + xlab("Time (s)") + ylab("Intensity (a.u.)")  +
    theme(panel.grid.major.x=element_line(), legend.position =  "right") +
    guides(color = guide_legend(nrow = 5, byrow = TRUE)) + 
    ggtitle("Unshifted")
  
  print(glycallplots)
}

#Moving all the background data so that the x-mid is at 0 sec
{#
  #Making a new list that contains each ROI with it's own Time column
  glycsinglesets <- list()
  
  for (i in 2:ncol(BK_GlycCARS)) {
    Time <- BK_GlycCARS$Time
    glycy <- BK_GlycCARS[i]
    glycsinglesets[[i]] <- data.frame(Time,glycy)
  }
  
  glycsinglesets <- compact(glycsinglesets)
  names(glycsinglesets) <- colnames(BK_GlycCARS[2:ncol(BK_GlycCARS)])
  
  #turn the xmid values into a list
  glycxmids <- sfGlyc[2,]
  glycxmids <- as.list(glycxmids)
  
  # use f function to take the "Time" column in the glycsinglesets list and subtracts the associated xmid 
  glycnewtimelist <- lapply(names(glycsinglesets), f, some = glycsinglesets, other = glycxmids)
  names(glycnewtimelist) <- colnames(BK_GlycCARS[2:ncol(BK_GlycCARS)])
  
  glycROIonlylist <- lapply(glycsinglesets, function(x) x[!(names(x) %in% c("Time"))])
  
  glyccombolist <- mapply(c, glycROIonlylist, glycnewtimelist, SIMPLIFY = FALSE)
  glyccombodf <- compact(glyccombolist)
  glyccombodf <- dplyr::bind_rows(glyccombodf)
  
  # Move Time column to front
  glyccombodf <- glyccombodf %>%
    select(Time, everything())
  
  meltglyccombodf <- reshape2::melt(glyccombodf, id.vars = "Time", value.name = "Intensity", variable.name = "ShiftROI")
  
  #Plotting shifted data
  shiftplots <- ggplot(meltglyccombodf, aes(Time,Intensity, col=ShiftROI)) + geom_point()
  print(shiftplots)
}

##Shifted fits
{#replaced calculated xmid values with 0
  #create new times data frame to generate lines for sigmoid model
  glycnewtimes <- compact(glycnewtimelist)
  glycnewtimes <- dplyr::bind_rows(glycnewtimes)
  
  glycyshiftfit <- list()
  x <- glycnewtimes$Time
  
  for (i in 1:ncol(sfGlyc)) {
    
    glycyshiftf <- (sfGlyc[1,i])/(1+exp((0-x)/sfGlyc[3,i]))
    
    glycyshiftfit[[i]] <- glycyshiftf
  }
  
  names(glycyshiftfit) <- colnames(BK_GlycCARS[2:ncol(BK_GlycCARS)])
  glycyshiftfit <- data.frame(t(sapply(glycyshiftfit,c)))
  glycyshiftfit <- as.data.frame(t(glycyshiftfit))
  rownames(glycyshiftfit) <- NULL
  #colnames(glycyfit) <- paste("Fit", colnames(glycyfit), sep = ".")
  glycyshiftfit$Time <- glycnewtimes$Time
  glycyshiftfit <- glycyshiftfit %>% select(Time, everything())
  
  
  glycmeltshiftfit <- reshape2::melt(glycyshiftfit, id.vars = "Time", value.name = "Intensity", variable.name = "ShiftFitROI")
  glycshiftfitplot <- ggplot(glycmeltshiftfit, aes(Time, Intensity, col = ShiftFitROI)) + geom_line()
  
  #plot the fitted lines with the actual data
  glycallshiftplots <- ggplot(glycmeltshiftfit, aes(Time,Intensity, col=ShiftFitROI)) + 
    geom_line() + geom_point(aes(Time, Intensity, col = ShiftROI), meltglyccombodf) +
    theme_classic() + xlab("Time (s)") + ylab("Intensity (a.u.)")  +
    theme(panel.grid.major.x=element_line(), legend.position =  "right") +
    guides(color = guide_legend(nrow = 5, byrow = TRUE)) + 
    ggtitle("Shifted")
  
  print(glycallshiftplots)
}

#put both plots on the same plot
{
  
  #png(file = "normalizedbackgroundglyc.png", units="in", width=7, height=5, res=300)
  
  glyccomboplot <- ggarrange(glycallplots, glycallshiftplots, common.legend = TRUE, legend="bottom")
  print(glyccomboplot)
  
  annotate_figure(comboplot, top = text_grob("Glycerol Background Intensity",
                                             color = "black", face = "bold", size = 14))
  #dev.off()

  
}



