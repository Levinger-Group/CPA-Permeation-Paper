#Modelling the response of cells based on exposure
#this is also called "allResponseCountingData.R" elsewhere
library(purrr)
library(tidyr)
library(dplyr)
library(readxl)
library(aod)

library(plotly)
library(RColorBrewer)
library(dplyr)

library(gtsummary)
library(gridExtra)


# data that is not binary ####
#importing data 
read_allsheets_excel <- function(filename, tibble = TRUE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

fullcounts <- rbind.data.frame(read_allsheets_excel(file.choose()))

#selecting just the response columns
#This is based off my own sheets and allowed me to choose the columns that were important
#I then re-named them to simpler names. 
responsecounts <- subset(fullcounts, select = c(2,3,4,5,6,7))
names(responsecounts)[1] <- "observable"
names(responsecounts)[2] <- "dead"
names(responsecounts)[3] <- "response"
names(responsecounts)[4] <- "Plas"
names(responsecounts)[5] <- "Deplas"
names(responsecounts)[6] <- "Treatment"

#here I added in a new column called "alive" by subtracting the number of dead cells from observable cells
responsecounts$alive <- responsecounts$observable - responsecounts$dead
responsecounts <- responsecounts[responsecounts$alive != 0, ]
alive_responsecounts <- subset(responsecounts, select = c(3,4,5,6,7))

# expanding the data into binary counts from: https://aosmith.rbind.io/2019/10/04/expanding-binomial-to-binary/
#in this data, 1 means the response occurred, 0 means it did not occur
binarycounts = pmap_dfr(alive_responsecounts, 
                        function(Treatment, alive, response, Plas, Deplas) {
                          data.frame(Treatment = Treatment,
                                     response = c( rep(1, response),
                                                   rep(0, alive - response) ),
                                     Plas = c( rep(1, Plas),
                                               rep(0, (alive) - Plas)),
                                     Deplas = c( rep(1, Deplas),
                                                 rep(0, alive - Deplas)))
                        }
)


# Alternate for already binary Data ####
#If you already have binary data in your excel sheet, use this instead of the above section.
bidata <- read_excel(file.choose()) # used Counting Cells_Final
bidata <- bidata[bidata$visible != 0, ]
bidata <- bidata[,c(3,4,5,6,7)]

# Actual model ####
# check that these new data match the original excel data
#This spits out the number of each identifier (response, plas, deplas)
#you can check the numbers against what you recorded in your excel sheet
xtabs(~Treatment + response, data = bidata)
xtabs(~Treatment + Plas, data = bidata)
xtabs(~Treatment + Deplas, data = bidata)

# make the models based on https://stats.idre.ucla.edu/r/dae/logit-regression/
bidata$Treatment <- factor(bidata$Treatment)
PredictResponse <- glm(response ~ Treatment, bidata[bidata$visible == 1,], family = "binomial")
#Make the model into a group of matrices
betaResponse <- summary(PredictResponse)
coefficientsResponse <- betaResponse[["coefficients"]]


#These are the predictors for plasmolysis and deplasmolysis. 
#They're set up to only include the the plasmolysis or deplasmolysis cells associated with response or plasmolysis respectively

PredictPlas <- glm(Plas ~ Treatment, data = bidata[bidata$visible == 1,], family = "binomial")
betaPlas <- summary(PredictPlas)
coefficientsPlas <- betaPlas[["coefficients"]]

PredictDeplas <- glm(Deplas ~ Treatment, data =bidata[bidata$Plas == 1,], family = "binomial")
betaDeplas <- summary(PredictDeplas)
coefficientsDeplas <- betaDeplas[["coefficients"]]

#physical notes on turning log odds into regular probability in lab notebook.
#Generally speaking a log odd is log[P/(1-P)] where P is the probability of success and 1-P is the probability of failure
# in Logit analysis, y = b_0 + b_i * X_i where, in this case, X_i is the treatment
# for this analysis, y = b_0 + b_1 * (TreatmentEG) + b_2 * (Treatmentglcycerol)
# b is the predicted log odds of each treatment, with b_0 being the predicted log odds of the DMSO treatment
# to interpret it more easily, change log odds to probability by doing exp{b_0 + b_i}/(1+exp{b_0 + b_i})

#pull out the beta values for each Treatment and each model. b_0 is for DMSo, b_1 is for EG, and b_2 is glycerol
b0Response <- coefficientsResponse[1,1]
b1Response <- coefficientsResponse[2,1]
b2Response <- coefficientsResponse[3,1]

b0Plas <- coefficientsPlas[1,1]
b1Plas <- coefficientsPlas[2,1]
b2Plas <- coefficientsPlas[3,1]

b0Deplas <- coefficientsDeplas[1,1]
b1Deplas <- coefficientsDeplas[2,1]
b2Deplas <- coefficientsDeplas[3,1]

#determine probability for each treatment and each model
ResponseDMSO <- exp(b0Response)/(1+exp(b0Response))
ResponseEG <- exp(b1Response + b0Response)/(1+exp(b1Response + b0Response))
ResponseGlycerol <- exp(b2Response + b0Response)/(1+exp(b2Response + b0Response))

PlasDMSO <- exp(b0Plas)/(1+exp(b0Plas))
PlasEG <- exp(b0Plas + b1Plas)/(1 +  exp(b0Plas + b1Plas))
PlasGlycerol <- exp(b0Plas + b2Plas)/(1 + exp(b0Plas + b2Plas))

DeplasDMSO <- exp(b0Deplas)/(1 + exp(b0Deplas))
DeplasEG <- exp(b0Deplas + b1Deplas)/(1 + exp(b0Deplas + b1Deplas))
DeplasGlycerol <- exp(b0Deplas + b2Deplas)/(1 + exp(b0Deplas + b2Deplas))

# Make these into a table of probabilities
table <- matrix(c(ResponseDMSO, PlasDMSO, DeplasDMSO, 
                  ResponseEG, PlasEG, DeplasEG, 
                  ResponseGlycerol, PlasGlycerol, DeplasGlycerol)*100, 
                ncol = 3, byrow = TRUE)
colnames(table) <- c("Response", "Plasmolysis", "Deplasmolysis")
rownames(table) <- c("DMSO", "Ethylene Glycol", "Glycerol")
table1 <- as.data.frame(table) %>% 
  round(1)

#Print results
table1

mypath <- file.path("D:", "Writing", "Cell response Paper", "Files for Sharing", "Bright Field Data", fsep="/",
                    paste("CellResponseModel.csv"))
write.csv(table1,mypath, row.names = TRUE)

# Goodness of Fit ####
library(ResourceSelection)
hoslem.test(bidata$response, fitted(PredictResponse))
hoslem.test(bidata$Plas, fitted(PredictPlas))
hoslem.test(bidata$Deplas, fitted(PredictDeplas))
# You want the p-value of the goodness of fit to be >0.5 to show there is no 
# difference between model and actual data