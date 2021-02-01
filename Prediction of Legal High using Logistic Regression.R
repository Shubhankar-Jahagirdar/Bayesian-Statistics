# MATH2269 Module 7 Codes

graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
library(ggplot2)
library(ggpubr)
library(ks)
library(rjags)
library(runjags)
source("DBDA2E-utilities.R") 
#setwd("~/Documents/MATH2269_Bayesian/2020/presentations/Module 7/Application2/")


#===============PRELIMINARY FUNCTIONS FOR POSTERIOR INFERENCES====================

smryMCMC_HD = function(  codaSamples , compVal = NULL,  saveName=NULL) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  paramName = colnames(mcmcMat)
  for ( pName in paramName ) {
    if (pName %in% colnames(compVal)){
      if (!is.na(compVal[pName])) {
        summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] , 
                                                          compVal = as.numeric(compVal[pName]) ))
      }
      else {
        summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] ) )
      }
    } else {
      summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] ) )
    }
  }
  rownames(summaryInfo) = paramName
  
  # summaryInfo = rbind( summaryInfo , 
  #                      "tau" = summarizePost( mcmcMat[,"tau"] ) )
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

plotMCMC_HD = function( codaSamples , data , xName="x" , yName="y", preds = FALSE ,
                     showCurve=FALSE ,  pairsPlot=FALSE , compVal = NULL, 
                     saveName=NULL , saveType="jpg" ) {
  # showCurve is TRUE or FALSE and indicates whether the posterior should
  #   be displayed as a histogram (by default) or by an approximate curve.
  # pairsPlot is TRUE or FALSE and indicates whether scatterplots of pairs
  #   of parameters should be displayed.
  #-----------------------------------------------------------------------------
  y = data[,yName]
  x = as.matrix(data[,xName])
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  # zbeta0 = mcmcMat[,"zbeta0"]
  # zbeta  = mcmcMat[,grep("^zbeta$|^zbeta\\[",colnames(mcmcMat))]
  # if ( ncol(x)==1 ) { zbeta = matrix( zbeta , ncol=1 ) }
  beta0 = mcmcMat[,"beta0"]
  beta  = mcmcMat[,grep("^beta$|^beta\\[",colnames(mcmcMat))]
  if ( ncol(x)==1 ) { beta = matrix( beta , ncol=1 ) }
  if (preds){
    pred = mcmcMat[,grep("^pred$|^pred\\[",colnames(mcmcMat))]
  } # Added by Demirhan
  guess = mcmcMat[,"guess"]
  #-----------------------------------------------------------------------------
  # Compute R^2 for credible parameters:
  YcorX = cor( y , x ) # correlation of y with each x predictor
  Rsq = beta %*% matrix( YcorX , ncol=1 )
  #-----------------------------------------------------------------------------
  if ( pairsPlot ) {
    # Plot the parameters pairwise, to see correlations:
    openGraph()
    nPtToPlot = 1000
    plotIdx = floor(seq(1,chainLength,by=chainLength/nPtToPlot))
    panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      usr = par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
    }
    pairs( cbind( beta0 , beta , tau )[plotIdx,] ,
           labels=c( "beta[0]" , 
                     paste0("beta[",1:ncol(beta),"]\n",xName) , 
                     expression(tau) ) , 
           lower.panel=panel.cor , col="skyblue" )
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"PostPairs",sep=""), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
  # Marginal histograms:
  
  decideOpenGraph = function( panelCount , saveName , finished=FALSE , 
                              nRow=2 , nCol=3 ) {
    # If finishing a set:
    if ( finished==TRUE ) {
      if ( !is.null(saveName) ) {
        saveGraph( file=paste0(saveName,ceiling((panelCount-1)/(nRow*nCol))), 
                   type=saveType)
      }
      panelCount = 1 # re-set panelCount
      return(panelCount)
    } else {
      # If this is first panel of a graph:
      if ( ( panelCount %% (nRow*nCol) ) == 1 ) {
        # If previous graph was open, save previous one:
        if ( panelCount>1 & !is.null(saveName) ) {
          saveGraph( file=paste0(saveName,(panelCount%/%(nRow*nCol))), 
                     type=saveType)
        }
        # Open new graph
        openGraph(width=nCol*7.0/3,height=nRow*2.0)
        layout( matrix( 1:(nRow*nCol) , nrow=nRow, byrow=TRUE ) )
        par( mar=c(4,4,2.5,0.5) , mgp=c(2.5,0.7,0) )
      }
      # Increment and return panel count:
      panelCount = panelCount+1
      return(panelCount)
    }
  }
  
  # Original scale:
  panelCount = 1
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( beta0 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(beta[0]) , main="Intercept", compVal = as.numeric(compVal["beta0"] ))
  for ( bIdx in 1:ncol(beta) ) {
    panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
    if (!is.na(compVal[paste0("beta[",bIdx,"]")])){
      histInfo = plotPost( beta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                           xlab=bquote(beta[.(bIdx)]) , main=xName[bIdx],
                           compVal = as.numeric(compVal[paste0("beta[",bIdx,"]")]))
    } else{
      histInfo = plotPost( beta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                           xlab=bquote(beta[.(bIdx)]) , main=xName[bIdx])
    }
  }
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( Rsq , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(R^2) , main=paste("Prop Var Accntd") , finished=FALSE )
  
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( guess , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab="Guess parameter" , main=paste("Prop Var Accntd") , finished=TRUE )
  
  panelCount = 1
  if ( pred){
    
    for ( pIdx in 1:ncol(pred) ) {
      panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
      histInfo = plotPost( pred[,pIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(pred[.(pIdx)]) , main=paste0("Prediction ",pIdx) ) 
    }
  }# Added by Demirhan
  # Standardized scale:
  panelCount = 1
  # panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
  # histInfo = plotPost( zbeta0 , cex.lab = 1.75 , showCurve=showCurve ,
  #                      xlab=bquote(z*beta[0]) , main="Intercept" )
  # for ( bIdx in 1:ncol(beta) ) {
  #   panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
  #   histInfo = plotPost( zbeta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
  #                        xlab=bquote(z*beta[.(bIdx)]) , main=xName[bIdx] )
  # }
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
  histInfo = plotPost( Rsq , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(R^2) , main=paste("Prop Var Accntd") )
  panelCount = decideOpenGraph( panelCount , finished=TRUE , saveName=paste0(saveName,"PostMargZ") )
  
  #-----------------------------------------------------------------------------
}

#===============PRELIMINARY FUNCTIONS FOR POSTERIOR INFERENCES====================


myData1 <- read.csv("DrugConsumptionFinal.csv")
head(myData1)
#myData1 <- myData1[,c(1,2,3,5,6,7,8,10)]
myData1$Gender <- as.numeric(as.factor(myData1$Gender)) - 1 # To get 0/1 instead of 1/2; Female = 0; Male = 1
head(myData1)
attach(myData1)
smp_siz = floor(0.75*nrow(myData1))
set.seed(123)
train_sample = sample(seq_len(nrow(myData1)),size = smp_siz)
trainData = myData1[train_sample,]
testData = myData1[-train_sample,]

#summary(myData2)

# Drop NAs! 
# myData2 = myData2[which(!is.na(myData2[,'Age'])),]

# Relace NAs with mean+error
#meanAge <- mean(myData2[,'Age'], na.rm = TRUE)
#myData2[which(is.na(myData2[,'Age'])),'Age'] <- rnorm(length(which(is.na(myData2[,'Age']))),meanAge,1)

#=== 1. Descriptive look ===
# Scatter plots
p1 <- ggplot(trainData, aes(x=X18.24, y = Legalh)) +
  geom_point()

p2 <- ggplot(trainData, aes(x=X25.34, y = Legalh)) +
  geom_point()
p3 <- ggplot(trainData, aes(x=X35.44, y = Legalh)) +
  geom_point()
p4 <- ggplot(trainData, aes(x=X45.54, y = Legalh)) +
  geom_point()
p5 <- ggplot(trainData, aes(x=X55.64, y = Legalh)) +
  geom_point()
p6 <- ggplot(trainData, aes(x=Gender, y = Legalh)) +
  geom_point()
p7 <- ggplot(trainData, aes(x=Nscore, y = Legalh)) +
  geom_point()
p8 <- ggplot(trainData, aes(x=Escore, y = Legalh)) +
  geom_point()
p9 <- ggplot(trainData, aes(x=Oscore, y = Legalh)) +
  geom_point()
p10 <- ggplot(trainData, aes(x=Ascore, y = Legalh)) +
  geom_point()
p11 <- ggplot(trainData, aes(x=Cscore, y = Legalh)) +
  geom_point()
p12 <- ggplot(trainData, aes(x=Impulsiveness, y = Legalh)) +
  geom_point()
p13 <- ggplot(trainData, aes(x=Sensation_seeking, y = Legalh)) +
  geom_point()

figure <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, nrow = 3, ncol = 2)
figure

#table(myData2[,c(2,3)])
#table(myData2[,c(2,4)])



# THE DATA.
y = trainData[,"Legalh"]
x = as.matrix(trainData[,2:14])

PredData = testData
head(PredData)
#PredData = PredData[,c(1,2,4,5,6,7,9)]
#PredData$Sex <- as.numeric(as.factor(PredData$Sex)) - 1
summary(PredData)

#PredData[which(is.na(PredData[,'Fare'])),'Fare'] <- mean(PredData[,'Fare'],na.rm = TRUE)
#PredData[which(is.na(PredData[,'Age'])),'Age'] <- mean(PredData[,'Age'],na.rm = TRUE)

# PredData = PredData[which(!is.na(PredData[,'Age'])),]
# PredData = PredData[which(!is.na(PredData[,'Fare'])),]
xPred = as.matrix(PredData[-1])

Nx = ncol(x)

# Specify the data in a list, for later shipment to JAGS:
dataList <- list(
  x = x ,
  y = y ,
  xPred = xPred ,
  Ntotal = length(y),
  Nx = Nx, 
  Npred = nrow(xPred)
)

# First run without initials!
initsList <- list(
  beta0 = 0,
  beta1 = 0,
  beta2 = 0,
  beta3 = 0,
  beta4 = 0,
  beta5 = 0,
  beta6 = 0,
  beta7 = 0,
  beta8 = 0,
  beta9 = 0,
  beta10 = 0,
  beta11 = 0,
  beta12 = 0,
  beta13 = 0
)

modelString = "
# data {
#   for ( j in 1:Nx ) {
#     xm[j]  <- mean(x[,j])
#     xsd[j] <-   sd(x[,j])
#     for ( i in 1:Ntotal ) {
#       zx[i,j] <- ( x[i,j] - xm[j] ) / xsd[j]
#     }
#   } 
# }

model {
  for ( i in 1:Ntotal ) {
    # In JAGS, ilogit is logistic:
    y[i] ~ dbern( mu[i] )
      mu[i] <- ( guess*(1/2) + (1.0-guess)*ilogit(beta0+sum(beta[1:Nx]*x[i,1:Nx])) )
  }
  # Priors vague on standardized scale:
  beta0 ~ dnorm( 0 , 1/2^2 )
  # non-informative run
   for ( j in 1:Nx ) {
     beta[j] ~ dnorm( 0 , 1/2^2 )
   }
  # Informative run
  #beta[1] ~ dnorm( 0, 1/4 ) #Pclass
  #beta[2] ~ dnorm( 0 , 1/12 ) #Sex
  #beta[3] ~ dnorm( 0 , 1/14 ) #Age
  #beta[4] ~ dnorm( 0 , 1/2^2 ) #SibSp
  #beta[5] ~ dnorm( 0 , 1/2^2 ) #Parch
  #beta[6] ~ dnorm( 0 , 1/14 ) #Fare
  guess ~ dbeta(1,9)
  # Transform to original scale:
  # beta[1:Nx] <- zbeta[1:Nx] / xsd[1:Nx]
  # beta0 <- zbeta0 - sum( zbeta[1:Nx] * xm[1:Nx] / xsd[1:Nx] )

  # Compute predictions at every step of the MCMC
  for ( k in 1:Npred){
    pred[k] <- ilogit(beta0 + sum(beta[1:Nx] * xPred[k,1:Nx]))
  }
}
" # close quote for modelString
# Write out modelString to a text file
writeLines( modelString , con="TEMPmodel.txt" )

# parameters = c( "zbeta0" , "beta0")
parameters = c( "beta0")
for ( i in 1:Nx){
  # parameters = c(parameters, paste0("zbeta[",i,"]"), paste0("beta[",i,"]")) 
  parameters = c(parameters, paste0("beta[",i,"]")) 
}
for ( i in 1:nrow(xPred)){
  parameters = c(parameters, paste0("pred[",i,"]")) 
}

parameters = c(parameters, "guess")
adaptSteps = 100  # Number of steps to "tune" the samplers
burnInSteps = 300
nChains = 3
thinSteps = 7
numSavedSteps = 100
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

runJagsOut <- run.jags( method="parallel" ,
                        model="TEMPmodel.txt" ,
                        monitor=parameters  ,
                        data=dataList ,
                        # inits=initsList ,
                        n.chains=nChains ,
                        adapt=adaptSteps ,
                        burnin=burnInSteps ,
                        sample=numSavedSteps ,
                        thin=thinSteps , summarise=FALSE , plots=FALSE )
codaSamples = as.mcmc.list( runJagsOut )

load(file="rEnvironment_80_20_1500_2000_3_120_1000.RData")

diagMCMC( codaSamples , parName="beta0" )
for ( i in 1:Nx){
  diagMCMC( codaSamples , parName=paste0("beta[",i,"]") )
}

# diagMCMC( codaSamples , parName="zbeta0" )
# for ( i in 1:Nx){
#   diagMCMC( codaSamples , parName=paste0("zbeta[",i,"]") )
# }
# for ( i in 1:nrow(xPred)){
#   diagMCMC( codaSamples , parName=paste0("pred[",i,"]") )
# }

compVal <- data.frame("beta0" = 15, "beta[1]" = 0, "beta[2]" = 0, "beta[3]" = 0, "beta[4]" =  0,  "beta[5]" =  0, 
                      "beta[6]" =  0,"beta[7]" = 0,"beta[8]" = 0, "beta[9]" = 0, "beta[10]" = 0, "beta[11]" = 0, "beta[12]" = 0,"beta[13]" = 0, check.names=FALSE)

summaryInfo <- smryMCMC_HD( codaSamples = codaSamples , compVal = compVal )
print(summaryInfo)

plotMCMC_HD( codaSamples = codaSamples , data = trainData, xName=c("X18.24","X25.34","X35.44","X45.54","X55.64","Gender","Nscore","Escore","Oscore","Ascore","Cscore","Impulsiveness","Sensation_seeking"),
          yName="Legalh", compVal = compVal, preds = FALSE)

# Parch and Fare are insignificant in non-informative run.
# Now let's run it with informative priors.
# With the informative run we find all the variabels significant. 


# Predictions for full records in training set
preds <- data.frame(PassengerId = PredData[,1], PredProb = summaryInfo[9:426,3], Survival = PredData[,1] )

threshold <- 0.5# summaryInfo[427,3]
preds[which(preds[,2]<threshold),3] <- 0
preds[which(preds[,2]>threshold),3] <- 1

predsSorted <- preds[order(preds$PassengerId),]
# write.csv(predsSorted,"predsSorted.csv")
table(preds$Survival)

actualSurvAll <- read.csv("actualSurvivalStatus.csv")
actualSurvival <- actualSurvAll[which(actualSurvAll$PassengerId %in% predsSorted$PassengerId),2]


# ============ Predictive check ============

confusionMatrix <- function(resp, pred){
  classRes <- data.frame(response = resp , predicted = pred)
  conf = xtabs(~ predicted + response, data = classRes)
  
  accuracy = sum(diag(conf))/sum(conf)
  accuracy
  precision = conf[1,1]/(conf[1,1]+conf[1,2])
  precision
  recall = conf[1,1]/(conf[1,1]+conf[2,1])
  recall
  Fscore = 2*((precision*recall)/(precision+recall))
  Fscore
  return(list(accuracy = accuracy, precision = precision, recall = recall, Fscore = Fscore))
}

confusionMatrix(resp = actualSurvival, pred = predsSorted[,3])


