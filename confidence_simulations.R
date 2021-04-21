#Category Confidence - simulations
#looking to find values of variance that maximise differences in confidence
#area of interest - between the means [200-500]

#load in packages
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(stringr)

####################################################
#start with just 1 estimated simulation
#use the values that were run in the experiment itself 70 and 100 sd
set.seed(1234)

#simulate distribution 1
narrow_mean <- 200
narrow_sd <- 70 
#simulate distribution 2
wide_mean <- 500
wide_sd <- 100 

#for n trials
n_trials <- 10000

#pick a distribution randomly
df <- data.frame(
  distribution=factor(rep(c("narrow", "wide"), each=n_trials)),
  location=round(c(rnorm(n_trials, mean=narrow_mean, sd=narrow_sd),
                   rnorm(n_trials, mean=wide_mean, sd=wide_sd)))
)

#resample out of range values ~~~~~~~~ DO THIS 
df <- subset(df, location<801 & location>0) #this just splices the distributions

#plot the distributions
a <- ggplot(data=df, aes(location)) +
  geom_density(aes(fill = distribution), alpha = 0.4)
a

#find the likelihood of each distribution at each location (5 pixel blocks of location)
#set these blocks around the midline of these 
f <- function(x) dnorm(x, m=narrow_mean, sd=narrow_sd)  - dnorm(x, m=wide_mean, sd=wide_sd)
midline <- uniroot(f, interval=c(narrow_mean, wide_mean))
midline <- round(midline$root)

locations <- seq.int(0, 795, 5)
likelihood <- matrix(data = 0, nrow = (length(locations)), ncol = 5) 
for (val in locations){
  val5 <- val +5
  ofInterest <- subset(df, location>val-1 & location<val5+1)
  ofInterestNarrow <- subset(df, location>val-1 & location<val5+1 & distribution == "narrow")
  totalLength <- length(ofInterest$location)
  lengthNarrow <- length(ofInterestNarrow$location)
  likelihood[((val+5)/5),1] <- signif(lengthNarrow/totalLength, digits = 5) #narrow
  likelihood[((val+5)/5),2] <- signif(1 - lengthNarrow/totalLength, digits = 5) #wide
}
likelihood[is.nan(likelihood)] <- 0 #replace any NaN value with 0 for computation

#make choice based on likelihood at that point 
locationBin <- c(1:160)

for (i in locationBin){
  if(likelihood[i,1]>likelihood[i,2]){
    likelihood[i,3] <- "narrow"
  }else if (likelihood[i,1]==likelihood[i,2]){
    if (i < 80){
      likelihood[i,3] <- "narrow"
    } else{
      likelihood[i,3] <- "wide"
    }
  }else{
    likelihood[i,3] <- "wide"
  }
  #make confidence report equal to the choice likelihood at that location
  if (likelihood[i,3] == "narrow"){
    likelihood[i,4] <- likelihood[i,1]#set confidence
  }else{
    likelihood[i,4] <- likelihood[i,2]#set confidence
  }
}

#look at difference in confidence between means of each distribution
likelihood[,5] <-seq.int(0, 795, 5) #set location bin value

dataframe <- as.data.frame(likelihood)
dataframe <- dataframe %>% 
  rename(
    narrowLikeli = V1,
    wideLikeli = V2,
    choice = V3,
    confidence = V4,
    location = V5
  )
dataframe$confidence <- signif(as.numeric(dataframe$confidence), digits = 6)
dataframe$location <- as.numeric(dataframe$location)
dataframe$conf <- dataframe$confidence / (as.numeric(dataframe$narrowLikeli) + as.numeric(dataframe$wideLikeli))
ggplot(data = dataframe, mapping = aes(location, conf))+ geom_point()

#find the difference in confidence reports between distributions at each step from the decision bound
midData <- subset(dataframe, location < 600 & location > 100) #pull out central values of interest
midBin <- which.min(abs(midData$location-midline)) #find the clostest 5pix bin to the midline
distance <- c(1:30)
confDiffsCurrent <- matrix(data = NA, nrow = length(distance), ncol = 1)
for(val in distance){
  narrowSide <-  subset(dataframe, location == (midData$location[midBin] - (val*5)))
  narrowConf <- narrowSide$confidence /(as.numeric(narrowSide$narrowLikeli) + as.numeric(narrowSide$wideLikeli))
  
  wideSide <- subset(dataframe, location == (midData$location[midBin] + ((val)*5)))
  wideConf <- wideSide$confidence /(as.numeric(wideSide$narrowLikeli) + as.numeric(wideSide$wideLikeli))
  
  confDiffsCurrent[val,1] <- abs(narrowConf- wideConf)
}

#visualise these differences
x <- confDiffsCurrent[,1]
y <- seq.int(0,145,5)

plot(y,x, xlab = "Distance from bound", ylab = "Difference in confidence", main = "Current Variance Pair") #going from midline to means (the differences drop off heading to the means)

# Plot histogram of x
hist(x,
     main = "Current Variance Pair",
     xlab = "Confidence Difference",
     xlim = c(0, 1))

###################################
#Now loop through many variance options and see which one is most different for the distributions 
#keeping sds within 0-800, so 2sd from the mean should still be inside the margins
narrow_sd <- seq.int(10, 80, 10) #narrow sd from 10-80
wide_sd <- seq.int(80, 150, 10) # wide sd from 50-120
variances <- c(1:64)
locations <- c(1:800)
distance <- c(1:30)
n_trials <- 10000

confidence <- matrix(data = 0, nrow = (length(locations)), ncol = length(narrow_sd)*length(wide_sd)) 
midline <- matrix(data = 0, nrow = 1, ncol = length(narrow_sd)*length(wide_sd))
confDiffs <- matrix(data = 0, nrow = length(distance), ncol = length(variances))

for(k in narrow_sd){
  for(j in wide_sd){
    #pick the distributions randomly
    set.seed(1234)
    df <- data.frame(
      distribution=factor(rep(c("narrow", "wide"), each=n_trials)),
      location=round(c(rnorm(n_trials, mean=narrow_mean, sd=k),
                       rnorm(n_trials, mean=wide_mean, sd=j)))
    )
    df <- subset(df, location<801 & location>0) #this just splices the distributions within screen frame
    
    index1 <- which(narrow_sd == k)
    index2 <- which(wide_sd == j)
    index <- ((8*index1)-8) + index2
    
    #find the midline location point of the distribution
    f <- function(x) dnorm(x, m=narrow_mean, sd=k)  - dnorm(x, m=wide_mean, sd=j)
    mid <- uniroot(f, interval=c(narrow_mean, wide_mean))
    midline[,index] <- round(mid$root)
    
    #find the likelihood of each distribution at each location (try every step)
    likelihood <-matrix(data = 0, nrow = (length(locations)), ncol = 5)
    for (val in locations){
      ofInterest <- subset(df, location>val-1 & location<val+1)
      ofInterestNarrow <- subset(df, location>val-1 & location<val+1 & distribution == "narrow")
      totalLength <- length(ofInterest$location)
      lengthNarrow <- length(ofInterestNarrow$location)
      likelihood[val,1] <- signif(lengthNarrow/totalLength, digits = 5) #narrow
      likelihood[val,2] <- signif(1 - lengthNarrow/totalLength, digits = 5) #wide
    }
    #replace any NaN value with 0 for computation
    likelihood[is.nan(likelihood)] <- 0
    
    #make choice based on likelihood at that point
    for (i in locations){
      if(likelihood[i,1]>likelihood[i,2]){
        likelihood[i,3] <- "narrow"
      }else if (likelihood[i,1]==likelihood[i,2]){
        if (i < 80){
          likelihood[i,3] <- "narrow"
        } else{
          likelihood[i,3] <- "wide"
        }
      }else{
        likelihood[i,3] <- "wide"
      }
      #make confidence report equal to the likelihood at that location
      if (likelihood[i,3] == "narrow"){
        likelihood[i,4] <- likelihood[i,1]#set confidence
      }else{
        likelihood[i,4] <- likelihood[i,2]#set confidence
      }
    }
    
    #look at difference in confidence between means of each distribution
    likelihood[,5] <- locations #set location bin value
    
    dataframe <- as.data.frame(likelihood)
    dataframe <- dataframe %>%
      rename(
        narrowLikeli = V1,
        wideLikeli = V2,
        choice = V3,
        confidence = V4,
        location = V5
      )
    dataframe$confidence <- signif(as.numeric(dataframe$confidence), digits = 6)
    dataframe$location <- as.numeric(dataframe$location)
    confidence[,index] <- dataframe$confidence
    
    #find the difference in confidence reports between distributions at each step from the decision bound
    for(dist in distance){
      narrowSide <-  subset(dataframe, location == (midline[,index] - (dist*5)))
      narrowConf <- narrowSide$confidence /(as.numeric(dataframe$narrowLikeli) + as.numeric(narrowSide$wideLikeli))
      
      wideSide <- subset(dataframe, location == (midline[,index] + ((dist)*5)))
      wideConf <- wideSide$confidence /(as.numeric(wideSide$narrowLikeli) + as.numeric(wideSide$wideLikeli))
      
      confDiffs[dist,index] <- abs(narrowConf- wideConf)
    }
  }
}


#visualise these differences
for (i in variances) {
  x <- confDiffs[,i]
  y <- seq.int(0,145,5)
  plot(y,x, xlab = "Distance from bound", ylab = "Difference in confidence", main = paste("Variance Pair", i), ylim =c(0, .5)) 
  #going from midline to means (the differences drop off heading to the means)
}
#visual inspection shows some pairs of variances values are obviously poor, as extremes 0-1 with no varation 
#(so either the same or completely opposite confidence reports i.e. no overlap of distributions)


#identify the columns with range of 0-1 in confidence differences (extremes with no intermediate variations)
confVals <- matrix(data = NA, nrow = length(distance), ncol = length(variances))
confDiffs[is.na(confDiffs)] <- 0
for (i in variances){
  x <- confDiffs[,i]
  if (sum(range(x)) != 1){
    confVals[,i] <- x
}
}

#remove columns of NA
not_all_na <- function(x) any(!is.na(x))
confVals <- as.data.frame(confVals)
confVals <- confVals %>% select(where(not_all_na))

confVals <- as.matrix(confVals)

#find the column with the highest median diff in conf values (mean gives same answer)
medians <- c(1:(length(confVals[1,])-1))
confMed <- matrix(data = NA, nrow = 1, ncol = length(medians))
for (i in medians) {
  x <- confVals[,i]
  confMed[i] <- median(x)
}
index <- which.max(confMed)#use this index to trace back
colName <-colnames(confVals)[index]
colName <- as.numeric(sub('.', '', colName)) #this is the numeric index for the variance pair we need to recreate

narrowIndex <- round(colName/length(narrow_sd))
narrowVariance <- narrow_sd[narrowIndex]
wideIndex <- length(wide_sd)*(((colName/length(narrow_sd)) -(narrowIndex)+1))
wideVariance <- wide_sd[wideIndex]

optimalVariance <- c(narrowVariance, wideVariance)
print(optimalVariance)
#if take median optimals are at (70,100) - if take median (80,110)
####################################################
#simulate with that set of values - see what the optimal really looks like
#distributions means and spreads must exist between 0-800 units
set.seed(1234)

#simulate distributions
narrow_sd <- optimalVariance[1] 
wide_sd <- optimalVariance[2] 

#pick from distribution randomly
df <- data.frame(
  distribution=factor(rep(c("narrow", "wide"), each=n_trials)),
  location=round(c(rnorm(n_trials, mean=narrow_mean, sd=narrow_sd),
                   rnorm(n_trials, mean=wide_mean, sd=wide_sd)))
)

#resample out of range values ~~~~~~~~ DO THIS 
df <- subset(df, location<801 & location>0) #this just splices the distributions

#plot the distributions
a <- ggplot(data=df, aes(location)) +
  geom_density(aes(fill = distribution), alpha = 0.4)
a

#find the likelihood of each distribution at each location (5 pixel blocks of location)
locations <- seq.int(0, 795, 5)
likelihood <- matrix(data = 0, nrow = (length(locations)), ncol = 5) 
for (val in locations){
  val5 <- val +5
  ofInterest <- subset(df, location>val-1 & location<val5+1)
  ofInterestNarrow <- subset(df, location>val-1 & location<val5+1 & distribution == "narrow")
  totalLength <- length(ofInterest$location)
  lengthNarrow <- length(ofInterestNarrow$location)
  likelihood[((val+5)/5),1] <- signif(lengthNarrow/totalLength, digits = 5) #narrow
  likelihood[((val+5)/5),2] <- signif(1 - lengthNarrow/totalLength, digits = 5) #wide
}

#replace any NaN value with 0 for computation
likelihood[is.nan(likelihood)] <- 0

#make choice based on likelihood at that point 
locationBin <- c(1:160)

for (i in locationBin){
  if(likelihood[i,1]>likelihood[i,2]){
    likelihood[i,3] <- "narrow"
  }else if (likelihood[i,1]==likelihood[i,2]){
    if (i < 80){
      likelihood[i,3] <- "narrow"
    } else{
      likelihood[i,3] <- "wide"
    }
  }else{
    likelihood[i,3] <- "wide"
  }
  #make confidence report equal to the likelihood at that location
  if (likelihood[i,3] == "narrow"){
    likelihood[i,4] <- likelihood[i,1]#set confidence
  }else{
    likelihood[i,4] <- likelihood[i,2]#set confidence
  }
}

#look at difference in confidence between means of each distribution
likelihood[,5] <-seq.int(0, 795, 5) #set location bin value

dataframe <- as.data.frame(likelihood)
dataframe <- dataframe %>% 
  rename(
    narrowLikeli = V1,
    wideLikeli = V2,
    choice = V3,
    confidence = V4,
    location = V5
  )
dataframe$confidence <- signif(as.numeric(dataframe$confidence), digits = 6)
dataframe$location <- as.numeric(dataframe$location)

dataFrame <- subset(dataframe, location<500 & location>195)
ggplot(data = dataFrame, mapping = aes(location, confidence),)+ 
  geom_point()

#find the difference in confidence reports between distributions at each step from the decision bound
midData <- subset(dataframe, location < 550 & location > 150)
midlineIndex <- which.min(midData$confidence)
midline<- midData$location[midlineIndex]

distance <- c(1:30)
confDiffsOptimal <- matrix(data = NA, nrow = length(distance), ncol = 1)
for(val in distance){
  narrowSide <-  subset(dataFrame, location == (midData$location[midline] - (val*5)))
  narrowConf <- narrowSide$confidence /(as.numeric(narrowSide$narrowLikeli) + as.numeric(narrowSide$wideLikeli))
  
  wideSide <- subset(dataFrame, location == (midData$location[midline] + ((val)*5)))
  wideConf <- wideSide$confidence /(as.numeric(wideSide$narrowLikeli) + as.numeric(wideSide$wideLikeli))
  
  confDiffsCurrent[val,1] <- abs(narrowConf- wideConf)
}

#visualise these differences
  x <- confDiffsOptimal[,1]
  y <- seq.int(0,145,5)
  plot(y,x, xlab = "Distance from bound", ylab = "Difference in confidence", main = "Optimal Pair") #going from midline to means (the differences drop off heading to the means)
  # Plot histogram of x
  hist(x,
       main = "Optimal Pair",
       xlab = "Confidence Difference",
       xlim = c(0, 1))

    
################################################################
#thinking about sampling: bins of 10 either side of this midline (320)

#starting inside the means - where we will take 90% of trials from
data <- subset(dataframe, location < 505 & location > 195)
head(data)  
#pick out the 10pixel bin slots from this
data <- subset(data, (location/10)%%1 == 0)

#have the likelihoods at each bin point of each distribution type
#use these to generate trials according to those probabilities at each point. 
k_trials <- c(1:10) #10 trial per bin location
locations <- c(1:length(data$location)) 
sim_trials <- matrix(data = 0, nrow = length(k_trials), ncol = length(data$location))


for (i in locations){
  for (j in k_trials){
    flip <- runif(1,0,1) #generate a random likelihood
    if(flip < data$confidence[i]){
      sim_trials[j,i] <- data$choice[i]
    } else {
      if (data$choice[i] == "narrow"){
        sim_trials[j,i] <- "wide"
      } else {
        sim_trials[j,i] <- "narrow"
      }
    }
  }
}

narrow_count <- sum(str_count(sim_trials, pattern = "narrow")) #138
wide_count <- sum(str_count(sim_trials, pattern = "wide")) #172
total_count <- length(str_count(sim_trials, pattern = "narrow")) #310

narrow_prop <- narrow_count/total_count  #~45% 
wide_prop <- wide_count/total_count  #~55%
#between the means it is ~5% more likely to be from the wide distribution than the narrow. 
#if we take this to be 90% of the overall trial then - we have:
#310 trials already simulated (90%) - 34.4 trial left (34 for even number)

narrow_inner <- .9*narrow_prop
wide_inner <- .9*wide_prop
#~40.5% of total as narrow and ~49.5% of total as wide. 
 
narrow_remain <- .5-narrow_inner 
wide_remain <- .5-wide_inner
#34 trials left of which 9.5% of 334 should be narrow. (32 and 2, 172, 174) aka 34 narrow and 0 wide 
#leaves us with 172 wide and 172 narrow


#what if 400 trials total? 310 already sampled, 90 left.
#.775 of trials sampled
narrow_inner <- .775*narrow_prop
wide_inner <- .775*wide_prop
narrow_remain <- .5-narrow_inner #15.5% of 400 trials to reach balance  = 62 narrow + 128 = 200
wide_remain <- .5-wide_inner #7% of 400 trials to reach balance = 28 wide narrow + 172 = 200




  
  