#-----------------------------------------------------
#---------------- Dependencies -----------------------
#-----------------------------------------------------
suppressPackageStartupMessages({
    library(plyr)
    library(dplyr)
    library(foreach)
    library(reshape2)
    library(ggplot2)
    library(doMC) # registerDoMC(cores = 2)
    source("read.mtx.R")
    source("write.mtx.R")
    source("normalised_errors.R")
    source("errorGP.R")
    source("bagGP.R")
    source("gasp_output/gaspjob.R")
    source("runBag.R")
})

#-----------------------------------------------------
#---------------- Data Inputs ------------------------
#-----------------------------------------------------
gpro <- read.csv("borehole.csv", header = TRUE)

# split data into list of data.frames
nVars <- ncol(gpro)
x.train <- split(gpro[,-nVars], (seq(nrow(gpro)) - 1) %/% 80)
y.train <- split(gpro[,nVars], (seq(nrow(gpro)) - 1) %/% 80)

# test set
gpro.test <- read.csv("borehole_test.csv", header = TRUE)
x.test <- select(gpro.test, -y)
y.test <- select(gpro.test, y)

# write test sets to file for gasp
write.mtx(x.test, file = "gasp_output/xp.mtx")
write.mtx(y.test, file = "gasp_output/yp.mtx")

#-----------------------------------------------------
#---------------- Single GP Fit ----------------------
#-----------------------------------------------------
N <- length(x.train)

singleGP <- foreach(m = 1:N, .combine = rbind) %dopar% 
{
    # set temporary wd
    setwd("gasp_output")

    # read training set and write to file for gasp
    yc <- data.frame(y = y.train[[m]])
    xc <- x.train[[m]]
    rownames(xc) <- 1:nrow(xc)
    write.mtx(xc, file = paste(paste("xc", m, sep=""), "mtx", sep="."))
    write.mtx(yc, file = paste(paste("yc", m, sep=""), "mtx", sep="."))
    
    # output a fit.out file indexed by m
    gaspfit(m=m)
  
    # make a system call to run gasp with fit.out
    fits <- paste0("./gasp", " ", "fit", m, ".gsp") 
    system(fits)

    # import predictions error
    read.predict <- paste0("ypr", m, ".mtx")
    singleGP_pred <- read.mtx(read.predict)[,1]
    
    # remove the output of each loop after importing the predictions
    system.remove <- paste0("rm ", "ypr", m, "* ", "fit", m, "* ", "xc", m, "* ", "yc", m, "*")
    system(system.remove)

    # return to original directory	
    setwd("..")
    
    # compute rmse and maxe
    RMSE.single <- rmse(singleGP_pred, y.test$y, yc$y, normalized = TRUE)
    MaxE.single <- maxe(singleGP_pred, y.test$y, yc$y, normalized = TRUE)
    
    return(list(pred = singleGP_pred, RMSE = RMSE.single, MaxE = MaxE.single))	
}

# normalised prediction error
errors.single <- data.frame(RMSE = ldply(singleGP[,2])[,2], MaxE = ldply(singleGP[,3])[,2]) %>%
			mutate(fits = paste("fit", 1:N, sep = ""))


#-----------------------------------------------------
#---------------- Bagged GP Fit ----------------------
#-----------------------------------------------------

#****** Inputs ********
n.iterate <- c(2, 3) 
n.size <- c(70,75)     
#**********************

# run bagging algorithm - number of iterations for each size
bagged.error <- runBag(x.train, y.train, y.test, iterations = n.iterate, size = n.size)


#-----------------------------------------------------
#--------------------- Plots -------------------------
#-----------------------------------------------------

# plot RMSE and MaxE for single GP fit
errors.single.plot <- errors.single %>% 
			melt(id.vars = "fits", variable.name = "error.type", 
				value.name = "error") %>%
			group_by(error.type) %>%
			ggplot(aes(y = error, x = error.type, colour = error.type)) + 
			geom_point() + facet_wrap(~error.type, scales = "free") +
			theme_bw()

# plot RMSE by size and iterations
rmse.plot <- bagged.error %>%
		select(-MaxE) %>%
		group_by(size, iterations) %>%
		ggplot(aes(x = iterations, y = RMSE)) +
		geom_boxplot(fill = "red") + facet_wrap(~size) +
		geom_boxplot(data = errors.single, fill = "red", aes(RMSE, x = "single")) + 
		scale_x_discrete(limits = c("single", n.iterate)) +	
		theme_bw()

# plot MaxE by size and iterations
maxe.plot <- bagged.error %>%
		select(-RMSE) %>%
		group_by(size, iterations) %>%
		ggplot(aes(x = iterations, y = MaxE)) +
		geom_boxplot(fill = "blue") + facet_wrap(~size) +
		geom_boxplot(data = errors.single, fill = "blue", aes(MaxE, x = "single")) +
		scale_x_discrete(limits = c("single", n.iterate)) +	
		theme_bw()


#-----------------------------------------------------
#---------------- print outputs ----------------------
#-----------------------------------------------------

# print plots
print(rmse.plot)
print(maxe.plot)

# write files to wd
write.csv(errors.single, file = "error_single.csv")
write.csv(bagged.error, file = "error_bagged.csv")


