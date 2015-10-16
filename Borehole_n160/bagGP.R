bagGP <- function(x, y, iterations, size, replace = FALSE, seed = 1)
# bagGP function performs bootstrap aggregation (Bagging) on the data inputs
# x and y. The sample size and number of iterations are inputed via the 
# function parameters. For deterministic Bagging, replacement is set to FALSE.
{
    # dependencies
    require(plyr)
    require(dplyr)
    require(foreach)
    
    # coerce into data.frame
    df <- data.frame(x, y = y)

    predictions <- foreach(m = 1:iterations, .combine = rbind) %dopar% 
	    {
		# set temporary wd
		setwd("gasp_output")		

		# repeatibility
		set.seed(m+seed)
		
		# randomly sample from data.frame
		sample.df <- sample_n(df, size, replace = replace)

		# write sample.df to .mtx for cpp call
		yc <- data.frame(y = sample.df$y)
		xc <- sample.df[,-ncol(sample.df)]
		rownames(xc) <- 1:nrow(xc)
		write.mtx(yc, file = paste0("yc", m, ".mtx"))
		write.mtx(xc, file = paste0("xc", m, ".mtx"))

		# create fit.gsp file for C input
		gaspfit(m = m)

		# fit GaSP to sample.df
		fits <- paste0("./gasp ", "fit", m, ".gsp")
		system(fits)

		# read prediction outputs from gasp
		ypred <- read.mtx(paste0("ypr", m, ".mtx"))
		pred <- ypred[,1]
		se <- ypred[,2]

		system.remove <- paste0("rm ", "ypr", m, "* ", "fit", m, "* ", 
					"xc", m, "* ", "yc", m, "*")
		system(system.remove)
		
		# revert back to original wd
		setwd("..")

		return(list(pred, se))
	    }

    return(predictions)
}
