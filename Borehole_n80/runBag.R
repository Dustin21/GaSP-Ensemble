runBag <- function(x.train, y.train, y.test, iterations = c(5, 10), size = c(30, 35))
{
# Run GP bagging procedure over a vector of iterations according to a vector of sizes.
# x.train and y.train should be lists of training data.frames, as what would be
# produced from a Latin hypercube sampling scheme. Output is a data.frame of the
# RMSE and MaxE for each iteration, size, and train set in the x.train list.
 
    # dependencies
    require(dplyr)
    require(foreach)
    source("bagGP.R")
    source("errorGP.R")
    
    loop.bag <- foreach(i = 1:length(y.train)) %do%
    {
	# loop over bagging procedure according to iterations, size and train set 
	errors <- errorGP(x.train[[i]], y.train[[i]], y.test, iterations, size)

	# create data.frame with a row i indicating the train
	# set for the current run
	n.train.errors <- data.frame(errors, train_set = i)
	
	return(list(n.train.errors))
    }

    # recursively rbind loop.bag list
    bag.out <- ldply(loop.bag[[1]])
    for(i in 2:length(y.train))
    {
	bag.out <- rbind(bag.out, ldply(loop.bag[[i]]))
    }

    # transform to data.frame and factor variables
    bag.out <- data.frame(bag.out) %>% 
		mutate(RMSE, MaxE, 
			iterations = as.factor(iterations), 
			size = as.factor(size), 
			train_set= as.factor(train_set))

    return(bag.out)
}
