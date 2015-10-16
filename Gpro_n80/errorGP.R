cumulativeError <- function(gatherGP, iterations, y.test, y.train)
# Input a bagGP object output as demonstrated in the function iterationError,
# then aggregate subgroups of the total bootstrap sample runs with the
# sub-group sizes indicated by the corresponding element of the iterations vector.
{
    cumError <- foreach(k = 1:length(iterations), .combine = rbind) %do%
    {
	# gather results from bagGP output and weight predictions according to
	# the inverse of the variance.
	cumResults <- gatherGP %>%
		      filter(.id %in%  paste("result", 1:iterations[k], sep = ".")) %>%
		      group_by(run) %>%
		      summarise(pred.weight = sum(pred/se^2), se2 = sum(1/se^2)) %>%
		      mutate(pred.test = pred.weight/se2) %>%
		      select(run, pred.test)
    
	# compute the RMSE and MaxE via the rmse and maxe functions
	errors <- cumResults %>%
		  summarise(RMSE = rmse(ypred = pred.test, ytest = y.test, ytrain = y.train),
			    MaxE = maxe(ypred = pred.test, ytest = y.test, ytrain = y.train))

	return(list(errors))
    }
   
    return(cumError)
}


iterationError <- function(x, y, y.test, iterations = 25, size, replace = FALSE, seed = 1)
# Run bagGP over the number of times indicated by the max of a vector of iterations. The 
# iterations are then grouped cumulatively according to each member of the iterations
# vector.  
{
    # dependencies
    require(plyr)
    require(dplyr)
    source("bagGP.R")
    source("normalised_errors.R")
    
    # perform bootstrap aggregation according to max of iterations
    bagGP_output <- bagGP(x, y, max(iterations), size, replace, seed)
	
    # cumulatively aggregate bootstrap samples according according to number of
    # iterations. If iterations = c(2,5), then aggregate two runs, then aggregate
    # 5 runs. Finally structure bagGP output into clean data.frame
    gather_output <- melt(ldply(bagGP_output[,1]), id.vars = ".id") %>%
			    data.frame(melt(ldply(bagGP_output[,2]), id.vars = ".id")) %>%
			    select(-variable.1, -.id.1) %>%
			    setNames(c(".id", "run", "pred", "se"))

    # apply RMSE and MaxE to each bootstrapped group indicated above
    cumError <- cumulativeError(gather_output, iterations, y.test, y.train = y)
    
    return(cumError)
}

errorGP <- function(x, y, y.test, iterations = c(5, 10), size = c(5, 10), 
		replace = FALSE, seed = 1)
# loop the functions above over size and output a data.frame with factor variables
# according to training set, iterations and size.
{
    error_by_size <- foreach(m = 1:length(size), .combine = rbind) %do%
    {
	# compute iteration error for each element of the size vector
	error_by_iterations <- iterationError(x, y, y.test, iterations, size[m], 
					replace, seed) %>% ldply
	
	return(list(error_by_iterations))
    }
    
    # conform data.frame to desired output
    error.df <- ldply(error_by_size) %>%
		data.frame(iterations = rep(iterations, length(size)),
			   size = rep(size, each = length(iterations)))
		 
    return(error.df)
}
