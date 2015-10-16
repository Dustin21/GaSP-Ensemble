rmse <- function(ypred, ytest, ytrain, normalized = TRUE)
	# 2013.02.20
{
	rmse <- sqrt(sum((ytest - ypred)^2, na.rm = TRUE) / sum(!is.na(ytest)))
	if (normalized)
	{
		rmse.mean <- sqrt(sum((ytest - mean(ytrain))^2, na.rm = TRUE) /
				sum(!is.na(ytest)))
		rmse <- rmse / rmse.mean
	}
	return(rmse)
}

maxe <- function(ypred, ytest, ytrain, normalized = TRUE)
	# 2014.03.28
{
	maxe <- max(abs(ytest - ypred), na.rm = TRUE)
	if (normalized == TRUE)
	{
		maxe.mean <- max(abs(ytest - mean(ytrain)), na.rm = TRUE)
		maxe <- maxe / maxe.mean
	}
	return(maxe)
}
