source("Gpro_func.R")

# set dimensions for gpro
d <- 4
#Set the function Name here 
funcname <- 'GProtein'

#Read in the base design.
xf <- read.table('XD4N40.txt',header=F)   
colnames(xf) <- paste('x',1:d,sep='')

#Compute test set data
xp <- data.matrix(read.table('test4.csv',header=F,sep=',')[,1:4])
colnames(xp) <- colnames(xf)
yp <- data.frame(y = do.call(funcname,list(xp)))

# Permutation on the base design
ocn <- colnames(xf)
cp.m <- data.frame(read.csv("order_gpro.csv",  h=F)[, 1:4])

# purmute training data 24 times
gpro_train <- as.data.frame(matrix(NA, ncol = d + 1, nrow = dim(xf)[1]*24))
for(i in 1:24) {
	cp <- as.numeric(cp.m[i,])
	xf_train = data.matrix(xf[,cp])
	yf_train = data.frame(y = do.call(funcname, list(xf_train)))
	train = data.frame(xf_train, y = yf_train)
	gpro_train[(((i-1)*40)+1):(i*40), 1:5] <- train
}
colnames(gpro_train) = c(ocn, "y")  


write.table(gpro_train, file = "gpro.csv", sep = ",", row.names = F)
	

