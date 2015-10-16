source("ptw_func.R")
source("read.mtx.R")

# set dimensions for gpro
d <- 11
#Set the function Name here 
funcname <- 'ptw'

#Read in the base design.
xf <- read.mtx('XD11N220.mtx')   
colnames(xf) <- paste('x',1:d,sep='')

#Compute test set data
xp <- data.matrix(read.table('test11.csv',header=F,sep=',')[,1:d])
colnames(xp) <- colnames(xf)
yp <- data.frame(y = do.call(funcname,list(xp)))

ptw.test <- data.frame(xp, y = yp)
write.table(ptw.test, file = "ptw_test.csv", sep = ",", row.names = FALSE)

# Permutation on the base design
ocn <- colnames(xf)
cp.m <- data.frame(read.csv("order_ptw.csv",  h=F)[, 1:d])

# purmute training data 24 times
Ptrain <- as.data.frame(matrix(NA, ncol = d + 1, nrow = dim(xf)[1]*25))
for(i in 1:25) {
	cp <- as.numeric(cp.m[i,])
	xf_train = data.matrix(xf[,cp])
	yf_train = data.frame(y = do.call(funcname, list(xf_train)))
	train = data.frame(xf_train, y = yf_train)
	Ptrain[(((i-1)*220)+1):(i*220), 1:12] <- train
}
colnames(Ptrain) = c(ocn, "y")  


write.table(Ptrain, file = "ptw.csv", sep = ",", row.names = F)
