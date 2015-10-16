source("borehole_func.R")

# set dimensions for gpro
d <- 8
#Set the function Name here 
funcname <- 'bore.function'

#Read in the base design.
xf <- read.table('XD8N40.txt',header=F)   
colnames(xf) <- paste('x',1:d,sep='')

#Compute test set data
xp <- data.matrix(read.table('test8.csv',header=F,sep=',')[,1:d])
colnames(xp) <- colnames(xf)
yp <- data.frame(y = do.call(funcname,list(xp)))

borehole.test <- data.frame(xp, y = yp)
write.table(borehole.test, file = "borehole_test.csv", sep = ",", row.names = FALSE)

# Permutation on the base design
ocn <- colnames(xf)
cp.m <- data.frame(read.csv("order_borehole.csv",  h=F)[, 1:d])

# purmute training data 24 times
train <- as.data.frame(matrix(NA, ncol = d + 1, nrow = dim(xf)[1]*25))
for(i in 1:25) {
	cp <- as.numeric(cp.m[i,])
	xf_train = data.matrix(xf[,cp])
	yf_train = data.frame(y = do.call(funcname, list(xf_train)))
	train = data.frame(xf_train, y = yf_train)
	gpro_train[(((i-1)*40)+1):(i*40), 1:9] <- train
}
colnames(train) = c(ocn, "y")  


write.table(train, file = "borehole.csv", sep = ",", row.names = F)
