write.mtx <- function(x, filename, col.name = "x", page.width = 80)
# Write matrix / data frame x to filename in blocked ACED/GaSP format.
# Row and column labels are preserved if present.
# If there are no row labels, then 1, 2, ... will be used.
# If there are no column labels, then x1, x2, ... will be
# used, unless another name is given for col.name.
# 1997.01.24: replaced with new version.
# 2001.11.08: Renamed write.mx (from write.mat).
# 2002.10.17: Renamed write.mtx (from write.mx).
{
     # Make sure x is a matrix.
     x <- as.matrix(x)

     nrows    <- nrow(x)
     ncolumns <- ncol(x)

     row.labels <- dimnames(x)[[1]]
     if (length(row.labels) == 0)
          row.labels <- c("----", "Case", "----", 1:nrows)
     else
          row.labels <- c("----", "Case", "----", row.labels)

     col.labels <- dimnames(x)[[2]]
     if (length(col.labels) == 0)
          col.labels <- paste(col.name, 1:ncolumns, sep="")

     spaces <- c(rep("   ", ncolumns))

     mat <- rbind(spaces, col.labels, spaces, x)
     dimnames(mat) <- list(row.labels, c(rep(" ", ncol(mat))))

     sink(filename)
     options(length=nrow(mat)+10)
     # Prevent repeats of column labels.
     temp <- options(length = nrow(mat) + 1, width = page.width)

     print(mat,quote=F)
     sink()
     cat(filename, ": ", nrow(mat) - 3, " rows and ", ncol(mat),
               " columns written.\n", sep = "")

     # Restore length.
     options(temp)
}


