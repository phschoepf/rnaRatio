#read input
input <- read.delim(file.choose(new=FALSE))

filtered <- input[which(input[, 3] != "NA"), ]
pearval <- cor(filtered$MYC, filtered$BASP1)
write.table(c(filtered[1,1], pearval), output.csv, append=TRUE)