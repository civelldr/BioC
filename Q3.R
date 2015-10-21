## ExpressionSet => eSet (ex: 1 matrix is methylated channel, another matrix is the unmethylated channel)
#
#
# (1)
sample_5 <- ALL[,5] # get just sample #5
mean(exprs(sample_5))


# (6)
mean(airway$avgLength)

# (7)
sample_3 <- airway[,3]
counts <- assay(sample_3, "counts")
sum(counts >= 1)
