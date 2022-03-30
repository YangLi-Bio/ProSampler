##########################################################
#                                                        #
#             Convert p-values to z-scores               #
#                                                        #
##########################################################

library(pbmcapply)

multi.n <- 1e-09 # the multiplication factor

res.df <- Reduce("rbind", pbmclapply (1:1000, function(i) {
  n <- i * multi.n
  zscore <- qnorm(n, lower.tail = F)
  c(n = i, pvalue = n, zscore = zscore)
}, mc.cores = detectCores()))

colnames(res.df) <- c("approx", "pvalue", "zscore")

write.csv(res.df, "pvalue_and_zscore.csv", row.names = F)