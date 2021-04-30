rm(list=ls())

calc_avg <- function(dat, ffile, first_2_cols){

  y0 <- median(dat$HON0_avg)
  se0 <- median(dat$HON0_se)
  m0 <- median(dat$HON0_median)

  y2 <- median(dat$HON2_avg)
  se2 <- median(dat$HON2_se)
  m2 <- median(dat$HON2_median)

  y3 <- median(dat$HON3_avg)
  se3 <- median(dat$HON3_se)
  m3 <- median(dat$HON3_median)

  sstr <- sprintf("%s (median)\t1\t1\tNA\t%.3f\t%.3f\t%.3f\tNA\tNA\t%.3f\t%.3f\t%.3f\tNA\tNA\t%.3f\t%.3f\t%.3f\tNA\tNA", 
    first_2_cols, y3,se3,m3, y2,se2,m2, y0,se0,m0)

  write(sstr, file=ffile, append=TRUE)

 }


args <- commandArgs(TRUE)
input <- args[1]
output <- args[2]
dat0 <- read.table(args[1], header = T, sep = '\t')

clade <- dat0$clade[1]
calc_avg(dat0[dat0$classification=="charge",], output, sprintf("charge\t%s", clade))
calc_avg(dat0[dat0$classification=="MY",], output, sprintf("MY\t%s", clade))

