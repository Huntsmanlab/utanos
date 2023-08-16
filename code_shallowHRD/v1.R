##### Find Threshold 1 #####

cat("Find Threshold 1... \n")

X = read.table(paste0(outputPath,"/",NAMEEE,"_ratio_median_gathered.txt"), sep = "\t", header = TRUE)

X = X[which(X$chr != 23),]


# true value

B <- data.frame(ratio_file_tsv)
B=B[,-1]
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")


GRanges_object_ratio_file = makeGRangesFromDataFrame(B,keep.extra.columns = TRUE, ignore.strand = TRUE, seqinfo = NULL,
                                                     seqnames.field = "chr",
                                                     start.field = "start", end.field = "end")
i=1
L=dim(X)[1]

while (i < L){
  gr=GRanges(seqnames=c(X[i,1]),
             ranges=IRanges(start=c(X[i,3]),end=c(X[i,4])),
             strand=c("*")) # to overlap
  subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
  X[i,5] = median(subsetGRobject$ratio)
  i=i+1
}


X = X[which(X[,6] > 2999999),]

L=dim(X)[1]
X=data.matrix(X)

test=c()
i=1
for (i in 1:L){
  v=i+1
  for (y in v:L-1){
    test = c(test, abs(X[i,5] - X[y,5]))}
}  


list_THR = c()
c = 0

nb_simulations = 100000


while (c < nb_simulations){
  n = length(test)
  n
  
  rand_n = sample(2:n,1)
  rand_n
  
  sampled_test = sample(test,rand_n)
  
  ## First maxima
  
  maxx = localMaxima(density(sampled_test)$y)[1]
  First_local_maxima = density(sampled_test)$x[maxx]
  
  
  ## First minima 
  
  minx = localMinima(density(sampled_test)$y)
  First_local_minima = density(sampled_test)$x[minx][which(density(sampled_test)$x[minx] > First_local_maxima)][1]
  
  
  ## Second maxima
  
  maxx = localMaxima(density(sampled_test)$y)[2]
  Second_local_maxima = density(sampled_test)$x[maxx]
  
  if (is.na(Second_local_maxima) == FALSE){
    if (abs((First_local_minima-First_local_maxima)-(Second_local_maxima-First_local_minima)) < 0.10){
      list_THR = c(list_THR,First_local_minima)
      test_optimized = sampled_test
    }
  }
  c = c + 1
}



if (is.null(list_THR) == TRUE){
  
  ## First maxima
  
  maxx = localMaxima(density(test)$y)[1]
  First_maxima = density(test)$x[maxx]
  
  ## First minima
  
  minx = localMinima(density(test)$y)
  First_local_minima = density(test)$x[minx][which(density(test)$x[minx] > First_maxima)][1]
  
  h <- hns(test, deriv.order = 2)
  den <- kdde(test, h=h, deriv.order = 2)
  
  flections<-c()
  for(i in 2:length(den$estimate)){
    if(sign(den$estimate[i])!=sign(den$estimate[i-1])){
      flections<-c(flections,i)
    }
  }
  
  All_inflexion_point = den$x[flections][which(den$x[flections] > First_maxima)]
  Threshold_change_sign_second_derivative_2nd_value = sort(All_inflexion_point)[1]
  
  Threshold = min(Threshold_change_sign_second_derivative_2nd_value, First_local_minima)
  
}


if (is.na(mean(list_THR)) == FALSE){
  number_positive = length(list_THR)
  
  Threshold = mean(list_THR)
  
  Threshold_first_round_simulation = Threshold
}