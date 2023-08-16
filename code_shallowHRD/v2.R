##### NEW Find Threshold 2 #####

cat("Find new Threshold 2... \n")


graphe_V_tab = data.frame(graphe_V_tab)

graphe_V_tab$size = graphe_V_tab$end - graphe_V_tab$start + 1

copy_graphe_V_tab = graphe_V_tab
copy_graphe_V_tab = copy_graphe_V_tab[,-8]
copy_graphe_V_tab = copy_graphe_V_tab[,-7]
copy_graphe_V_tab = copy_graphe_V_tab[,-1]

copy_graphe_V_tab$size <- copy_graphe_V_tab$end - copy_graphe_V_tab$start + 1

X = copy_graphe_V_tab

Y = X

Y = Y[which(Y$chr != 23),]


L=dim(Y)[1]
X=data.matrix(Y)

test=c()
i=1
for (i in 1:L){
  v=i+1
  for (y in v:L-1){
    test = c(test, abs(Y[i,5] - Y[y,5]))}
}  



### multiple test THR

list_THR = c()
list_max2 = c()
c = 0


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
    if (abs((First_local_minima-First_local_maxima)-(Second_local_maxima-First_local_minima)) < 0.05){
      list_THR = c(list_THR,First_local_minima)
      list_max2 = c(list_max2, Second_local_maxima)
      test_optimized = sampled_test
    }
  }
  c = c + 1
}


### Faire second inflexion point if no THR found


if (is.null(list_THR) == TRUE){
  
  ## First maxima
  
  maxx = localMaxima(density(test)$y)[1]
  First_maxima = density(test)$x[maxx]
  
  
  cat("First_minima... \n")
  
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
  
  Threshold_second_round_simulation = Threshold
  
  MAX2 = mean(list_max2)
}



###  Classical THR + graphe

X = copy_graphe_V_tab
# size > 3Mb
X = X[which(X[,6] > 2999999),]

# Segments are defined as ‘large’ if size > (Q3-Q1)/2
X = X[which(X[,6] > ((quantile(X[,6])[[4]] - quantile(X[,6])[[2]])/2)),]  


X = X[which(X$chr != 23),]

L=dim(X)[1]
X=data.matrix(X)

test=c()
i=1
for (i in 1:L){
  v=i+1
  for (y in v:L-1){
    test = c(test, abs(X[i,5] - X[y,5]))}
}  

## First maxima

maxx = localMaxima(density(test)$y)[1]
First_local_maxima = density(test)$x[maxx]


## First minima 

minx = localMinima(density(test)$y)
First_local_minima = density(test)$x[minx][which(density(test)$x[minx] > First_local_maxima)][1]


## Second maxima

maxx = localMaxima(density(test)$y)[2]
Second_local_maxima = density(test)$x[maxx]


if (is.na(Second_local_maxima) == FALSE){
  if (abs((First_local_minima-First_local_maxima)-(Second_local_maxima-First_local_minima)) < 0.10){
    if (First_local_minima < Threshold){
      Threshold = First_local_minima
      Threshold_global_view = First_local_minima
      MAX2 = Second_local_maxima
    }
  }
}