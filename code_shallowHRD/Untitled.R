line_largest_index = which.max(cop_tmp[,7])             # biggest segment
largest_index = cop_tmp[line_largest_index,1]           # associated index
ratio_largest = cop_tmp[line_largest_index,6]           # median_ratio segment
ratio_largest = ratio_largest + 0.00001                 # avoid picking same value

close_index = cop_tmp[which.min(abs(cop_tmp[,6] - ratio_largest)),1] #they find the index of the segments whose ratio is closes to ratio_largest
closest_index = c(largest_index, close_index)

## Find all value that matches the two closest index fit with THR 

L=dim(cop_tmp)[1]

# so essentially we are going to iterate through all segments, and if a segment
# has validation 2, then that means that the difference in ratio medians between the segment
# and the largest/closest is less than the threshold.
# if this is the case we're adding the index of the segment to closest_index
for (i in 1:L){
  validation = 0
  n = length(closest_index)
  for (g in 1:n){
    if (abs(cop_tmp[i,6] - cop_tmp[which((cop_tmp[,1]==closest_index[g])),6]) < THR){   
      validation = validation + 1}
  }
  if (validation == n){
    closest_index = c(closest_index, cop_tmp[i,1])}}

closest_index = unique(closest_index)