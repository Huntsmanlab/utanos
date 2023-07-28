if (c+2 > L_3mb){         # anticipate end of file
  gr=GRanges(seqnames=c(tmp_3mb[c+1,2]),
             ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c+1,5])),
             strand=c("*")) # to overlap
  subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
  
  tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_3mb[c+1,1], tmp_3mb[c+1,2], tmp_3mb[c+1,3], tmp_0.1_3mb[i,4], tmp_3mb[c+1,5], median(subsetGRobject$ratio),
                                  tmp_3mb[c+1,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c+1,8])) 
  i=i+1
  c=L_3mb+1
  # add data after c
  # essentially merging small and next_large
}
else{
  # if not end of file
  gr=GRanges(seqnames=c(tmp_3mb[c+1,2]),
             ranges=IRanges(start=c(tmp_0.1_3mb[i,4]),end=c(tmp_3mb[c+1,5])),
             strand=c("*")) # to overlap
  subsetGRobject = subsetByOverlaps(GRanges_object_ratio_file, gr)
  
  tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_3mb[c+1,1], tmp_3mb[c+1,2], tmp_3mb[c+1,3], tmp_0.1_3mb[i,4], tmp_3mb[c+1,5], median(subsetGRobject$ratio),
                                  tmp_3mb[c+1,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c+1,8]) , tmp_3mb[(c+2):L_3mb,])  
  i=i+1
  c=L_3mb+1
  # same as before, but between c and c+2
  
}