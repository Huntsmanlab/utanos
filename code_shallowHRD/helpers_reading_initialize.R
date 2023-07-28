#' Finds the indices of big segments that are close ratio_median-wise
#'
#' @description
#' Starting from the largest segment, we find big segments if their ratio_median
#' difference is less than the the given threshold. We save their indices.
#'
#' @param segments_copy A data frame containing segment data. For the first pass, these segments
#' are those with size >= 3mb.
#' @param thr A float: the threshold previously estimated via KDE. Used to determine whether
#' big segments are merged.
#'
#' @export

FindBigSegmentsIndices <- function(segments_copy, thr) {
  # Getting largest segment by size
  largest_segment_index = segments_copy[which.max(segments_copy[,7]), 1]
  largest_segment_ratio = segments_copy[which.max(segments_copy[,7]), 6] + 0.00001
  
  #Getting segment index whose ratio_median is closest to largest_segment
  closest_segment_index = segments_copy[which.min(abs(segments_copy[,6] - largest_segment_ratio)),1]
  
  segments_indices = c(largest_segment_index, closest_segment_index)
  
  #Finding segments with ratio_medians close to segments in segments_indices
  for (i in 1:dim(segments_copy)[1]) {
    num_of_times_close_to_segments = 0
    # Iterate through segments in segments_indices. We compare segment i against all of these
    # We add +1 to num every time segment i and the other segment have a ratio_median difference
    # less than the threshold
    for (j in 1:length(segments_indices)) {
      index = which(segments_copy[,1] == segments_indices[j])
      if (abs(copy_segments[i,6] - copy_segments[index,6]) < thr) {
        num_of_times_close_to_segments = num_of_times_close_to_segments + 1
      }
    }
    
    # If segment i is close to every single segment in the list, then we add segment i
    # to the list too. 
    if (num_of_times_close_to_segments == length(segments_indices)) {
      segments_indices = c(segments_indices, segments_copy[i,1])
    }
  }
  
  segments_indices = unique(segments_indices) # make sure there's no repeats
  segments_indices
}













## Function readSegmFile ##
#' TODO: write docs
ReadSegmFile<-function(seg_file_name){
  aa<-scan(seg_file_name,nlines=1,what="character",quiet = T)
  k<-1
  
  while((length(grep("#",aa))!=0)&&(grep("#",aa)==1)){
    aa<-scan(seg_file_name,skip=k,nlines=1,what="character",quiet = T)
    k<-k+1
  }
  
  tmp<-read.delim(seg_file_name,skip=k-1, header=T,sep="\t")
  tmp<-data.matrix(tmp)
  if(k==2){
    results<-getInfoSegm(scan(seg_file_name,nlines=1,what="character",quiet = T))
    homoConst<-results$homoConst
    sdH<-results$sdH	
    p_BAF<-results$p_BAF
    q_LRR<-results$q_LRR
    Delta<-results$Delta
  }else{homoConst<-1;sdH<-NA;p_BAF<-NA;q_LRR<-NA;Delta<-c(2,NA)}
  
  out<-list(tmp=tmp,homoConst=homoConst,sdH=sdH,p_BAF=p_BAF,q_LRR=q_LRR,Delta=Delta)
}

## getSegmentID : gather in the same name segment of the length  ##
#' TODO: write docs
GetSegmentID<-function(THR,tmp,c_chr,c_cn,c_conf){
  
  tmp[,c_conf]<-0;tmp[1,c_conf]<-1
  
  for(k in 2:dim(tmp)[[1]]){
    
    if(tmp[k,c_chr]==tmp[k-1,c_chr] & abs(tmp[k,c_cn]-tmp[k-1,c_cn])<THR){
      
      tmp[k,c_conf]<-tmp[k-1,c_conf]
      
    }else{tmp[k,c_conf]<-tmp[k-1,c_conf]+1}
  } 
  
  out<-tmp
  
}

#' TODO: write docs
GetInfoSegm <- function(infoString){
  
  tt<-unlist(strsplit(infoString,";"))
  hh<-grep("homoConst=",tt)
  if(length(hh)==1){homoConst<-as.numeric(sub("homoConst=","",tt[hh]))}else{homoConst<-1}
  hh<-grep("sdH=",tt)
  if(length(hh)==1){sdH<-as.numeric(sub("sdH=","",tt[hh]))}else{sdH<-NA}
  
  hh<-grep("p=",tt)
  if(length(hh)==1){p_BAF<-as.numeric(sub("p=","",tt[hh]))}else{p_BAF<-NA}
  
  hh<-grep("q=",tt)
  if(length(hh)==1){q_LRR<-as.numeric(sub("q=","",tt[hh]))}else{q_LRR<-NA}
  
  hh<-grep("2LR=",tt)
  if(length(hh)==1){Delta<-as.numeric(sub("2LR=","",tt[hh]))}else{Delta<-NA}
  Delta<-c(2,Delta)
  
  
  if((is.na(homoConst))||(!is.numeric(homoConst))){homoConst<-1
  }else{if((homoConst<0.8)||(homoConst>1)){homoConst<-1}
  }
  if((is.na(sdH))||(!is.numeric(sdH))){sdH<-NA
  }else{if((sdH<0.05)||(sdH>1)){sdH<-NA}
  }
  
  out<-list(homoConst=homoConst,sdH=sdH,p_BAF=p_BAF,q_LRR=q_LRR,Delta=Delta)
  
}