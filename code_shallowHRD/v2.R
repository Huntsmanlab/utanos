##### NEW Fast gathering 2 #####

cat("Find fast gathering 2... \n")

dataTable <- read.table(paste0(inputPath,"/",NAMEEE,".bam_ratio.txt"), header = TRUE)
dataTable = dataTable[,1:4]

dataTable[,1] <- gsub('X','23',dataTable[,1])

dataTable = dataTable[which(!dataTable[,1] == "Y"),] 

dataTable[,1] = as.numeric(as.character(dataTable[,1]))

dataTable = dataTable[order(dataTable[,1]),]

L = dim(dataTable)[1]

size_window = dataTable[2,2] - dataTable[1,2]

dataTable = cbind(1:L, dataTable[,1], dataTable[,2], dataTable[,2]+size_window-1, 
                  dataTable[,3], dataTable[,4]) 

dataTable = as.data.frame(dataTable)
colnames(dataTable) = c("feature", "chromosome", "start", "end", "ratio", "ratio_median")

dataTable = dataTable[which(!dataTable$ratio == -1),]
dataTable = dataTable[which(!dataTable$ratio_median == -1),]

# log2 transformation
dataTable[,5] = log2(dataTable[,5])
dataTable[,6] = log2(dataTable[,6])

ratio <- data.frame(dataTable)

short_size_window = size_window/1000

X=ratio

ratio_file_tsv = ratio # save for later

X=X[,-1] # feature
X=X[,-4] 


## Remove spurious regions based on telomere and centromere from UCSC

X <- X[which(!(X[,1]==1 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==1 & X[,2]>=121535434 & X[,3]<=124535434)),]
X <- X[which(!(X[,1]==1 & X[,2]>=249240621 & X[,3]<=249250621)),]
X <- X[which(!(X[,1]==2 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==2 & X[,2]>=92326171 & X[,3]<=95326171)),]
X <- X[which(!(X[,1]==2 & X[,2]>=243189373 & X[,3]<=243199373)),]
X <- X[which(!(X[,1]==3 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==3 & X[,2]>=90504854 & X[,3]<=93504854)),]
X <- X[which(!(X[,1]==3 & X[,2]>=198012430 & X[,3]<=198022430)),]
X <- X[which(!(X[,1]==4 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==4 & X[,2]>=49660117 & X[,3]<=52660117)),]
X <- X[which(!(X[,1]==4 & X[,2]>=191144276 & X[,3]<=191154276)),]
X <- X[which(!(X[,1]==5 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==5 & X[,2]>=46405641 & X[,3]<=9405641)),]
X <- X[which(!(X[,1]==5 & X[,2]>=180905260 & X[,3]<=180915260)),]
X <- X[which(!(X[,1]==6 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==6 & X[,2]>=58830166 & X[,3]<=61830166)),]
X <- X[which(!(X[,1]==6 & X[,2]>=171105067 & X[,3]<=171115067)),]
X <- X[which(!(X[,1]==7 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==7 & X[,2]>=58054331 & X[,3]<=61054331)),]
X <- X[which(!(X[,1]==7 & X[,2]>=159128663 & X[,3]<=159138663)),]
X <- X[which(!(X[,1]==8 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==8 & X[,2]>=43838887 & X[,3]<=46838887)),]
X <- X[which(!(X[,1]==8 & X[,2]>=146354022 & X[,3]<=146364022)),]
X <- X[which(!(X[,1]==9 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==9 & X[,2]>=47367679 & X[,3]<=50367679)),]
X <- X[which(!(X[,1]==9 & X[,2]>=141203431 & X[,3]<=141213431)),]
X <- X[which(!(X[,1]==10 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==10 & X[,2]>=39254935 & X[,3]<=42254935)),]
X <- X[which(!(X[,1]==10 & X[,2]>=135524747 & X[,3]<=135534747)),]
X <- X[which(!(X[,1]==11 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==11 & X[,2]>=51644205 & X[,3]<=54644205)),]
X <- X[which(!(X[,1]==11 & X[,2]>=134996516 & X[,3]<=135006516)),]
X <- X[which(!(X[,1]==12 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==12 & X[,2]>=34856694 & X[,3]<=37856694)),]
X <- X[which(!(X[,1]==12 & X[,2]>=133841895 & X[,3]<=133851895)),]
X <- X[which(!(X[,1]==13 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==13 & X[,2]>=16000000 & X[,3]<=19000000)),]
X <- X[which(!(X[,1]==13 & X[,2]>=115159878 & X[,3]<=115169878)),]
X <- X[which(!(X[,1]==14 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==14 & X[,2]>=16000000 & X[,3]<=19000000)),]
X <- X[which(!(X[,1]==14 & X[,2]>=107339540 & X[,3]<=107349540)),]
X <- X[which(!(X[,1]==15 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==15 & X[,2]>=17000000 & X[,3]<=20000000)),]
X <- X[which(!(X[,1]==15 & X[,2]>=102521392 & X[,3]<=102531392)),]
X <- X[which(!(X[,1]==16 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==16 & X[,2]>=35335801 & X[,3]<=38335801)),]
X <- X[which(!(X[,1]==16 & X[,2]>=90344753 & X[,3]<=90354753)),]
X <- X[which(!(X[,1]==17 & X[,2]>=22263006 & X[,3]<=25263006)),]
X <- X[which(!(X[,1]==18 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==18 & X[,2]>=15460898 & X[,3]<=18460898)),]
X <- X[which(!(X[,1]==18 & X[,2]>=78067248 & X[,3]<=78077248)),]
X <- X[which(!(X[,1]==19 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==19 & X[,2]>=24681782 & X[,3]<=27681782)),]
X <- X[which(!(X[,1]==19 & X[,2]>=59118983 & X[,3]<=59128983)),]
X <- X[which(!(X[,1]==20 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==20 & X[,2]>=26369569 & X[,3]<=29369569)),]
X <- X[which(!(X[,1]==20 & X[,2]>=63015520 & X[,3]<=63025520)),]
X <- X[which(!(X[,1]==21 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==21 & X[,2]>=11288129 & X[,3]<=14288129)),]
X <- X[which(!(X[,1]==21 & X[,2]>=48119895 & X[,3]<=48129895)),]
X <- X[which(!(X[,1]==22 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==22 & X[,2]>=13000000 & X[,3]<=16000000)),]
X <- X[which(!(X[,1]==22 & X[,2]>=51294566 & X[,3]<=51304566)),]

X <- X[which(!(X[,1]==23 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==23 & X[,2]>=58100000 & X[,3]<=63800000)),]
X <- X[which(!(X[,1]==23 & X[,2]>=155260560 & X[,3]<=155270560)),]


## Add chromosome arm 

options(show.error.messages = FALSE)

X$chr_arm <- rep(0,nrow(X))
colnames(X) <- c("chr", "start", "end", "ratio_median", "chr_arm")
X=X[c("chr", "chr_arm", "start", "end", "ratio_median")]

X[,3]=as.numeric(as.character(X[,3])) 
X[,4]=as.numeric(as.character(X[,4])) 

tt <- X[,1] == 1 & X[,3] < 125000000 & X[,4] < 125000000
X[tt,2] = 1
tt <- X[,1] == 1 & X[,3] > 125000000 & X[,4] > 125000000
X[tt,2] = 2
tt <- X[,1] == 1 & X[,3] < 125000000 & X[,4] > 125000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 125000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 125000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])  
rownames(X) <- NULL 

tt <- X[,1] == 2 & X[,3] < 93300000 & X[,4] < 93300000
X[tt,2] = 1
tt <- X[,1] == 2 & X[,3] > 93300000 & X[,4] > 93300000
X[tt,2] = 2
tt <- X[,1] == 2 & X[,3] < 93300000 & X[,4] > 93300000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 93300000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 93300001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 3 & X[,3] < 91000000 & X[,4] < 91000000
X[tt,2] = 1
tt <- X[,1] == 3 & X[,3] > 91000000 & X[,4] > 91000000
X[tt,2] = 2
tt <- X[,1] == 3 & X[,3] < 91000000 & X[,4] > 91000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 91000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 91000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 4 & X[,3] < 50400000 & X[,4] < 50400000
X[tt,2] = 1
tt <- X[,1] == 4 & X[,3] > 50400000 & X[,4] > 50400000
X[tt,2] = 2
tt <- X[,1] == 4 & X[,3] < 50400000 & X[,4] > 50400000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 50400000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 50400001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 5 & X[,3] < 48400000 & X[,4] < 48400000
X[tt,2] = 1
tt <- X[,1] == 5 & X[,3] > 48400000 & X[,4] > 48400000
X[tt,2] = 2
tt <- X[,1] == 5 & X[,3] < 48400000 & X[,4] > 48400000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 48400000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 48400001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 6 & X[,3] < 61000000 & X[,4] < 61000000
X[tt,2] = 1
tt <- X[,1] == 6 & X[,3] > 61000000 & X[,4] > 61000000
X[tt,2] = 2
tt <- X[,1] == 6 & X[,3] < 61000000 & X[,4] > 61000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 61000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 61000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 7 & X[,3] < 59900000 & X[,4] < 59900000
X[tt,2] = 1
tt <- X[,1] == 7 & X[,3] > 59900000 & X[,4] > 59900000
X[tt,2] = 2
tt <- X[,1] == 7 & X[,3] < 59900000 & X[,4] > 59900000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 59900000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 59900001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 8 & X[,3] < 45600000 & X[,4] < 45600000
X[tt,2] = 1
tt <- X[,1] == 8 & X[,3] > 45600000 & X[,4] > 45600000
X[tt,2] = 2
tt <- X[,1] == 8 & X[,3] < 45600000 & X[,4] > 45600000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 45600000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 45600001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 9 & X[,3] < 49000000 & X[,4] < 49000000
X[tt,2] = 1
tt <- X[,1] == 9 & X[,4] > 49000000 & X[,4] > 49000000
X[tt,2] = 2
tt <- X[,1] == 9 & X[,3] < 49000000 & X[,4] > 49000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 49000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 49000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 10 & X[,3] < 40200000 & X[,4] < 40200000
X[tt,2] = 1
tt <- X[,1] == 10 & X[,3] > 40200000 & X[,4] > 40200000
X[tt,2] = 2
tt <- X[,1] == 10 & X[,3] < 40200000 & X[,4] > 40200000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 40200000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 40200001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 11 & X[,3] < 53700000 & X[,4] < 53700000
X[tt,2] = 1
tt <- X[,1] == 11 & X[,3] > 53700000 & X[,4] > 53700000
X[tt,2] = 2
tt <- X[,1] == 11 & X[,3] < 53700000 & X[,4] > 53700000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 53700000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 53700001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 12 & X[,3] < 35800000 & X[,4] < 35800000
X[tt,2] = 1
tt <- X[,1] == 12 & X[,3] > 35800000 & X[,4] > 35800000
X[tt,2] = 2
tt <- X[,1] == 12 & X[,3] < 35800000 & X[,4] > 35800000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 35800000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 35800001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 13 & X[,3] < 17900000 & X[,4] < 17900000
X[tt,2] = 1
tt <- X[,1] == 13 & X[,3] > 17900000 & X[,4] > 17900000
X[tt,2] = 2
tt <- X[,1] == 13 & X[,3] < 17900000 & X[,4] > 17900000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 17900000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 17900001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 14 & X[,3] < 17600000 & X[,4] < 17600000
X[tt,2] = 1
tt <- X[,1] == 14 & X[,3] > 17600000 & X[,4] > 17600000
X[tt,2] = 2
tt <- X[,1] == 14 & X[,3] < 17600000 & X[,4] > 17600000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 17600000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 17600001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 15 & X[,3] < 19000000 & X[,4] < 19000000
X[tt,2] = 1
tt <- X[,1] == 15 & X[,3] > 19000000 & X[,4] > 19000000
X[tt,2] = 2
tt <- X[,1] == 15 & X[,3] < 19000000 & X[,4] > 19000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 19000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 19000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 16 & X[,3] < 36600000 & X[,4] < 36600000
X[tt,2] = 1
tt <- X[,1] == 16 & X[,3] > 36600000 & X[,4] > 36600000
X[tt,2] = 2
tt <- X[,1] == 16 & X[,3] < 36600000 & X[,4] > 36600000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 36600000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 36600001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 17 & X[,3] < 24000000 & X[,4] < 24000000
X[tt,2] = 1
tt <- X[,1] == 17 & X[,3] > 24000000 & X[,4] > 24000000
X[tt,2] = 2
tt <- X[,1] == 17 & X[,3] < 24000000 & X[,4] > 24000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 24000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 24000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 18 & X[,3] < 17200000 & X[,4] < 17200000
X[tt,2] = 1
tt <- X[,1] == 18 & X[,3] > 17200000 & X[,4] > 17200000
X[tt,2] = 2
tt <- X[,1] == 18 & X[,3] < 17200000 & X[,4] > 17200000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 17200000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 17200001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 19 & X[,3] < 26500000 & X[,4] < 26500000
X[tt,2] = 1
tt <- X[,1] == 19 & X[,3] > 26500000 & X[,4] > 26500000
X[tt,2] = 2
tt <- X[,1] == 19 & X[,3] < 26500000 & X[,4] > 26500000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 26500000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 26500001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 20 & X[,3] < 27500000 & X[,4] < 27500000
X[tt,2] = 1
tt <- X[,1] == 20 & X[,3] > 27500000 & X[,4] > 27500000
X[tt,2] = 2
tt <- X[,1] == 20 & X[,3] < 27500000 & X[,4] > 27500000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 27500000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 27500001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 21 & X[,3] < 13200000 & X[,4] < 13200000
X[tt,2] = 1
tt <- X[,1] == 21 & X[,3] > 13200000 & X[,4] > 13200000
X[tt,2] = 2
tt <- X[,1] == 21 & X[,3] < 13200000 & X[,4] > 13200000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 13200000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 13200001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 22 & X[,3] < 14700000 & X[,4] < 14700000
X[tt,2] = 1
tt <- X[,1] == 22 & X[,3] > 14700000 & X[,4] > 14700000
X[tt,2] = 2
tt <- X[,1] == 22 & X[,3] < 14700000 & X[,4] > 14700000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 14700000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 14700001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 23 & X[,3] < 60950000 & X[,4] < 60950000
X[tt,2] = 1
tt <- X[,1] == 23 & X[,3] > 60950000 & X[,4] > 60950000
X[tt,2] = 2
tt <- X[,1] == 23 & X[,3] < 60950000 & X[,4] > 60950000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 60950000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 60950000, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

X[,2]<-X[,1]+(X[,2]-1)/2      

tt=which(X$ratio_median==-Inf)

if (length(tt) >= 1){
  X=X[-tt,]  
}

X=as.matrix(X)

L = dim(X)[1]
A = matrix(0, ncol=6, nrow=L)

i=1    
c=1    

while (i < L){
  if(X[i,2]==X[i+1,2]){             
    if (X[i,5]==X[i+1,5]){          
      n=1
      while (X[i,2]==X[i+n,2] & X[i,5]==X[i+n,5] & i+n < L){
        n=n+1
      }
      if(i+n == L){
        A[c,]=c(X[i,1], X[i,2], X[i,3], X[i+n,4], X[i,5], X[i+n,4]-X[i,3]+1) 
        i=i+n
        c=c+1
      }
      else{
        A[c,]=c(X[i,1], X[i,2], X[i,3], X[i+n-1,4], X[i,5], X[i+n-1,4]-X[i,3]+1) 
        i=i+n
        c=c+1
      }
    }
    else{
      A[c,]=c(X[i,1], X[i,2], X[i,3], X[i,4], X[i,5], X[i,4]-X[i,3]+1)
      i=i+1
      c=c+1
    }
  }
  else{
    A[c,]=c(X[i,1], X[i,2], X[i,3], X[i,4], X[i,5], X[i,4]-X[i,3]+1)
    i=i+1                                        
    c=c+1                         
  }  
}

A=subset(A, A[,1] != 0)
rownames(A) <- NULL 
colnames(A) <- c("chr", "chr_arm", "start", "end", "ratio_median", "size")


write.table(A, file = paste0(outputPath,"/",NAMEEE,"_ratio_median_gathered.txt"), sep = "\t", row.names = FALSE)

options(show.error.messages = TRUE)

