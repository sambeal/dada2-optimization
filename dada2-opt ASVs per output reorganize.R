#ASV analysis of SINE dada2 optimization
#1 <-read.csv("dada2-optimization/organized_ASVs/output1_ASVs.csv")
#2 <-read.csv("dada2-optimization/organized_ASVs/output2_ASVs.csv")
#3 <-read.csv("dada2-optimization/organized_ASVs/output3_ASVs.csv")
#4 <-read.csv("dada2-optimization/organized_ASVs/output4_ASVs.csv")
#5 <-read.csv("dada2-optimization/organized_ASVs/output5_ASVs.csv")
#6 <-read.csv("dada2-optimization/organized_ASVs/output6_ASVs.csv")
#7 <-read.csv("dada2-optimization/organized_ASVs/output7_ASVs.csv")
#8 <-read.csv("dada2-optimization/organized_ASVs/output8_ASVs.csv")
#9 <-read.csv("dada2-optimization/organized_ASVs/output9_ASVs.csv")


#pushed after output-9 on 23.03.28


#setwd("/Users/samanthabeal/Documents/MSc/Bioinformatics")
#getwd()
#setwd("..")

#load data----
output9 <-read.csv("dada2-optimization/output9/output9_ASVs.csv")

#rename 1st dada column
names(output9)[1] <- "Sequence"
#add names to the column - will need this for the haplotype network
output9$ASV <- 1:nrow(output9)

#move to first location
library(dplyr)
output9 <- output9 %>% relocate(ASV, .before = Sequence)

#add 'seq' to each seq number 
output9$ASV <- sub("^", "seq", output9$ASV)

#rename columns to only be fish and index
names(output9) = gsub(pattern = "X", replacement = "", x = names(output9))

#pull out individual fish----
#190222-A1----
A1190222o9 <- output9[c(1:5)]
names(A1190222o9)[3] <- "PCR.1"
names(A1190222o9)[4] <- "PCR.2"
names(A1190222o9)[5] <- "PCR.3"

#determine if the seq is detected in all 3 reps
A1190222o9$Consistently_detected <- ifelse(A1190222o9$PCR.1==0 | A1190222o9$PCR.2== 0 | A1190222o9$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
A1190222o9_true <- subset(A1190222o9, Consistently_detected == "true",
                          select=c("ASV","Sequence","PCR.1","PCR.2","PCR.3"))
A1190222o9_true$Fish <- "A1190222"
A1190222o9_true$Output <- "output9"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
A1190222o9sum <- as.data.frame(colSums(A1190222o9[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A1190222o9sum <- as.data.frame(t(A1190222o9sum))

#add column of avg depth across all PCR reps
A1190222o9sum$Avg.depth <-apply(A1190222o9sum,1,mean)
#add columns of total depth
A1190222o9sum$Total.depth <-apply(A1190222o9sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
A1190222o9sum$Num.ASVs <- sum(A1190222o9$PCR.1>0 | A1190222o9$PCR.2>0 | A1190222o9$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
A1190222o9sum$Output <- "output9"

#add column of if in all 3 reps
A1190222o9sum$All.reps <- "No"

#only true ASVs
A1190222o9_truesum <- as.data.frame(colSums(A1190222o9_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A1190222o9_truesum <- as.data.frame(t(A1190222o9_truesum))

#add column of avg depth across all PCR reps
A1190222o9_truesum $Avg.depth <-apply(A1190222o9_truesum,1,mean)

#add columns of total depth
A1190222o9_truesum $Total.depth <-apply(A1190222o9_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
A1190222o9_truesum$Num.ASVs <- sum(A1190222o9_true$PCR.1>0 | A1190222o9_true$PCR.2>0 | A1190222o9_true$PCR.3>0)

#add column of output (so can plot with others for comparison)
A1190222o9_truesum$Output <- "output9"

#add column of if in all 3 reps
A1190222o9_truesum$All.reps <- "Yes"

#add data frames together
A1190222_output9 <- rbind(A1190222o9_truesum, A1190222o9sum)

#add column of fish ID
A1190222_output9$Fish <- "190922-A1"

#190222-A2----
A2190222o9 <- output9[c(1:2,6:8)]
names(A2190222o9)[3] <- "PCR.1"
names(A2190222o9)[4] <- "PCR.2"
names(A2190222o9)[5] <- "PCR.3"

#determine if the seq is detected in all 3 reps
A2190222o9$Consistently_detected <- ifelse(A2190222o9$PCR.1==0 | A2190222o9$PCR.2== 0 | A2190222o9$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
A2190222o9_true <- subset(A2190222o9, Consistently_detected == "true",
                          select=c("ASV","Sequence","PCR.1","PCR.2","PCR.3"))
A2190222o9_true$Fish <- "A2190222"
A2190222o9_true$Output <- "output9"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
A2190222o9sum <- as.data.frame(colSums(A2190222o9[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A2190222o9sum <- as.data.frame(t(A2190222o9sum))

#add column of avg depth across all PCR reps
A2190222o9sum$Avg.depth <-apply(A2190222o9sum,1,mean)
#add columns of total depth
A2190222o9sum$Total.depth <-apply(A2190222o9sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
A2190222o9sum$Num.ASVs <- sum(A2190222o9$PCR.1>0 | A2190222o9$PCR.2>0 | A2190222o9$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
A2190222o9sum$Output <- "output9"

#add column of if in all 3 reps
A2190222o9sum$All.reps <- "No"

#only true ASVs
A2190222o9_truesum <- as.data.frame(colSums(A2190222o9_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A2190222o9_truesum <- as.data.frame(t(A2190222o9_truesum))

#add column of avg depth across all PCR reps
A2190222o9_truesum $Avg.depth <-apply(A2190222o9_truesum,1,mean)

#add columns of total depth
A2190222o9_truesum $Total.depth <-apply(A2190222o9_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
A2190222o9_truesum$Num.ASVs <- sum(A2190222o9_true$PCR.1>0 | A2190222o9_true$PCR.2>0 | A2190222o9_true$PCR.3>0)

#add column of pipeline (so can plot with different outputs for comparison)
A2190222o9_truesum$Output <- "output9"

#add column of if in all 3 reps
A2190222o9_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
A2190222_output9 <- rbind(A2190222o9_truesum, A2190222o9sum)

#add column of fish ID
A2190222_output9$Fish <- "190922-A2"

#190222-H1----
H1190222o9 <- output9[c(1:2,9:11)]
names(H1190222o9)[3] <- "PCR.1"
names(H1190222o9)[4] <- "PCR.2"
names(H1190222o9)[5] <- "PCR.3"

#determine if the seq is detected in all 3 reps
H1190222o9$Consistently_detected <- ifelse(H1190222o9$PCR.1==0 | H1190222o9$PCR.2== 0 | H1190222o9$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
H1190222o9_true <- subset(H1190222o9, Consistently_detected == "true",
                          select=c("ASV","Sequence","PCR.1","PCR.2","PCR.3"))
H1190222o9_true$Fish <- "H1190222"
H1190222o9_true$Output <- "output9"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
H1190222o9sum <- as.data.frame(colSums(H1190222o9[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
H1190222o9sum <- as.data.frame(t(H1190222o9sum))

#add column of avg depth across all PCR reps
H1190222o9sum$Avg.depth <-apply(H1190222o9sum,1,mean)
#add columns of total depth
H1190222o9sum$Total.depth <-apply(H1190222o9sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
H1190222o9sum$Num.ASVs <- sum(H1190222o9$PCR.1>0 | H1190222o9$PCR.2>0 | H1190222o9$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
H1190222o9sum$Output <- "output9"

#add column of if in all 3 reps
H1190222o9sum$All.reps <- "No"

#only true ASVs
H1190222o9_truesum <- as.data.frame(colSums(H1190222o9_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
H1190222o9_truesum <- as.data.frame(t(H1190222o9_truesum))

#add column of avg depth across all PCR reps
H1190222o9_truesum $Avg.depth <-apply(H1190222o9_truesum,1,mean)

#add columns of total depth
H1190222o9_truesum $Total.depth <-apply(H1190222o9_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
H1190222o9_truesum$Num.ASVs <- sum(H1190222o9_true$PCR.1>0 | H1190222o9_true$PCR.2>0 | H1190222o9_true$PCR.3>0)

#add column of pipeline (so can plot with different outputs for comparison)
H1190222o9_truesum$Output <- "output9"

#add column of if in all 3 reps
H1190222o9_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
H1190222_output9 <- rbind(H1190222o9_truesum, H1190222o9sum)

#add column of fish ID
H1190222_output9$Fish <- "190922-H1"

#190222-M1----
M1190222o9 <- output9[c(1:2,12:14)]
names(M1190222o9)[3] <- "PCR.1"
names(M1190222o9)[4] <- "PCR.2"
names(M1190222o9)[5] <- "PCR.3"

#determine if the seq is detected in all 3 reps
M1190222o9$Consistently_detected <- ifelse(M1190222o9$PCR.1==0 | M1190222o9$PCR.2== 0 | M1190222o9$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
M1190222o9_true <- subset(M1190222o9, Consistently_detected == "true",
                          select=c("ASV","Sequence","PCR.1","PCR.2","PCR.3"))
M1190222o9_true$Fish <- "M1190222"
M1190222o9_true$Output <- "output9"


#compare size of two dfs (# ASVs and total reads)
#all ASVs
M1190222o9sum <- as.data.frame(colSums(M1190222o9[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
M1190222o9sum <- as.data.frame(t(M1190222o9sum))

#add column of avg depth across all PCR reps
M1190222o9sum$Avg.depth <-apply(M1190222o9sum,1,mean)
#add columns of total depth
M1190222o9sum$Total.depth <-apply(M1190222o9sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
M1190222o9sum$Num.ASVs <- sum(M1190222o9$PCR.1>0 | M1190222o9$PCR.2>0 | M1190222o9$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
M1190222o9sum$Output <- "output9"

#add column of if in all 3 reps
M1190222o9sum$All.reps <- "No"

#only true ASVs
M1190222o9_truesum <- as.data.frame(colSums(M1190222o9_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
M1190222o9_truesum <- as.data.frame(t(M1190222o9_truesum))

#add column of avg depth across all PCR reps
M1190222o9_truesum $Avg.depth <-apply(M1190222o9_truesum,1,mean)

#add columns of total depth
M1190222o9_truesum $Total.depth <-apply(M1190222o9_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
M1190222o9_truesum$Num.ASVs <- sum(M1190222o9_true$PCR.1>0 | M1190222o9_true$PCR.2>0 | M1190222o9_true$PCR.3>0)

#add column of pipeline (so can plot with different outputs for comparison)
M1190222o9_truesum$Output <- "output9"

#add column of if in all 3 reps
M1190222o9_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
M1190222_output9 <- rbind(M1190222o9_truesum, M1190222o9sum)

#add column of fish ID
M1190222_output9$Fish <- "190922-M1"

#190918-A1----
A1190918o9 <- output9[c(1:2,15:17)]
names(A1190918o9)[3] <- "PCR.1"
names(A1190918o9)[4] <- "PCR.2"
names(A1190918o9)[5] <- "PCR.3"

#determine if the seq is detected in all 3 reps
A1190918o9$Consistently_detected <- ifelse(A1190918o9$PCR.1==0 | A1190918o9$PCR.2== 0 | A1190918o9$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
A1190918o9_true <- subset(A1190918o9, Consistently_detected == "true",
                          select=c("ASV","Sequence","PCR.1","PCR.2","PCR.3"))
A1190918o9_true$Fish <- "A1190918"
A1190918o9_true$Output <- "output9"
#compare size of two dfs (# ASVs and total reads)
#all ASVs
A1190918o9sum <- as.data.frame(colSums(A1190918o9[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A1190918o9sum <- as.data.frame(t(A1190918o9sum))

#add column of avg depth across all PCR reps
A1190918o9sum$Avg.depth <-apply(A1190918o9sum,1,mean)
#add columns of total depth
A1190918o9sum$Total.depth <-apply(A1190918o9sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
A1190918o9sum$Num.ASVs <- sum(A1190918o9$PCR.1>0 | A1190918o9$PCR.2>0 | A1190918o9$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
A1190918o9sum$Output <- "output9"

#add column of if in all 3 reps
A1190918o9sum$All.reps <- "No"

#only true ASVs
A1190918o9_truesum <- as.data.frame(colSums(A1190918o9_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A1190918o9_truesum <- as.data.frame(t(A1190918o9_truesum))

#add column of avg depth across all PCR reps
A1190918o9_truesum $Avg.depth <-apply(A1190918o9_truesum,1,mean)

#add columns of total depth
A1190918o9_truesum $Total.depth <-apply(A1190918o9_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
A1190918o9_truesum$Num.ASVs <- sum(A1190918o9_true$PCR.1>0 | A1190918o9_true$PCR.2>0 | A1190918o9_true$PCR.3>0)

#add column of pipeline (so can plot with different outputs for comparison)
A1190918o9_truesum$Output <- "output9"

#add column of if in all 3 reps
A1190918o9_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
A1190918_output9 <- rbind(A1190918o9_truesum, A1190918o9sum)

#add column of fish ID
A1190918_output9$Fish <- "190918-A1"

#190918-A3----
A3190918o9 <- output9[c(1:2,18:20)]
names(A3190918o9)[3] <- "PCR.1"
names(A3190918o9)[4] <- "PCR.2"
names(A3190918o9)[5] <- "PCR.3"

#determine if the seq is detected in all 3 reps
A3190918o9$Consistently_detected <- ifelse(A3190918o9$PCR.1==0 | A3190918o9$PCR.2== 0 | A3190918o9$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
A3190918o9_true <- subset(A3190918o9, Consistently_detected == "true",
                          select=c("ASV","Sequence","PCR.1","PCR.2","PCR.3"))
A3190918o9_true$Fish <- "A3190918"
A3190918o9_true$Output <- "output9"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
A3190918o9sum <- as.data.frame(colSums(A3190918o9[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A3190918o9sum <- as.data.frame(t(A3190918o9sum))

#add column of avg depth across all PCR reps
A3190918o9sum$Avg.depth <-apply(A3190918o9sum,1,mean)
#add columns of total depth
A3190918o9sum$Total.depth <-apply(A3190918o9sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
A3190918o9sum$Num.ASVs <- sum(A3190918o9$PCR.1>0 | A3190918o9$PCR.2>0 | A3190918o9$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
A3190918o9sum$Output <- "output9"

#add column of if in all 3 reps
A3190918o9sum$All.reps <- "No"

#only true ASVs
A3190918o9_truesum <- as.data.frame(colSums(A3190918o9_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A3190918o9_truesum <- as.data.frame(t(A3190918o9_truesum))

#add column of avg depth across all PCR reps
A3190918o9_truesum $Avg.depth <-apply(A3190918o9_truesum,1,mean)

#add columns of total depth
A3190918o9_truesum $Total.depth <-apply(A3190918o9_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
A3190918o9_truesum$Num.ASVs <- sum(A3190918o9_true$PCR.1>0 | A3190918o9_true$PCR.2>0 | A3190918o9_true$PCR.3>0)

#add column of pipeline (so can plot with different outputs for comparison)
A3190918o9_truesum$Output <- "output9"

#add column of if in all 3 reps
A3190918o9_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
A3190918_output9 <- rbind(A3190918o9_truesum, A3190918o9sum)

#add column of fish ID
A3190918_output9$Fish <- "190918-A3"

#190918-A4----
A4190918o9 <- output9[c(1:2,21:23)]
names(A4190918o9)[3] <- "PCR.1"
names(A4190918o9)[4] <- "PCR.2"
names(A4190918o9)[5] <- "PCR.3"

#determine if the seq is detected in all 3 reps
A4190918o9$Consistently_detected <- ifelse(A4190918o9$PCR.1==0 | A4190918o9$PCR.2== 0 | A4190918o9$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
A4190918o9_true <- subset(A4190918o9, Consistently_detected == "true",
                          select=c("ASV","Sequence","PCR.1","PCR.2","PCR.3"))
A4190918o9_true$Fish <- "A4190918"
A4190918o9_true$Output <- "output9"


#compare size of two dfs (# ASVs and total reads)
#all ASVs
A4190918o9sum <- as.data.frame(colSums(A4190918o9[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A4190918o9sum <- as.data.frame(t(A4190918o9sum))

#add column of avg depth across all PCR reps
A4190918o9sum$Avg.depth <-apply(A4190918o9sum,1,mean)
#add columns of total depth
A4190918o9sum$Total.depth <-apply(A4190918o9sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
A4190918o9sum$Num.ASVs <- sum(A4190918o9$PCR.1>0 | A4190918o9$PCR.2>0 | A4190918o9$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
A4190918o9sum$Output <- "output9"

#add column of if in all 3 reps
A4190918o9sum$All.reps <- "No"

#only true ASVs
A4190918o9_truesum <- as.data.frame(colSums(A4190918o9_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A4190918o9_truesum <- as.data.frame(t(A4190918o9_truesum))

#add column of avg depth across all PCR reps
A4190918o9_truesum $Avg.depth <-apply(A4190918o9_truesum,1,mean)

#add columns of total depth
A4190918o9_truesum $Total.depth <-apply(A4190918o9_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
A4190918o9_truesum$Num.ASVs <- sum(A4190918o9_true$PCR.1>0 | A4190918o9_true$PCR.2>0 | A4190918o9_true$PCR.3>0)

#add column of pipeline (so can plot with different outputs for comparison)
A4190918o9_truesum$Output <- "output9"

#add column of if in all 3 reps
A4190918o9_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
A4190918_output9 <- rbind(A4190918o9_truesum, A4190918o9sum)

#add column of fish ID
A4190918_output9$Fish <- "190918-A4"

#190918-C4----
C4190918o9 <- output9[c(1:2,24:26)]
names(C4190918o9)[3] <- "PCR.1"
names(C4190918o9)[4] <- "PCR.2"
names(C4190918o9)[5] <- "PCR.3"

#determine if the seq is detected in all 3 reps
C4190918o9$Consistently_detected <- ifelse(C4190918o9$PCR.1==0 | C4190918o9$PCR.2== 0 | C4190918o9$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
C4190918o9_true <- subset(C4190918o9, Consistently_detected == "true",
                          select=c("ASV","Sequence","PCR.1","PCR.2","PCR.3"))
C4190918o9_true$Fish <- "C4190918"
C4190918o9_true$Output <- "output9"


#compare size of two dfs (# ASVs and total reads)
#all ASVs
C4190918o9sum <- as.data.frame(colSums(C4190918o9[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
C4190918o9sum <- as.data.frame(t(C4190918o9sum))

#add column of avg depth across all PCR reps
C4190918o9sum$Avg.depth <-apply(C4190918o9sum,1,mean)
#add columns of total depth
C4190918o9sum$Total.depth <-apply(C4190918o9sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
C4190918o9sum$Num.ASVs <- sum(C4190918o9$PCR.1>0 | C4190918o9$PCR.2>0 | C4190918o9$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
C4190918o9sum$Output <- "output9"

#add column of if in all 3 reps
C4190918o9sum$All.reps <- "No"

#only true ASVs
C4190918o9_truesum <- as.data.frame(colSums(C4190918o9_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
C4190918o9_truesum <- as.data.frame(t(C4190918o9_truesum))

#add column of avg depth across all PCR reps
C4190918o9_truesum $Avg.depth <-apply(C4190918o9_truesum,1,mean)

#add columns of total depth
C4190918o9_truesum $Total.depth <-apply(C4190918o9_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
C4190918o9_truesum$Num.ASVs <- sum(C4190918o9_true$PCR.1>0 | C4190918o9_true$PCR.2>0 | C4190918o9_true$PCR.3>0)

#add column of pipeline (so can plot with different outputs for comparison)
C4190918o9_truesum$Output <- "output9"

#add column of if in all 3 reps
C4190918o9_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
C4190918_output9 <- rbind(C4190918o9_truesum, C4190918o9sum)

#add column of fish ID
C4190918_output9$Fish <- "190918-C4"

#190918-D4----
D4190918o9 <- output9[c(1:2,27:29)]
names(D4190918o9)[3] <- "PCR.1"
names(D4190918o9)[4] <- "PCR.2"
names(D4190918o9)[5] <- "PCR.3"

#determine if the seq is detected in all 3 reps
D4190918o9$Consistently_detected <- ifelse(D4190918o9$PCR.1==0 | D4190918o9$PCR.2== 0 | D4190918o9$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
D4190918o9_true <- subset(D4190918o9, Consistently_detected == "true",
                          select=c("ASV","Sequence","PCR.1","PCR.2","PCR.3"))
D4190918o9_true$Fish <- "D4190918"
D4190918o9_true$Output <- "output9"


#compare size of two dfs (# ASVs and total reads)
#all ASVs
D4190918o9sum <- as.data.frame(colSums(D4190918o9[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
D4190918o9sum <- as.data.frame(t(D4190918o9sum))

#add column of avg depth across all PCR reps
D4190918o9sum$Avg.depth <-apply(D4190918o9sum,1,mean)
#add columns of total depth
D4190918o9sum$Total.depth <-apply(D4190918o9sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
D4190918o9sum$Num.ASVs <- sum(D4190918o9$PCR.1>0 | D4190918o9$PCR.2>0 | D4190918o9$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
D4190918o9sum$Output <- "output9"

#add column of if in all 3 reps
D4190918o9sum$All.reps <- "No"

#only true ASVs
D4190918o9_truesum <- as.data.frame(colSums(D4190918o9_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
D4190918o9_truesum <- as.data.frame(t(D4190918o9_truesum))

#add column of avg depth across all PCR reps
D4190918o9_truesum$Avg.depth <-apply(D4190918o9_truesum,1,mean)

#add columns of total depth
D4190918o9_truesum$Total.depth <-apply(D4190918o9_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
D4190918o9_truesum$Num.ASVs <- sum(D4190918o9_true$PCR.1>0 | D4190918o9_true$PCR.2>0 | D4190918o9_true$PCR.3>0)

#add column of output (so can plot with different outputs for comparison)
D4190918o9_truesum$Output <- "output9"

#add column of if in all 3 reps
D4190918o9_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
D4190918_output9 <- rbind(D4190918o9_truesum, D4190918o9sum)

#add column of fish ID
D4190918_output9$Fish <- "190918-D4"

#190918-E7----
E7190918o9 <- output9[c(1:2,30:32)]
names(E7190918o9)[3] <- "PCR.1"
names(E7190918o9)[4] <- "PCR.2"
names(E7190918o9)[5] <- "PCR.3"

#determine if the seq is detected in all 3 reps
E7190918o9$Consistently_detected <- ifelse(E7190918o9$PCR.1==0 | E7190918o9$PCR.2== 0 | E7190918o9$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
E7190918o9_true <- subset(E7190918o9, Consistently_detected == "true",
                          select=c("ASV","Sequence","PCR.1","PCR.2","PCR.3"))
E7190918o9_true$Fish <- "E7190918"
E7190918o9_true$Output <- "output9"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
E7190918o9sum <- as.data.frame(colSums(E7190918o9[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
E7190918o9sum <- as.data.frame(t(E7190918o9sum))

#add column of avg depth across all PCR reps
E7190918o9sum$Avg.depth <-apply(E7190918o9sum,1,mean)
#add columns of total depth
E7190918o9sum$Total.depth <-apply(E7190918o9sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
E7190918o9sum$Num.ASVs <- sum(E7190918o9$PCR.1>0 | E7190918o9$PCR.2>0 | E7190918o9$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
E7190918o9sum$Output <- "output9"

#add column of if in all 3 reps
E7190918o9sum$All.reps <- "No"

#only true ASVs
E7190918o9_truesum <- as.data.frame(colSums(E7190918o9_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
E7190918o9_truesum <- as.data.frame(t(E7190918o9_truesum))

#add column of avg depth across all PCR reps
E7190918o9_truesum $Avg.depth <-apply(E7190918o9_truesum,1,mean)

#add columns of total depth
E7190918o9_truesum $Total.depth <-apply(E7190918o9_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
E7190918o9_truesum$Num.ASVs <- sum(E7190918o9_true$PCR.1>0 | E7190918o9_true$PCR.2>0 | E7190918o9_true$PCR.3>0)

#add column of pipeline (so can plot with different outputs for comparison)
E7190918o9_truesum$Output <- "output9"

#add column of if in all 3 reps
E7190918o9_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
E7190918_output9 <- rbind(E7190918o9_truesum, E7190918o9sum)

#add column of fish ID
E7190918_output9$Fish <- "190918-E7"

#190918-F2----
F2190918o9 <- output9[c(1:2,33:35)]
names(F2190918o9)[3] <- "PCR.1"
names(F2190918o9)[4] <- "PCR.2"
names(F2190918o9)[5] <- "PCR.3"

#determine if the seq is detected in all 3 reps
F2190918o9$Consistently_detected <- ifelse(F2190918o9$PCR.1==0 | F2190918o9$PCR.2== 0 | F2190918o9$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
F2190918o9_true <- subset(F2190918o9, Consistently_detected == "true",
                          select=c("ASV","Sequence","PCR.1","PCR.2","PCR.3"))
F2190918o9_true$Fish <- "F2190918"
F2190918o9_true$Output <- "output9"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
F2190918o9sum <- as.data.frame(colSums(F2190918o9[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
F2190918o9sum <- as.data.frame(t(F2190918o9sum))

#add column of avg depth across all PCR reps
F2190918o9sum$Avg.depth <-apply(F2190918o9sum,1,mean)
#add columns of total depth
F2190918o9sum$Total.depth <-apply(F2190918o9sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
F2190918o9sum$Num.ASVs <- sum(F2190918o9$PCR.1>0 | F2190918o9$PCR.2>0 | F2190918o9$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
F2190918o9sum$Output <- "output9"

#add column of if in all 3 reps
F2190918o9sum$All.reps <- "No"

#only true ASVs
F2190918o9_truesum <- as.data.frame(colSums(F2190918o9_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
F2190918o9_truesum <- as.data.frame(t(F2190918o9_truesum))

#add column of avg depth across all PCR reps
F2190918o9_truesum $Avg.depth <-apply(F2190918o9_truesum,1,mean)

#add columns of total depth
F2190918o9_truesum $Total.depth <-apply(F2190918o9_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
F2190918o9_truesum$Num.ASVs <- sum(F2190918o9_true$PCR.1>0 | F2190918o9_true$PCR.2>0 | F2190918o9_true$PCR.3>0)

#add column of pipeline (so can plot with different outputs for comparison)
F2190918o9_truesum$Output <- "output9"

#add column of if in all 3 reps
F2190918o9_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
F2190918_output9 <- rbind(F2190918o9_truesum, F2190918o9sum)

#add column of fish ID
F2190918_output9$Fish <- "190918-F2"

#190922-D2----
D2190922o9 <- output9[c(1:2,36:38)]
names(D2190922o9)[3] <- "PCR.1"
names(D2190922o9)[4] <- "PCR.2"
names(D2190922o9)[5] <- "PCR.3"

#determine if the seq is detected in all 3 reps
D2190922o9$Consistently_detected <- ifelse(D2190922o9$PCR.1==0 | D2190922o9$PCR.2== 0 | D2190922o9$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
D2190922o9_true <- subset(D2190922o9, Consistently_detected == "true",
                          select=c("ASV","Sequence","PCR.1","PCR.2","PCR.3"))
D2190922o9_true$Fish <- "D2190922"
D2190922o9_true$Output <- "output9"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
D2190922o9sum <- as.data.frame(colSums(D2190922o9[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
D2190922o9sum <- as.data.frame(t(D2190922o9sum))

#add column of avg depth across all PCR reps
D2190922o9sum$Avg.depth <-apply(D2190922o9sum,1,mean)
#add columns of total depth
D2190922o9sum$Total.depth <-apply(D2190922o9sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
D2190922o9sum$Num.ASVs <- sum(D2190922o9$PCR.1>0 | D2190922o9$PCR.2>0 | D2190922o9$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
D2190922o9sum$Output <- "output9"

#add column of if in all 3 reps
D2190922o9sum$All.reps <- "No"

#only true ASVs
D2190922o9_truesum <- as.data.frame(colSums(D2190922o9_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
D2190922o9_truesum <- as.data.frame(t(D2190922o9_truesum))

#add column of avg depth across all PCR reps
D2190922o9_truesum $Avg.depth <-apply(D2190922o9_truesum,1,mean)

#add columns of total depth
D2190922o9_truesum $Total.depth <-apply(D2190922o9_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
D2190922o9_truesum$Num.ASVs <- sum(D2190922o9_true$PCR.1>0 | D2190922o9_true$PCR.2>0 | D2190922o9_true$PCR.3>0)

#add column of pipeline (so can plot with different outputs for comparison)
D2190922o9_truesum$Output <- "output9"

#add column of if in all 3 reps
D2190922o9_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
D2190922_output9 <- rbind(D2190922o9_truesum, D2190922o9sum)

#add column of fish ID
D2190922_output9$Fish <- "190922-D2"

#190922-F2----
F2190922o9 <- output9[c(1:2,39:41)]
names(F2190922o9)[3] <- "PCR.1"
names(F2190922o9)[4] <- "PCR.2"
names(F2190922o9)[5] <- "PCR.3"

#determine if the seq is detected in all 3 reps
F2190922o9$Consistently_detected <- ifelse(F2190922o9$PCR.1==0 | F2190922o9$PCR.2== 0 | F2190922o9$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
F2190922o9_true <- subset(F2190922o9, Consistently_detected == "true",
                          select=c("ASV","Sequence","PCR.1","PCR.2","PCR.3"))
F2190922o9_true$Fish <- "F2190922"
F2190922o9_true$Output <- "output9"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
F2190922o9sum <- as.data.frame(colSums(F2190922o9[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
F2190922o9sum <- as.data.frame(t(F2190922o9sum))

#add column of avg depth across all PCR reps
F2190922o9sum$Avg.depth <-apply(F2190922o9sum,1,mean)
#add columns of total depth
F2190922o9sum$Total.depth <-apply(F2190922o9sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
F2190922o9sum$Num.ASVs <- sum(F2190922o9$PCR.1>0 | F2190922o9$PCR.2>0 | F2190922o9$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
F2190922o9sum$Output <- "output9"

#add column of if in all 3 reps
F2190922o9sum$All.reps <- "No"

#only true ASVs
F2190922o9_truesum <- as.data.frame(colSums(F2190922o9_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
F2190922o9_truesum <- as.data.frame(t(F2190922o9_truesum))

#add column of avg depth across all PCR reps
F2190922o9_truesum $Avg.depth <-apply(F2190922o9_truesum,1,mean)

#add columns of total depth
F2190922o9_truesum $Total.depth <-apply(F2190922o9_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
F2190922o9_truesum$Num.ASVs <- sum(F2190922o9_true$PCR.1>0 | F2190922o9_true$PCR.2>0 | F2190922o9_true$PCR.3>0)

#add column of pipeline (so can plot with different outputs for comparison)
F2190922o9_truesum$Output <- "output9"

#add column of if in all 3 reps
F2190922o9_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
F2190922_output9 <- rbind(F2190922o9_truesum, F2190922o9sum)

#add column of fish ID
F2190922_output9$Fish <- "190922-F2"

#190922-G2----
G2190922o9 <- output9[c(1:2,42:44)]
names(G2190922o9)[3] <- "PCR.1"
names(G2190922o9)[4] <- "PCR.2"
names(G2190922o9)[5] <- "PCR.3"

#determine if the seq is detected in all 3 reps
G2190922o9$Consistently_detected <- ifelse(G2190922o9$PCR.1==0 | G2190922o9$PCR.2== 0 | G2190922o9$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
G2190922o9_true <- subset(G2190922o9, Consistently_detected == "true",
                          select=c("ASV","Sequence","PCR.1","PCR.2","PCR.3"))
G2190922o9_true$Fish <- "G2190922"
G2190922o9_true$Output <- "output9"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
G2190922o9sum <- as.data.frame(colSums(G2190922o9[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
G2190922o9sum <- as.data.frame(t(G2190922o9sum))

#add column of avg depth across all PCR reps
G2190922o9sum$Avg.depth <-apply(G2190922o9sum,1,mean)
#add columns of total depth
G2190922o9sum$Total.depth <-apply(G2190922o9sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
G2190922o9sum$Num.ASVs <- sum(G2190922o9$PCR.1>0 | G2190922o9$PCR.2>0 | G2190922o9$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
G2190922o9sum$Output <- "output9"

#add column of if in all 3 reps
G2190922o9sum$All.reps <- "No"

#only true ASVs
G2190922o9_truesum <- as.data.frame(colSums(G2190922o9_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
G2190922o9_truesum <- as.data.frame(t(G2190922o9_truesum))

#add column of avg depth across all PCR reps
G2190922o9_truesum $Avg.depth <-apply(G2190922o9_truesum,1,mean)

#add columns of total depth
G2190922o9_truesum $Total.depth <-apply(G2190922o9_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
G2190922o9_truesum$Num.ASVs <- sum(G2190922o9_true$PCR.1>0 | G2190922o9_true$PCR.2>0 | G2190922o9_true$PCR.3>0)

#add column of pipeline (so can plot with different outputs for comparison)
G2190922o9_truesum$Output <- "output9"

#add column of if in all 3 reps
G2190922o9_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
G2190922_output9 <- rbind(G2190922o9_truesum, G2190922o9sum)

#add column of fish ID
G2190922_output9$Fish <- "190922-G2"

#190922-H2----
H2190922o9 <- output9[c(1:2,45:47)]
names(H2190922o9)[3] <- "PCR.1"
names(H2190922o9)[4] <- "PCR.2"
names(H2190922o9)[5] <- "PCR.3"

#determine if the seq is detected in all 3 reps
H2190922o9$Consistently_detected <- ifelse(H2190922o9$PCR.1==0 | H2190922o9$PCR.2== 0 | H2190922o9$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
H2190922o9_true <- subset(H2190922o9, Consistently_detected == "true",
                          select=c("ASV","Sequence","PCR.1","PCR.2","PCR.3"))
H2190922o9_true$Fish <- "H2190922"
H2190922o9_true$Output <- "output9"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
H2190922o9sum <- as.data.frame(colSums(H2190922o9[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
H2190922o9sum <- as.data.frame(t(H2190922o9sum))

#add column of avg depth across all PCR reps
H2190922o9sum$Avg.depth <-apply(H2190922o9sum,1,mean)
#add columns of total depth
H2190922o9sum$Total.depth <-apply(H2190922o9sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
H2190922o9sum$Num.ASVs <- sum(H2190922o9$PCR.1>0 | H2190922o9$PCR.2>0 | H2190922o9$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
H2190922o9sum$Output <- "output9"

#add column of if in all 3 reps
H2190922o9sum$All.reps <- "No"

#only true ASVs
H2190922o9_truesum <- as.data.frame(colSums(H2190922o9_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
H2190922o9_truesum <- as.data.frame(t(H2190922o9_truesum))

#add column of avg depth across all PCR reps
H2190922o9_truesum $Avg.depth <-apply(H2190922o9_truesum,1,mean)

#add columns of total depth
H2190922o9_truesum $Total.depth <-apply(H2190922o9_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
H2190922o9_truesum$Num.ASVs <- sum(H2190922o9_true$PCR.1>0 | H2190922o9_true$PCR.2>0 | H2190922o9_true$PCR.3>0)

#add column of pipeline (so can plot with different outputs for comparison)
H2190922o9_truesum$Output <- "output9"

#add column of if in all 3 reps
H2190922o9_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
H2190922_output9 <- rbind(H2190922o9_truesum, H2190922o9sum)

#add column of fish ID
H2190922_output9$Fish <- "190922-H2"

#190922-J2----
J2190922o9 <- output9[c(1:2,48:50)]
names(J2190922o9)[3] <- "PCR.1"
names(J2190922o9)[4] <- "PCR.2"
names(J2190922o9)[5] <- "PCR.3"

#determine if the seq is detected in all 3 reps
J2190922o9$Consistently_detected <- ifelse(J2190922o9$PCR.1==0 | J2190922o9$PCR.2== 0 | J2190922o9$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
J2190922o9_true <- subset(J2190922o9, Consistently_detected == "true",
                          select=c("ASV","Sequence","PCR.1","PCR.2","PCR.3"))
J2190922o9_true$Fish <- "J2190922"
J2190922o9_true$Output <- "output9"


#compare size of two dfs (# ASVs and total reads)
#all ASVs
J2190922o9sum <- as.data.frame(colSums(J2190922o9[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
J2190922o9sum <- as.data.frame(t(J2190922o9sum))

#add column of avg depth across all PCR reps
J2190922o9sum$Avg.depth <-apply(J2190922o9sum,1,mean)
#add columns of total depth
J2190922o9sum$Total.depth <-apply(J2190922o9sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
J2190922o9sum$Num.ASVs <- sum(J2190922o9$PCR.1>0 | J2190922o9$PCR.2>0 | J2190922o9$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
J2190922o9sum$Output <- "output9"

#add column of if in all 3 reps
J2190922o9sum$All.reps <- "No"

#only true ASVs
J2190922o9_truesum <- as.data.frame(colSums(J2190922o9_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
J2190922o9_truesum <- as.data.frame(t(J2190922o9_truesum))

#add column of avg depth across all PCR reps
J2190922o9_truesum $Avg.depth <-apply(J2190922o9_truesum,1,mean)

#add columns of total depth
J2190922o9_truesum $Total.depth <-apply(J2190922o9_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
J2190922o9_truesum$Num.ASVs <- sum(J2190922o9_true$PCR.1>0 | J2190922o9_true$PCR.2>0 | J2190922o9_true$PCR.3>0)

#add column of pipeline (so can plot with different outputs for comparison)
J2190922o9_truesum$Output <- "output9"

#add column of if in all 3 reps
J2190922o9_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
J2190922_output9 <- rbind(J2190922o9_truesum, J2190922o9sum)

#add column of fish ID
J2190922_output9$Fish <- "190922-J2"






#put it all all together ----
all_output9 <- rbind(A1190222_output9,A2190222_output9,H1190222_output9,M1190222_output9,
                     A1190918_output9,A3190918_output9,A4190918_output9,C4190918_output9,
                     D4190918_output9,E7190918_output9,F2190918_output9,D2190922_output9,
                     F2190922_output9,G2190922_output9,H2190922_output9,J2190922_output9)

#reorganizae
rownames(all_output9)<-c(1:32)
library(tidyverse)
all_output9 <- all_output9 %>% relocate(Fish, .before = PCR.1)

#save for future analysis
write_csv(all_output9, "dada2-optimization/organized_ASVs/output9_ASVanalysis.csv")


#change data types
all_output9$Fish <- as.factor(all_output9$Fish)
all_output9$Output <- as.factor(all_output9$Output)
all_output9$All.reps <- as.factor(all_output9$All.reps)


#averages----
#want to add a label to each bar showing how many ASVs this average was taken from
avgASV <- aggregate(Num.ASVs ~ Fish + All.reps, data = all_output9, FUN = mean)
avgdepth <- aggregate(Avg.depth ~ Fish + All.reps, data = all_output9, FUN = mean)

#plot with bars for individual fish
ggplot(avgdepth, aes(fill=All.reps, y=Avg.depth, x=Fish)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Pipeline") + ylab("Average read depth") + 
  labs(fill= "ASV detected in all PCR replicates?") +
  geom_text(aes(label = round(avgASV$Num.ASVs, digits = 0)), position=position_dodge(width=0.9), vjust=-0.25, size=3)


#average across ALL fish
avgallASV <- aggregate(Num.ASVs ~ All.reps, data = all_output9, FUN = mean)
avgalldepth <- aggregate(Avg.depth ~ All.reps, data = all_output9, FUN = mean)

ggplot(avgalldepth, aes(fill=All.reps, y=Avg.depth, x=All.reps)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("ASV detected in all PCR replicates?") + ylab("Average read depth") + 
  labs(fill= "ASV detected in all PCR replicates?") +
  ggtitle("ASV analysis of output9") +
  guides(fill="none") +
  geom_text(aes(label = round(avgallASV$Num.ASVs, digits = 0)), position=position_dodge(width=0.9), vjust=-0.25, size=3)


#total depth-----
#want to add a label to each bar showing how many ASVs this average was taken from
numASV <- aggregate(Num.ASVs ~ All.reps, data = all_output9, FUN = sum)
totaldepth <- aggregate(Total.depth ~ All.reps, data = all_output9, FUN = sum)

#plot with bars for individual fish
ggplot(totaldepth, aes(fill=All.reps, y=Total.depth, x=All.reps)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none") +
  ggtitle("Total depth and ASV analysis of output9") +
  xlab("ASV detected in all PCR replicates?") + ylab("Total read depth") + 
  labs(fill= "ASV detected in all PCR replicates?") +
  geom_text(aes(label = round(numASV$Num.ASVs, digits = 0)), position=position_dodge(width=0.9), vjust=-0.25, size=3)

#the number above the bar is the TOTAL number of ASVs the depth came from,
# NOT the number of unique ASVs

#only looking at consistently detected ASVs----
alltrueo9 <- rbind(A1190222o9_true,A2190222o9_true,H1190222o9_true,M1190222o9_true,
                   A1190918o9_true,A3190918o9_true,A4190918o9_true,C4190918o9_true,
                   D4190918o9_true,E7190918o9_true,F2190918o9_true,D2190922o9_true,
                   F2190922o9_true,G2190922o9_true,H2190922o9_true,J2190922o9_true)

#save for future analysis
write_csv(alltrueo9, "dada2-optimization/organized_ASVs/output9_TrueASVs.csv")

uniqueASVo9 <- unique(alltrueo9$ASV)
length(unique(alltrueo9$ASV))

#1: 14 ASVs seen in all FISH (but not all of these ASVs are present in EACH fish)
#2: 9 ASVs
#3: "
#4: "
#5: "
#6: 28
#7: "
#8: 9
#9: 28


#Total depth of unique ASVs retained across whole pipeline

#add column of total depth across all PCR reps
alltrueo9$Total.depth <-apply(alltrueo9[,3:5],1,sum)
sum(alltrueo9$Total.depth)

#total depth of unique ASVs summed across all individuals
#1: 1,165,763
#2: 1,158,009
#3: 1,142,352
#4:1161563
#5:1162159
#6:1158141
#7:1158141
#8:157007
#9:169241


#if want to see TOTAL depth per fish:
aggregate(Total.depth ~ Fish, data = alltrueo9, FUN = sum)


#################################3

#this is not that informative, but interesting to see how depth per fish is
#total depth retained = most informative metric
#avg depth of unique ASVs per fish (avg. across PCR reps) summed across all individuals
#aka total depth of the average depth per fish

summary(alltrueo9)
alltrueo9$ASV <- as.factor(alltrueo9$ASV)
alltrueo9$Fish <- as.factor(alltrueo9$Fish)
alltrueo9$Output <- as.factor(alltrueo9$Output)

#add column of avg depth across all PCR reps
alltrueo9$Avg.depth <-apply(alltrueo9[,3:5],1,mean)
sum(alltrueo9$Avg.depth)

#1: 388,587.7
#2: 386,003
#3: 380,784
#4:387187.7
#5:387386.3
#6:386047
#7:386047
#8:52335.67
#9:56413.67


#if want to see AVG depth per fish:
aggregate(Avg.depth ~ Fish, data = alltrueo9, FUN = sum)







