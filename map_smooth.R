#!/usr/bin/env Rscript 


######## make function MapToRates ########

mapToRates <- function(chrMapFile, PLOT = F){
  chr_length <- dim(chrMapFile)[1]
  rates <- numeric(length = chr_length)

  for(i in 1:chr_length-1){
    rates[i] <- (chrMapFile$H[i+1]-chrMapFile$H[i])/(chrMapFile$SNP[i+1]-chrMapFile$SNP[i])*1e6 # multiply by 1e6 to convert bp to Mb
  }
  rates[chr_length] <- 0 # the last rate should be zero since there is no next marker - this was the case in the drosophila and human maps

  if(PLOT==T){
    plot(chrMapFile$SNP, rates,type="l", ylim = c(0,20))
  }

  return(rates)
}


setwd("~/popSim/arabidopsis_recomb_maps/arab_chr1/")
arab_chr1 <- list.files(pattern = "chr1")
## initialize vectors with first population then loop over others
P2_chr1 <- read.table("P2_chr1.txt", header = T)
rates_P2_chr1 <- mapToRates(P2_chr1)
rates_all_chr1 <- mapToRates(P2_chr1)
snp_all_chr1 <- P2_chr1$SNP
for(i in 2:length(arab_chr1)){
  input <- read.table(arab_chr1[i], header = T)
  rates <- mapToRates(input)
  rates_all_chr1 <- c(rates_all_chr1, mapToRates(input))
  snp_all_chr1 <- c(snp_all_chr1, input$SNP)

}

## make data frame of all snps and rates and order by snps
rates_snps_chr1 <- data.frame(snp_all_chr1,rates_all_chr1)
rates_snps_chr1_sorted <- rates_snps_chr1[order(rates_snps_chr1$snp_all_chr1),]

## generate loess curve -- span parameter is important for determining smoothing, 0.09 seems reasonable
loess_chr1 <- loess(rates_snps_chr1_sorted$rates_all_chr1 ~ rates_snps_chr1_sorted$snp_all_chr1, data=rates_snps_chr1_sorted, span=0.09)
coords = seq(50000,tail(sort(snp_all_chr1),1),50000)
smoothed_chr1 <- round(predict(loess_chr1,newdata=coords),6)
chr1_loess <- data.frame(coords,smoothed_chr1)
colnames(chr1_loess) <- c("Position(bp)", "Rate(cM/Mb)")
chr1_loess <- cbind(Chromosome = "chr1",chr1_loess, "Map(cM)" = "0")
chr1_loess$`Rate(cM/Mb)`[dim(chr1_loess)[1]] <- 0 # this is by definition 0 yet loess will give it a nonzero value (sometimes negative!!)
write.table(x = chr1_loess, file = "~/popSim/arabidopsis_recomb_maps/arab_chr1_map_loess.txt", row.names = F, quote = F)

###### CHROMOSOME 2 #####
setwd("~/popSim/arabidopsis_recomb_maps/arab_chr2/")
arab_chr2 <- list.files(pattern = "chr2")
## initialize vectors with first population then loop over others
P2_chr2 <- read.table("P2_chr2.txt", header = T)
rates_P2_chr2 <- mapToRates(P2_chr2)
rates_all_chr2 <- mapToRates(P2_chr2)
snp_all_chr2 <- P2_chr2$SNP
for(i in 2:length(arab_chr2)){
  input <- read.table(arab_chr2[i], header = T)
  rates <- mapToRates(input)
  rates_all_chr2 <- c(rates_all_chr2, mapToRates(input))
  snp_all_chr2 <- c(snp_all_chr2, input$SNP)
}
## make data frame of all snps and rates and order by snps
rates_snps_chr2 <- data.frame(snp_all_chr2,rates_all_chr2)
rates_snps_chr2_sorted <- rates_snps_chr2[order(rates_snps_chr2$snp_all_chr2),]
## several rates are VERY high after SNP 6809099-- remove anything greater than 50??
max(rates_snps_chr2_sorted$rates_all_chr2)
rates_snps_chr2_sorted <- subset(rates_snps_chr2_sorted, rates_snps_chr2_sorted$rates_all_chr2 < 50)
## generate loess curve -- span parameter is important for determining smoothing, 0.09 seems reasonable
loess_chr2 <- loess(rates_snps_chr2_sorted$rates_all_chr2 ~ rates_snps_chr2_sorted$snp_all_chr2, data=rates_snps_chr2_sorted, span=0.09)
coords = seq(50000,tail(sort(snp_all_chr2),1),50000)
smoothed_chr2 <- round(predict(loess_chr2,newdata=coords),6)
## make data frame of snps and loess estimates and remove redundant lines to get final map
chr2_loess <- data.frame(coords,smoothed_chr2)
colnames(chr2_loess) <- c("Position(bp)", "Rate(cM/Mb)")
chr2_loess <- cbind(Chromosome = "chr2",chr2_loess, "Map(cM)" = "0")
chr2_loess$`Rate(cM/Mb)`[dim(chr2_loess)[1]] <- 0
write.table(x = chr2_loess, file = "~/popSim/arabidopsis_recomb_maps/arab_chr2_map_loess.txt", row.names = F, quote = F)

###### CHROMOSOME 3 #####
setwd("~/popSim/arabidopsis_recomb_maps/arab_chr3/")
arab_chr3 <- list.files(pattern = "chr3")
## initialize vectors with first population then loop over others
P2_chr3 <- read.table("P2_chr3.txt", header = T)
rates_P2_chr3 <- mapToRates(P2_chr3)
rates_all_chr3 <- mapToRates(P2_chr3)
snp_all_chr3 <- P2_chr3$SNP
for(i in 2:length(arab_chr3)){
  input <- read.table(arab_chr3[i], header = T)
  rates <- mapToRates(input)
  rates_all_chr3 <- c(rates_all_chr3, mapToRates(input))
  snp_all_chr3 <- c(snp_all_chr3, input$SNP)
}
## make data frame of all snps and rates and order by snps
rates_snps_chr3 <- data.frame(snp_all_chr3,rates_all_chr3)
rates_snps_chr3_sorted <- rates_snps_chr3[order(rates_snps_chr3$snp_all_chr3),]
## several rates are VERY high after SNP 6809099-- remove anything greater than 50??
max(rates_snps_chr3_sorted$rates_all_chr3)
rates_snps_chr3_sorted <- subset(rates_snps_chr3_sorted, rates_snps_chr3_sorted$rates_all_chr3 < 50)
## generate loess curve -- span parameter is important for determining smoothing, 0.09 seems reasonable
loess_chr3 <- loess(rates_snps_chr3_sorted$rates_all_chr3 ~ rates_snps_chr3_sorted$snp_all_chr3, data=rates_snps_chr3_sorted, span=0.09)
coords = seq(50000,tail(sort(snp_all_chr3),1),50000)
smoothed_chr3 <- round(predict(loess_chr3,newdata=coords),6)
## make data frame of snps and loess estimates and remove redundant lines to get final map
chr3_loess <- data.frame(coords,smoothed_chr3)
colnames(chr3_loess) <- c("Position(bp)", "Rate(cM/Mb)")
chr3_loess <- cbind(Chromosome = "chr3",chr3_loess, "Map(cM)" = "0")
chr3_loess$`Rate(cM/Mb)`[dim(chr3_loess)[1]] <- 0
write.table(x = chr3_loess, file = "~/popSim/arabidopsis_recomb_maps/arab_chr3_map_loess.txt", row.names = F, quote = F)

###### CHROMOSOME 4 #####
setwd("~/popSim/arabidopsis_recomb_maps/arab_chr4/")
arab_chr4 <- list.files(pattern = "chr4")
## initialize vectors with first population then loop over others
P2_chr4 <- read.table("P2_chr4.txt", header = T)
rates_P2_chr4 <- mapToRates(P2_chr4)
rates_all_chr4 <- mapToRates(P2_chr4)
snp_all_chr4 <- P2_chr4$SNP
for(i in 2:length(arab_chr4)){
  input <- read.table(arab_chr4[i], header = T)
  rates <- mapToRates(input)
  rates_all_chr4 <- c(rates_all_chr4, mapToRates(input))
  snp_all_chr4 <- c(snp_all_chr4, input$SNP)
}
## make data frame of all snps and rates and order by snps
rates_snps_chr4 <- data.frame(snp_all_chr4,rates_all_chr4)
rates_snps_chr4_sorted <- rates_snps_chr4[order(rates_snps_chr4$snp_all_chr4),]
## several rates are VERY high after SNP 6809099-- remove anything greater than 50??
max(rates_snps_chr4_sorted$rates_all_chr4)
rates_snps_chr4_sorted <- subset(rates_snps_chr4_sorted, rates_snps_chr4_sorted$rates_all_chr4 < 50)
## generate loess curve -- span parameter is important for determining smoothing, 0.09 seems reasonable
loess_chr4 <- loess(rates_snps_chr4_sorted$rates_all_chr4 ~ rates_snps_chr4_sorted$snp_all_chr4, data=rates_snps_chr4_sorted, span=0.09)
coords = seq(50000,tail(sort(snp_all_chr4),1),50000)
smoothed_chr4 <- round(predict(loess_chr4,newdata=coords),6)
## make data frame of snps and loess estimates and remove redundant lines to get final map
chr4_loess <- data.frame(coords,smoothed_chr4)
colnames(chr4_loess) <- c("Position(bp)", "Rate(cM/Mb)")
chr4_loess <- cbind(Chromosome = "chr4",chr4_loess, "Map(cM)" = "0")
chr4_loess$`Rate(cM/Mb)`[dim(chr4_loess)[1]] <- 0
write.table(x = chr4_loess, file = "~/popSim/arabidopsis_recomb_maps/arab_chr4_map_loess.txt", row.names = F, quote = F)

###### CHROMOSOME 5 #####
setwd("~/popSim/arabidopsis_recomb_maps/arab_chr5/")
arab_chr5 <- list.files(pattern = "chr5")
## initialize vectors with first population then loop over others
P2_chr5 <- read.table("P2_chr5.txt", header = T)
rates_P2_chr5 <- mapToRates(P2_chr5)
rates_all_chr5 <- mapToRates(P2_chr5)
snp_all_chr5 <- P2_chr5$SNP
for(i in 2:length(arab_chr5)){
  input <- read.table(arab_chr5[i], header = T)
  rates <- mapToRates(input)
  rates_all_chr5 <- c(rates_all_chr5, mapToRates(input))
  snp_all_chr5 <- c(snp_all_chr5, input$SNP)
}
## make data frame of all snps and rates and order by snps
rates_snps_chr5 <- data.frame(snp_all_chr5,rates_all_chr5)
rates_snps_chr5_sorted <- rates_snps_chr5[order(rates_snps_chr5$snp_all_chr5),]
## several rates are VERY high after SNP 6809099-- remove anything greater than 50??
max(rates_snps_chr5_sorted$rates_all_chr5)
rates_snps_chr5_sorted <- subset(rates_snps_chr5_sorted, rates_snps_chr5_sorted$rates_all_chr5 < 50)
## generate loess curve -- span parameter is important for determining smoothing, 0.09 seems reasonable
loess_chr5 <- loess(rates_snps_chr5_sorted$rates_all_chr5 ~ rates_snps_chr5_sorted$snp_all_chr5, data=rates_snps_chr5_sorted, span=0.09)
coords = seq(50000,tail(sort(snp_all_chr5),1),50000)
smoothed_chr5 <- round(predict(loess_chr5,newdata=coords),6)
## make data frame of snps and loess estimates and remove redundant lines to get final map
chr5_loess <- data.frame(coords,smoothed_chr5)
colnames(chr5_loess) <- c("Position(bp)", "Rate(cM/Mb)")
chr5_loess <- cbind(Chromosome = "chr5",chr5_loess, "Map(cM)" = "0")
chr5_loess$`Rate(cM/Mb)`[dim(chr5_loess)[1]] <- 0
write.table(x = chr5_loess, file = "~/popSim/arabidopsis_recomb_maps/arab_chr5_map_loess.txt", row.names = F, quote = F)


