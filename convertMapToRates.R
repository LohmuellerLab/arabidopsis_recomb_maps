
## These scripts are for converting the Arabidopsis thaliana recombination rates from Salome et al. 2012 to genetic maps
## the file from Salome 2012 that was used was their "maps K H CF.xls" found on dryad https://datadryad.org/resource/doi:10.5061/dryad.v655ns36
## this file includes map distances from crosses of 20 A. thaliana lines
## these scripts convert the map distances to rates and then average across the 20 crosses using loess

setwd("~/Box/UCLA/PopSim/Arabidopsis_recomb_map")


###### Example using human map that already has rates and map distances as sanity check ########
human_map <- read.table("HapmapII_GRCh37_RecombinationHotspots/genetic_map_GRCh37_chr1.txt", header=T)
head(human_map)

human_length <-  dim(human_map)[1]

hum_rates <- numeric(length = human_length)

for(i in 1:human_length-1){
  hum_rates[i] <- (human_map$Map.cM.[i+1]-human_map$Map.cM.[i])/(human_map$Position.bp.[i+1]-human_map$Position.bp.[i])*1e6
}
hum_rates[human_length] <- 0


human_map$rates_calc <- hum_rates

head(hum_rates)
tail(hum_rates)

start <- human_length-1000
end <- human_length

plot(hum_rates[start:end], type ="l", lwd = 4)
lines(human_map$Rate.cM.Mb.[start:end], col=2, lty=2, lwd=2)


# note that their rates differ somewhat from the ones I compute because of rounding error - especially problematic when markers are 1 bp apart and genetic distance is very small
# this becomes a bit of a problem for arabidopsis for a few markers that are close together and lead to very high estimates of recombination rates (>1000 cM/Mb)
diff <- human_map$Rate.cM.Mb.-hum_rates
human_map$diff <- diff
plot(diff)




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






###### run on Arabidopsis files #########
## note that the input data for this is the "maps K H CF.xls" file from the SI from Salome et al. 2012: https://datadryad.org/resource/doi:10.5061/dryad.v655ns36
## the data for each cross/chromosome was copy/pasted into individual text files which can then be read in by chromosome
## there were a few peculiarities with the data that were corrected, for example some SNPs were repeated or out of order for two of the chr1 crosses, and the P7 cross was missing some data for chr 2


#### CHROMOSOME 1 #####
setwd("~/Box/UCLA/PopSim/Arabidopsis_recomb_map/arab_chr1/")
arab_chr1 <- list.files(pattern = "chr1")
arab_chr1

## plot each pop rates individually if wanted
#par(mfrow = c(3,3))
for(i in 1:length(arab_chr1)){
  input <- read.table(arab_chr1[i], header = T)
  rates <- mapToRates(input, PLOT=T)
}

## initialize vectors with first population then loop over others
P2_chr1 <- read.table("P2_chr1.txt", header = T)
rates_P2_chr1 <- mapToRates(P2_chr1)
rates_all_chr1 <- mapToRates(P2_chr1)
snp_all_chr1 <- P2_chr1$SNP

plot(P2_chr1$SNP, rates_P2_chr1, type = "l", ylim=c(0,15))
for(i in 2:length(arab_chr1)){
  input <- read.table(arab_chr1[i], header = T)
  rates <- mapToRates(input)
  lines(input$SNP,rates, col=i)
  rates_all_chr1 <- c(rates_all_chr1, mapToRates(input))
  snp_all_chr1 <- c(snp_all_chr1, input$SNP)

}

## make data frame of all snps and rates and order by snps
rates_snps_chr1 <- data.frame(snp_all_chr1,rates_all_chr1)
rates_snps_chr1_sorted <- rates_snps_chr1[order(rates_snps_chr1$snp_all_chr1),]

## generate loess curve -- span parameter is important for determining smoothing, 0.09 seems reasonable
loess_chr1 <- loess(rates_snps_chr1_sorted$rates_all_chr1 ~ rates_snps_chr1_sorted$snp_all_chr1, data=rates_snps_chr1_sorted, span=0.09)
smoothed_chr1 <- round(predict(loess_chr1),6)

## make data frame of snps and loess estimates and remove redundant lines to get final map
chr1_loess <- data.frame(rates_snps_chr1_sorted$snp_all_chr1,smoothed_chr1)
chr1_loess <- unique(chr1_loess)
colnames(chr1_loess) <- c("Position(bp)", "Rate(cM/Mb)")
chr1_loess <- cbind(Chromosome = "chr1",chr1_loess, "Map(cM)" = "0")
head(chr1_loess)
chr1_loess$`Rate(cM/Mb)`[dim(chr1_loess)[1]] <- 0 # this is by definition 0 yet loess will give it a nonzero value (sometimes negative!!)

plot(rates_snps_chr1_sorted$snp_all_chr1, rates_snps_chr1_sorted$rates_all_chr1, pch=19, col = "gray", xlab="Position (Mb)", ylab = "Rate (cM/Mb)", main = "Chromosome 1")
lines(smoothed_chr1, x = rates_snps_chr1_sorted$snp_all_chr1, col = "blue", lwd=3)
abline(a=mean(chr1_loess$`Rate(cM/Mb)`),b=0, lwd=3, lty=2)

#write file with recomb map and save plot
setwd("~/Box/UCLA/PopSim/Arabidopsis_recomb_map")
write.table(x = chr1_loess, file = "arab_chr1_map_loess.txt", row.names = F, quote = F)

png("arab_chr1_plot_loess.png", width=8, height=6, units = "in", res=600)
plot(rates_snps_chr1_sorted$snp_all_chr1, rates_snps_chr1_sorted$rates_all_chr1, pch=19, col = "gray", xlab="Position (Mb)", ylab = "Rate (cM/Mb)", main = "Chromosome 1")
lines(smoothed_chr1, x = rates_snps_chr1_sorted$snp_all_chr1, col = "blue", lwd=3)
abline(a=mean(chr1_loess$`Rate(cM/Mb)`),b=0, lwd=3, lty=2)
dev.off()




###### CHROMOSOME 2 #####
setwd("~/Box/UCLA/PopSim/Arabidopsis_recomb_map/arab_chr2/")
arab_chr2 <- list.files(pattern = "chr2")
arab_chr2

## plot each pop rates individually - good for checking for errors
#par(mfrow = c(3,3))
for(i in 1:length(arab_chr2)){
  input <- read.table(arab_chr2[i], header = T)
  rates <- mapToRates(input, PLOT=T)
}

## initialize vectors with first population then loop over others
P2_chr2 <- read.table("P2_chr2.txt", header = T)
rates_P2_chr2 <- mapToRates(P2_chr2)
rates_all_chr2 <- mapToRates(P2_chr2)
snp_all_chr2 <- P2_chr2$SNP

plot(P2_chr2$SNP, rates_P2_chr2, type = "l", ylim=c(0,35))
for(i in 2:length(arab_chr2)){
  input <- read.table(arab_chr2[i], header = T)
  rates <- mapToRates(input)
  lines(input$SNP,rates, col=i)
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
smoothed_chr2 <- round(predict(loess_chr2),6)

## make data frame of snps and loess estimates and remove redundant lines to get final map
chr2_loess <- data.frame(rates_snps_chr2_sorted$snp_all_chr2,smoothed_chr2)
chr2_loess <- unique(chr2_loess)
colnames(chr2_loess) <- c("Position(bp)", "Rate(cM/Mb)")
chr2_loess <- cbind(Chromosome = "chr2",chr2_loess, "Map(cM)" = "0")
head(chr2_loess)
chr2_loess$`Rate(cM/Mb)`[dim(chr2_loess)[1]] <- 0 # this is by definition 0 yet loess will give it a nonzero value (sometimes negative!!)

plot(rates_snps_chr2_sorted$snp_all_chr2, rates_snps_chr2_sorted$rates_all_chr2, pch=19, col = "gray", xlab="Position (Mb)", ylab = "Rate (cM/Mb)", main = "Chromosome 2")
lines(smoothed_chr2, x = rates_snps_chr2_sorted$snp_all_chr2, col = "blue", lwd=3)
abline(a=mean(chr2_loess$`Rate(cM/Mb)`),b=0, lwd=3, lty=2)

#write file with recomb map and save plot
setwd("~/Box/UCLA/PopSim/Arabidopsis_recomb_map")
write.table(x = chr2_loess, file = "arab_chr2_map_loess.txt", row.names = F, quote = F)

png("arab_chr2_plot_loess.png", width=8, height=6, units = "in", res=600)
plot(rates_snps_chr2_sorted$snp_all_chr2, rates_snps_chr2_sorted$rates_all_chr2, pch=19, col = "gray", xlab="Position (Mb)", ylab = "Rate (cM/Mb)", main = "Chromosome 2")
lines(smoothed_chr2, x = rates_snps_chr2_sorted$snp_all_chr2, col = "blue", lwd=3)
abline(a=mean(chr2_loess$`Rate(cM/Mb)`),b=0, lwd=3, lty=2)
dev.off()







###### CHROMOSOME 3 #####
setwd("~/Box/UCLA/PopSim/Arabidopsis_recomb_map/arab_chr3/")
arab_chr3 <- list.files(pattern = "chr3")
arab_chr3

## plot each pop rates individually - good for checking for errors
#par(mfrow = c(3,3))
for(i in 1:length(arab_chr3)){
  input <- read.table(arab_chr3[i], header = T)
  rates <- mapToRates(input, PLOT=T)
}

## initialize vectors with first population then loop over others
P2_chr3 <- read.table("P2_chr3.txt", header = T)
rates_P2_chr3 <- mapToRates(P2_chr3)
rates_all_chr3 <- mapToRates(P2_chr3)
snp_all_chr3 <- P2_chr3$SNP

plot(P2_chr3$SNP, rates_P2_chr3, type = "l", ylim=c(0,35))
for(i in 2:length(arab_chr3)){
  input <- read.table(arab_chr3[i], header = T)
  rates <- mapToRates(input)
  lines(input$SNP,rates, col=i)
  rates_all_chr3 <- c(rates_all_chr3, mapToRates(input))
  snp_all_chr3 <- c(snp_all_chr3, input$SNP)
  
}

## make data frame of all snps and rates and order by snps
rates_snps_chr3 <- data.frame(snp_all_chr3,rates_all_chr3)
rates_snps_chr3_sorted <- rates_snps_chr3[order(rates_snps_chr3$snp_all_chr3),]

## several rates are VERY high after SNP 13495379-- remove anything greater than 50??
max(rates_snps_chr3_sorted$rates_all_chr3) 
rates_snps_chr3_sorted <- subset(rates_snps_chr3_sorted, rates_snps_chr3_sorted$rates_all_chr3 < 50)

## generate loess curve -- span parameter is important for determining smoothing, 0.09 seems reasonable
loess_chr3 <- loess(rates_snps_chr3_sorted$rates_all_chr3 ~ rates_snps_chr3_sorted$snp_all_chr3, data=rates_snps_chr3_sorted, span=0.09)
smoothed_chr3 <- round(predict(loess_chr3),6)

## make data frame of snps and loess estimates and remove redundant lines to get final map
chr3_loess <- data.frame(rates_snps_chr3_sorted$snp_all_chr3,smoothed_chr3)
chr3_loess <- unique(chr3_loess)
colnames(chr3_loess) <- c("Position(bp)", "Rate(cM/Mb)")
chr3_loess <- cbind(Chromosome = "chr3",chr3_loess, "Map(cM)" = "0")
head(chr3_loess)
chr3_loess$`Rate(cM/Mb)`[dim(chr3_loess)[1]] <- 0 # this is by definition 0 yet loess will give it a nonzero value (sometimes negative!!)

plot(rates_snps_chr3_sorted$snp_all_chr3, rates_snps_chr3_sorted$rates_all_chr3, pch=19, col = "gray", xlab="Position (Mb)", ylab = "Rate (cM/Mb)", main = "Chromosome 3")
lines(smoothed_chr3, x = rates_snps_chr3_sorted$snp_all_chr3, col = "blue", lwd=3)
abline(a=mean(chr3_loess$`Rate(cM/Mb)`),b=0, lwd=3, lty=2)

#write file with recomb map and save plot
setwd("~/Box/UCLA/PopSim/Arabidopsis_recomb_map")
write.table(x = chr3_loess, file = "arab_chr3_map_loess.txt", row.names = F, quote = F)

png("arab_chr3_plot_loess.png", width=8, height=6, units = "in", res=600)
plot(rates_snps_chr3_sorted$snp_all_chr3, rates_snps_chr3_sorted$rates_all_chr3, pch=19, col = "gray", xlab="Position (Mb)", ylab = "Rate (cM/Mb)", main = "Chromosome 3")
lines(smoothed_chr3, x = rates_snps_chr3_sorted$snp_all_chr3, col = "blue", lwd=3)
abline(a=mean(chr3_loess$`Rate(cM/Mb)`),b=0, lwd=3, lty=2)
dev.off()






###### CHROMOSOME 4 #####
setwd("~/Box/UCLA/PopSim/Arabidopsis_recomb_map/arab_chr4/")
arab_chr4 <- list.files(pattern = "chr4")
arab_chr4

## plot each pop rates individually - good for checking for errors
#par(mfrow = c(3,3))
for(i in 1:length(arab_chr4)){
  input <- read.table(arab_chr4[i], header = T)
  rates <- mapToRates(input, PLOT=T)
}

## initialize vectors with first population then loop over others
P2_chr4 <- read.table("P2_chr4.txt", header = T)
rates_P2_chr4 <- mapToRates(P2_chr4)
rates_all_chr4 <- mapToRates(P2_chr4)
snp_all_chr4 <- P2_chr4$SNP

plot(P2_chr4$SNP, rates_P2_chr4, type = "l", ylim=c(0,35))
for(i in 2:length(arab_chr4)){
  input <- read.table(arab_chr4[i], header = T)
  rates <- mapToRates(input)
  lines(input$SNP,rates, col=i)
  rates_all_chr4 <- c(rates_all_chr4, mapToRates(input))
  snp_all_chr4 <- c(snp_all_chr4, input$SNP)
  
}

## make data frame of all snps and rates and order by snps
rates_snps_chr4 <- data.frame(snp_all_chr4,rates_all_chr4)
rates_snps_chr4_sorted <- rates_snps_chr4[order(rates_snps_chr4$snp_all_chr4),]

## no snps over ~35 for chr4 but leaving this in anyway
max(rates_snps_chr4_sorted$rates_all_chr4) 
rates_snps_chr4_sorted <- subset(rates_snps_chr4_sorted, rates_snps_chr4_sorted$rates_all_chr4 < 50)

## generate loess curve -- span parameter is important for determining smoothing, 0.09 seems reasonable
loess_chr4 <- loess(rates_snps_chr4_sorted$rates_all_chr4 ~ rates_snps_chr4_sorted$snp_all_chr4, data=rates_snps_chr4_sorted, span=0.09)
smoothed_chr4 <- round(predict(loess_chr4),6)

## make data frame of snps and loess estimates and remove redundant lines to get final map
chr4_loess <- data.frame(rates_snps_chr4_sorted$snp_all_chr4,smoothed_chr4)
chr4_loess <- unique(chr4_loess)
colnames(chr4_loess) <- c("Position(bp)", "Rate(cM/Mb)")
chr4_loess <- cbind(Chromosome = "chr4",chr4_loess, "Map(cM)" = "0")
head(chr4_loess)
chr4_loess$`Rate(cM/Mb)`[dim(chr4_loess)[1]] <- 0 # this is by definition 0 yet loess will give it a nonzero value (sometimes negative!!)
chr4_loess$`Rate(cM/Mb)`[chr4_loess$`Rate(cM/Mb)` < 0] <- 0  ## two negative values were generated by loess near centromere - replace with 0

plot(rates_snps_chr4_sorted$snp_all_chr4, rates_snps_chr4_sorted$rates_all_chr4, pch=19, col = "gray", xlab="Position (Mb)", ylab = "Rate (cM/Mb)", main = "Chromosome 4")
lines(smoothed_chr4, x = rates_snps_chr4_sorted$snp_all_chr4, col = "blue", lwd=3)
abline(a=mean(chr4_loess$`Rate(cM/Mb)`),b=0, lwd=3, lty=2)

#write file with recomb map and save plot
setwd("~/Box/UCLA/PopSim/Arabidopsis_recomb_map")
write.table(x = chr4_loess, file = "arab_chr4_map_loess.txt", row.names = F, quote = F)

png("arab_chr4_plot_loess.png", width=8, height=6, units = "in", res=600)
plot(rates_snps_chr4_sorted$snp_all_chr4, rates_snps_chr4_sorted$rates_all_chr4, pch=19, col = "gray", xlab="Position (Mb)", ylab = "Rate (cM/Mb)", main = "Chromosome 4")
lines(smoothed_chr4, x = rates_snps_chr4_sorted$snp_all_chr4, col = "blue", lwd=3)
abline(a=mean(chr4_loess$`Rate(cM/Mb)`),b=0, lwd=3, lty=2)
dev.off()







###### CHROMOSOME 5 #####
setwd("~/Box/UCLA/PopSim/Arabidopsis_recomb_map/arab_chr5/")
arab_chr5 <- list.files(pattern = "chr5")
arab_chr5

## plot each pop rates individually - good for checking for errors
#par(mfrow = c(3,3))
for(i in 1:length(arab_chr5)){
  input <- read.table(arab_chr5[i], header = T)
  rates <- mapToRates(input, PLOT=T)
}

## initialize vectors with first population then loop over others
P2_chr5 <- read.table("P2_chr5.txt", header = T)
rates_P2_chr5 <- mapToRates(P2_chr5)
rates_all_chr5 <- mapToRates(P2_chr5)
snp_all_chr5 <- P2_chr5$SNP

plot(P2_chr5$SNP, rates_P2_chr5, type = "l", ylim=c(0,35))
for(i in 2:length(arab_chr5)){
  input <- read.table(arab_chr5[i], header = T)
  rates <- mapToRates(input)
  lines(input$SNP,rates, col=i)
  rates_all_chr5 <- c(rates_all_chr5, mapToRates(input))
  snp_all_chr5 <- c(snp_all_chr5, input$SNP)
  
}

## make data frame of all snps and rates and order by snps
rates_snps_chr5 <- data.frame(snp_all_chr5,rates_all_chr5)
rates_snps_chr5_sorted <- rates_snps_chr5[order(rates_snps_chr5$snp_all_chr5),]

## no snps over ~35 for chr5 but leaving this in anyway
max(rates_snps_chr5_sorted$rates_all_chr5) 
rates_snps_chr5_sorted <- subset(rates_snps_chr5_sorted, rates_snps_chr5_sorted$rates_all_chr5 < 50)

## generate loess curve -- span parameter is important for determining smoothing, 0.09 seems reasonable
loess_chr5 <- loess(rates_snps_chr5_sorted$rates_all_chr5 ~ rates_snps_chr5_sorted$snp_all_chr5, data=rates_snps_chr5_sorted, span=0.09)
smoothed_chr5 <- round(predict(loess_chr5),6)

## make data frame of snps and loess estimates and remove redundant lines to get final map
chr5_loess <- data.frame(rates_snps_chr5_sorted$snp_all_chr5,smoothed_chr5)
chr5_loess <- unique(chr5_loess)
colnames(chr5_loess) <- c("Position(bp)", "Rate(cM/Mb)")
chr5_loess <- cbind(Chromosome = "chr5",chr5_loess, "Map(cM)" = "0")
head(chr5_loess)
chr5_loess$`Rate(cM/Mb)`[dim(chr5_loess)[1]] <- 0 # this is by definition 0 yet loess will give it a nonzero value (sometimes negative!!)
chr5_loess$`Rate(cM/Mb)`[chr5_loess$`Rate(cM/Mb)` < 0] <- 0  ## two negative values were generated by loess near centromere - replace with 0

plot(rates_snps_chr5_sorted$snp_all_chr5, rates_snps_chr5_sorted$rates_all_chr5, pch=19, col = "gray", xlab="Position (Mb)", ylab = "Rate (cM/Mb)", main = "Chromosome 5")
lines(smoothed_chr5, x = rates_snps_chr5_sorted$snp_all_chr5, col = "blue", lwd=3)
abline(a=mean(chr5_loess$`Rate(cM/Mb)`),b=0, lwd=3, lty=2)

#write file with recomb map and save plot
setwd("~/Box/UCLA/PopSim/Arabidopsis_recomb_map")
write.table(x = chr5_loess, file = "arab_chr5_map_loess.txt", row.names = F, quote = F)

png("arab_chr5_plot_loess.png", width=8, height=6, units = "in", res=600)
plot(rates_snps_chr5_sorted$snp_all_chr5, rates_snps_chr5_sorted$rates_all_chr5, pch=19, col = "gray", xlab="Position (Mb)", ylab = "Rate (cM/Mb)", main = "Chromosome 5")
lines(smoothed_chr5, x = rates_snps_chr5_sorted$snp_all_chr5, col = "blue", lwd=3)
abline(a=mean(chr5_loess$`Rate(cM/Mb)`),b=0, lwd=3, lty=2)
dev.off()












