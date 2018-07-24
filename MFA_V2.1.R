# Written by Dr. A M Mahedi Hasan, Post-doctoral Research Associate, The Leach Lab, ICB, University of Edinburgh, UK.
# This script is under GPLv3 licensing criteria.
# This script is for normalising depth of coverage (A.K.A count data) data from exponentially growing culture against total depth of each dataset and corresponding wild type stationary phase data.

# As output, this script will generate the following data files (.txt) in the "normalisation_output" directory:
# - raw_data.txt
# - normalized_raw_data.txt
# - raw_stationary_data.txt
# - normalized_stationary_data.txt
# - average_normalized_stationary_data.txt
# - normalised_to_depth_&_stn.txt
# - 1000bp_fixed_window_data.txt
# - loess data with a defined span based on the 1000bp_fixed_window_data.txt
#	- along with that, in the terminal, only a few end lines of datapoints will be printed with corresponding data heading(s).

# First, the depth of coverage data from the exponentially growing culture of a strain is normalized by the total depth of each data set.
# The depth or coverage data is read from the files starting with "depth" and seperated by "_" like "depth_DL4184Aara.txt". Please follow this naming pattern. The last part should be indicating the experimental details (meta information) without any gap followed by ".txt"

# On the other hand, the depth files from the stationary phase culture should be started with "stn", like "stn_depth_DL4184Aara.txt"
# MOST IMPORTANTLY, if you are missing the depth of coverage data from stationary phase culture, please comment-out the from line 75 to line 120. 


# The reference genome file should given in .fasta format and the name of the file should start with "genome", like "genome_ref.fasta".







print("The calculation started at:")
print(Sys.time())

# creating an output directory called "normalisation_output"
dir.create("normalisation_output", showWarnings = TRUE, recursive = FALSE)


# loading/installing nececssary package(s)
if (!require(seqinr)) install.packages('seqinr')
library(seqinr)

# dealing with genome_length
genome_file <- list.files(pattern = "^genome")
genome_length <- length(read.fasta(genome_file)[[1]])



# For the logarithmic-phase data:
files <- list.files(pattern="^depth")

data <- data.frame(position=seq(from=1,to=genome_length))

# The files are being read and name of the strains and their experimental details is being extracted here.
for(i in 1:length(files)){
  temp <- read.table(files[i])
  temp[,1] <- NULL
  log_name_split <- strsplit(as.character(files[i]),"_")[[1]]
  colnames(temp) <- c("position",strsplit(log_name_split[length(log_name_split)],".txt")[[1]])
  data <- merge(data,temp,by="position",all.x=TRUE)
}

print("Raw log-phase data")
print(tail(data))
write.table(data,"./normalisation_output/log-phase_raw_data.txt")

for(j in 2:ncol(data)){
  data[,j] <- data[,j]/sum(data[,j],na.rm=T)
}

print("Log-phase data normalised to total depth")
print(tail(data))
write.table(data,"./normalisation_output/normalized_log-phase_data.txt")


# For the stationary-phase data:
stn_files <- list.files(pattern="^stn")
stn_data <- data.frame(position=seq(from=1,to=genome_length))
for(i in 1:length(stn_files)){
  temp <- read.table(stn_files[i])
  temp[,1] <- NULL
  stn_name_split <- strsplit(as.character(stn_files[i]),"_")[[1]]
  colnames(temp) <- c("position",strsplit(stn_name_split[length(stn_name_split)],".txt")[[1]])
  stn_data <- merge(stn_data,temp,by="position",all.x=TRUE)
}

print("Raw stationary-phase data")
print(tail(stn_data))
write.table(stn_data,"./normalisation_output/stationary-phase_raw_data.txt")

for(j in 2:ncol(stn_data)){
  stn_data[,j] <- stn_data[,j]/sum(stn_data[,j],na.rm=T)
}

print("Stationary-phase data normalised to total depth")
print(tail(stn_data))
write.table(stn_data,"./normalisation_output/normalized_stationary-phase_data.txt")

# Generating average of stationary data:
average_stn <- data.frame(matrix(0,ncol=1,nrow=dim(stn_data)[1]))
for(k in 2:ncol(stn_data)){
  average_stn <- average_stn + stn_data[,k]
}
average_stn <- average_stn/(dim(stn_data)[2]-1)
average_stn <- data.frame(cbind(stn_data[,1]),average_stn)
colnames(average_stn) <- c("position","average_stn")

print("Average stationary-phase data normalised to depth")
print(tail(average_stn))
write.table(average_stn,"./normalisation_output/average_normalized_stationary-phase_data.txt")

# Normalising the logarithmic-phase data with the averaged stationary_phase data:
data <- merge(data,average_stn,by="position",all.x=TRUE)
for(x in 2:(ncol(data)-1)){
  data[,x] <- data[,x]/data[,ncol(data)]
}
data[,ncol(data)] <- NULL

print("Log-phase data normalised to total depth & stationary-phase depth")
print(tail(data))
write.table(data,"./normalisation_output/log-phase_data_normalised_to_depth_&_stn.txt")


# Generating 1000 bp fixed window average data:
position <- data.frame(position=seq(from=1,to=genome_length))
data2 <- merge(position,data,by="position",all.x=TRUE)
FW_data <- data.frame(matrix(NA,ncol=dim(data2)[2],nrow=floor((dim(data2)[1])/1000)))
colnames(FW_data) <- colnames(data2)
for(j in seq(from=1,to=floor(length(data2[,1])/1000))){
  FW_data[j,1] <- (j*1000)-500
}
for(i in 2:dim(data2)[2]){
  for(j in seq(from=1,to=floor(length(data2[,i])/1000)))
  {
    FW_data[j,i] <- mean(data2[((j*1000)-(1000-1)):(j*1000),i])
  }
}
print("1Kb Fixed_window data")
print(tail(FW_data))
write.table(FW_data,"./normalisation_output/1Kb_fixed_window data.txt")


# calculating the loess data

loess_data <- FW_data
for(i in 2:ncol(loess_data)){
  loess_data[,i] <- NA
}
for(i in 2:ncol(loess_data)){
  y.loess <- loess(y~x, span=0.05, data.frame(x=FW_data[,1],y=FW_data[,i]))
  loess_data[,i] <- predict(y.loess, data.frame(x=FW_data[,1]))
}
write.table(loess_data, "./normalisation_output/loess_1Kb_fixed_winodw_data.txt")



print("The calculation ended at:")
print(Sys.time())
