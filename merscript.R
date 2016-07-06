options(warn=0)
setwd("C:/Users/wsaoillab/Desktop/Ola Data Analysis/R/Oil/Merichem") #Set working directory

#Function for converting NA values to 0 in sample sets 
na.zero <- function (x)
{
  x[is.na(x)] <- 0
  return(x)
}

#Function for repeating rows
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
#Function for repeating columns 
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

#Function thats converts Inf values to 0 in sample sets
inf.zero <- function(x)
{
  x[is.infinite(as.matrix(x))] <- 0
  return(x) 
}

#Function that converts string to associated variable name
var.name <- function (str.name)
{
  variable <- eval(parse(text = str.name))
  return(variable)
}

#Acquire and Declare sample sets Mer(1 to 9)
for(i in 1:9)
{
  #assign name based on iteration and upload data from working directory
  assign(paste("mer",i, sep = ""),read.csv(paste("NAs-50ppm-M05-0",i,"Nas-O2-Brife.csv", sep = ""), 
  header = T, na.strings = c("N/A", "NA", "")))
  #converting n/a values to 0 with na.zero function through assign function
  assign(paste("mer",i, sep = ""), na.zero(eval(as.symbol(paste("mer",i, sep = "")))))
} 

ars <- data.frame(mer1$Component.Name) #Generates the Area Ratio Subtract(ARS) Spreadsheet
names(ars)[1] <- "Sample_Name"

#Fills in ARS with Merichem sample data
for (i in 1:9) 
{
  mer <- var.name(paste("mer",i, sep = "")) 
  ars <- cbind(ars, mer$Area) #Adds Merichem data to ars dataframe in iterations
  
  colnames(ars)[i+1] <- paste("Mer",i, sep = "") #Labels columns with appropriate name
}

#Determining ratio of merichem sample value to internal standard on ARS 
for (i in 1:9)
{
  ars[2:440,i+1] <- ars[2:440, i+1]/ars[1,i+1]
}

#Determining blank substracted merichem sample value on ARS
for (i in 1:9)
{
  ars[2:435,i+1] <- ars[2:435,i+1] - area_ratio_sub[which(substr(area_ratio_sub[2:695,1], 10, 11) <= 40)+1,2] #Subtracting Std-Mean from sample set
  ars[436:442,i+1] <- ars[436:442,i+1] - area_ratio_sub[696:702,2] #Subtracting Std-Mean from Surrogates 
}

#Converting all negative values to Zero
for (i in 1:9)
{
  ars[which(ars[,i+1] < 0),i+1] <- 0
}  

msumdata <- colSums(ars[2:435,2:10]) #Collecting the sums of each merichem sample column

#Sum for each carbon number value in merichem sample data
for (i in 1:35)
{
  if (i < 5){
    zero = "0"
  } else {
    zero = NULL
  } #Modifying variable that allows integers to represented as "01" for values less than 10
 
  cnum <- which(substr(ars[2:435,1], 10, 11) == paste(zero,i+5, sep = "")) # Determines the rows of each carbon number from C6 to C40
  cnumsum <- colSums(ars[cnum, 2:10])
  
  if (i == 1) {
    mcsumdata <- cnumsum #transfers sample data
    cnumsum <- NULL #nullify sample data to avoid error in rbind 
  } #Applies for the first iteration only to allow rbind in loop
  
  mcsumdata <- rbind(mcsumdata, cnumsum)
}
rownames(mcsumdata) <- 6:40
mcsumdata <- data.frame(mcsumdata)

#Sum for each z-value in merichem sample data
for (i in seq(0,24,2))
{
  if (i < 9){
    zero = "0"
  } else {
    zero = NULL
  } #Modifying variable that allows integers to represented as "01" for values less than 10
  
  znum <- which(substr(ars[2:435,1], 6, 7) == paste(zero,i, sep = "")) # Determines the rows of each z-value
  znumsum <- colSums(ars[znum,2:10])
  
  if (i == 0) {
    mzsumdata <- znumsum #transfers sample data
    znumsum <- NULL #nullify sample data to avoid error in rbind 
  } #Applies for the first iteration only to allow rbind in loop
  
  mzsumdata <- rbind(mzsumdata, znumsum)
}
rownames(mzsumdata) <- seq(0,24,2)
mzsumdata <- data.frame(mzsumdata)

#O2-Z-00-C06 - O2-Z-24-C60 divided by its respective sum z val
for (i in seq(0,24,2))
{
  if (i < 9){
    zero = "0"
  } else {
    zero = NULL
  } #Modifying variable that allows integers to be represented as "01" for values less than 10
  
  znum <- which(substr(ars[2:435,1], 6, 7) == paste(zero,i, sep = "")) # Determines the rows of each z-value
  newzdata <- ars[znum+1,2:10]/mzsumdata[(4/8)*i+1,]*100
  
  if (i == 0) {
    mdp <- newzdata
    newzdata <- NULL
  }
  
  mdp <- rbind(mdp, newzdata)
}
rownames(mdp) <- ars[2:435,1]

mdp <- na.zero(mdp) #Remove NA values from mdp
mdp <- inf.zero(mdp) #Remove Inf values from mdp

#Carbon number divided by c data sum
mcp <- mcsumdata/rep.row(msumdata, nrow(mcsumdata))*100
mcp <- data.frame(mcp)

#Z value divided by c data sum
mzp <- mzsumdata/rep.row(msumdata, nrow(mzsumdata))*100
mzp <- data.frame(mzp)


#Graphical representation of merichem data
mz <- matrix(c(mdp[1:434,2]), nrow = 13, ncol = 35, byrow = T)
rownames(z) <- seq(0, 24, by = 2)
colnames(z) <- (6:40)
mx <- (seq(0, 24, by = 2))
my <- (6:40)

colours <- c("blue3", "darkred", "green", "darkviolet", "darkslategray4", "goldenrod4", "deepskyblue3", "indianred3", "limegreen", "magenta4", "royalblue1", "salmon", "skyblue")

m <- matrix(rep(seq(13), each=35), ncol = 35, nrow = 13, byrow = TRUE)
hist3D (mx, my, mz, colvar = m, scale = T, colkey = F, 
        bty = "g", phi = 45,  theta = 45, ltheta = 60, lphi = 90, contour = FALSE,
        xlab = "Z-Value", ylab = "Carbon #", zlab = "Intensity %", main = "Oil Data: 56-IR",
        border = "black", shade = 0,
        ticktype = "detailed", space = 0.6, d = 2, cex.axis = 0.8, axis = T, nticks = 13,
        names.arg = c("0","2", "4", "6", "8", "10","12", "14", "16", "18", "20", "22", "24"))

