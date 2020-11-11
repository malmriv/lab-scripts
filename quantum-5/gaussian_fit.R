#Load pracma library for finding peaks
library(pracma)
library(errors)

files = list.files(path="./data/")
hydrogen = files[1:4]
helium = files[5:7]
mercury = files[8:12]
unknown = files[13:15]

#Read each file
data = read.csv(paste("./data/",unknown[3],sep=""),header=FALSE)
exptime = as.numeric(data[18,2])

#Prepare dataset and plot
data = data[-(1:23),]
lambda = as.numeric(data$V1)
intensity = as.numeric(data$V2)
intensity = intensity/exptime
plot(lambda,intensity,type="l",xlab="lambda (nm)",ylab="relative intensity (adim.)",
     main="Relative intensities found.")

#Smooth data in order to avoid saturated-signal issues
polinom=loess(intensity ~ lambda, span=0.02)
intensity=predict(polinom,lambda)

#Find peaks
N = 15
peaks=findpeaks(intensity,sortstr=T,npeaks=N,
                minpeakheight=max(intensity)/20)


center = peaks[,2]
start = peaks[,3]
end = peaks[,4]
for(i in 1:N) print(paste("pos = ",peaks[center[i]],", lambda =",
                          lambda[center[i]],", I = ",intensity[center[i]]))


#Perform gaussian fit
for(i in 1:N) {
  lam = lambda[start[i]:end[i]]
  inten = intensity[start[i]:end[i]]
  gaussianfit = nls(inten ~ a*exp(-(lam-b)^2/(2*c^2)),
          start=list(a=intensity[center[i]],b=lambda[center[i]],c=0.1))
  coef_a = summary(gaussianfit)$coefficients[1,1]
  errors(coef_a) = summary(gaussianfit)$coefficients[1,2]
  coef_b = summary(gaussianfit)$coefficients[2,1]
  errors(coef_b) = summary(gaussianfit)$coefficients[2,2]
  coef_c = summary(gaussianfit)$coefficients[3,1]
  errors(coef_c) = summary(gaussianfit)$coefficients[3,2]
  
  int = sqrt(2*pi)*coef_a*abs(coef_c)
  
  print(paste(lambda[center[i]],as.numeric(int),errors(int)))
  abline(v=lambda[center[i]],col="red")
 
  write.table(paste(lambda[center[i]],as.numeric(int),errors(int)),
  file="./unknown.txt",append=T,quote=F,col.names=F,row.names=F,sep="\t",dec=".")
}

normalise = read.table("./unknown.txt")
pos = which(normalise$V2 == max(normalise$V2))
normalise$V2 = normalise$V2/normalise$V2[pos]
normalise$V3 = normalise$V3/normalise$V2[pos]
write.table(normalise, file="./normalised.txt", append=T,quote=F,col.names=F,
            row.names=F,sep="\t",dec=".")




