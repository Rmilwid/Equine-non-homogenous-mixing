# plot % reduction for different intervention strategies
#import prevalence lists
setwd("C:/Users/a/Documents/PhD/thesis/thesis doc/Thesis articles/comparing mutiple farms influenza/homogenous mixing matrix")
quar0vac0<-read.csv("quar_0_vac_0_homogenous matrix.csv")
quar0vac25<-read.csv("quar_0_vac_25_homogenous matrix.csv")
quar0vac50<-read.csv("quar_0_vac_50_homogenous matrix.csv")
quar0vac75<-read.csv("quar_0_vac_75_homogenous matrix.csv")
quar0vac100<-read.csv("quar_0_vac_100_homogenous matrix.csv")
quar25vac0<-read.csv("quar_25_vac_0_homogenous matrix.csv")
quar25vac25<-read.csv("quar_25_vac_25_homogenous matrix.csv")
quar25vac50<-read.csv("quar_25_vac_50_homogenous matrix.csv")
quar25vac75<-read.csv("quar_25_vac_75_homogenous matrix.csv")
quar25vac100<-read.csv("quar_25_vac_100_homogenous matrix.csv")
quar50vac0<-read.csv("quar_50_vac_0_homogenous matrix.csv")
quar50vac25<-read.csv("quar_50_vac_25_homogenous matrix.csv")
quar50vac50<-read.csv("quar_50_vac_50_homogenous matrix.csv")
quar50vac75<-read.csv("quar_50_vac_75_homogenous matrix.csv")
quar50vac100<-read.csv("quar_50_vac_100_homogenous matrix.csv")
quar75vac0<-read.csv("quar_75_vac_0_homogenous matrix.csv")
quar75vac25<-read.csv("quar_75_vac_25_homogenous matrix.csv")
quar75vac50<-read.csv("quar_75_vac_50_homogenous matrix.csv")
quar75vac75<-read.csv("quar_75_vac_75_homogenous matrix.csv")
quar75vac100<-read.csv("quar_75_vac_100_homogenous matrix.csv")
quar100vac0<-read.csv("quar_100_vac_0_homogenous matrix.csv")
quar100vac25<-read.csv("quar_100_vac_25_homogenous matrix.csv")
quar100vac50<-read.csv("quar_100_vac_50_homogenous matrix.csv")
quar100vac75<-read.csv("quar_100_vac_75_homogenous matrix.csv")
quar100vac100<-read.csv("quar_100_vac_100_homogenous matrix.csv")

# #percet reuction at vac rate i=I_i/I_0
# tiff("percent reduction constant quarantine.tiff", height=18, width= 18, units= "cm", compression = "lzw", res=300)
# par(mfrow=c(3,2))
# plot(1-(quar0vac25$i.num/quar0vac0$i.num), type="l", col="red",ylim=range(0:1),lwd=2,xlim=range(0:20),
#      lty=1,xlab="Time (days)",ylab="Percent reduction", main="Proportion quarantined =0%")
# lines(1-(quar0vac50$i.num/quar0vac0$i.num), type="l", col="blue",lty=2,lwd=2)
# lines(1-(quar0vac75$i.num/quar0vac0$i.num), type="l", col="green",lty=3,lwd=2)
# lines(1-(quar0vac100$i.num/quar0vac0$i.num), type="l", col="black",lty=4,lwd=2)
# grid(NULL,NULL,col="black")
# 
# plot(1-(quar25vac25$i.num/quar25vac0$i.num), type="l", col="red",ylim=range(0:1),xlim=range(0:20),
#      lty=1,lwd=2,xlab="Time (days)",ylab="Percent reduction", main="Proportion quarantined =25%")
# lines(1-(quar25vac50$i.num/quar25vac0$i.num), type="l", col="blue",lty=2,lwd=2)
# lines(1-(quar25vac75$i.num/quar25vac0$i.num), type="l", col="green",lty=3,lwd=2)
# lines(1-(quar25vac100$i.num/quar25vac0$i.num), type="l", col="black",lty=4,lwd=2)
# grid(NULL, NULL,col="black")
# 
# plot(1-(quar50vac25$i.num/quar50vac0$i.num), type="l", col="red",ylim=range(0:1),xlim=range(0:20),
#      lty=1,lwd=2,xlab="Time (days)",ylab="Percent reduction", main="Proportion quarantined =50%")
# lines(1-(quar50vac50$i.num/quar50vac0$i.num), type="l", col="blue",lty=2,lwd=2)
# lines(1-(quar50vac75$i.num/quar50vac0$i.num), type="l", col="green",lty=3,lwd=2)
# lines(1-(quar50vac100$i.num/quar50vac0$i.num), type="l", col="black",lty=4,lwd=2)
# grid(NULL, NULL,col="black")
# 
# plot(1-(quar75vac25$i.num/quar75vac0$i.num), type="l", col="red",ylim=range(0:1),xlim=range(0:20),
#      lty=1,lwd=2,xlab="Time (days)",ylab="Percent reduction", main="Proportion quarantined =75%")
# lines(1-(quar75vac50$i.num/quar75vac0$i.num), type="l", col="blue",lty=2,lwd=2)
# lines(1-(quar75vac75$i.num/quar75vac0$i.num), type="l", col="green",lty=3,lwd=2)
# lines(1-(quar75vac100$i.num/quar75vac0$i.num), type="l", col="black",lty=4,lwd=2)
# grid(NULL, NULL,col="black")
# 
# plot(1-(quar100vac25$i.num/quar100vac0$i.num), type="l", col="red",ylim=range(0:1),xlim=range(0:20),
#      lty=1,lwd=2,xlab="Time (days)",ylab="Percent reduction", main="Proportion quarantined =100%")
# lines(1-(quar100vac50$i.num/quar100vac0$i.num), type="l", col="blue",lty=2,lwd=2)
# lines(1-(quar100vac75$i.num/quar100vac0$i.num), type="l", col="green",lty=3,lwd=2)
# lines(1-(quar100vac100$i.num/quar100vac0$i.num), type="l", col="black",lty=4,lwd=2)
# grid(NULL, NULL,col="black")
# 
# plot(0,type='n',axes=FALSE,ann=FALSE)
# legend("topleft",legend=c("vac=25","vac=50","vac=75","vac=100"),col=c("red","blue","green","black"),lty=c(1,2,3,4),lwd=c(2,2,2,2))
# 
# dev.off()
# #percet reuction at quar rate i=I_i/I_0
# tiff("percent reduction constant vac.tiff", height=18, width= 18, units= "cm", compression = "lzw", res=300)
# par(mfrow=c(3,2))
# 
# 
# plot(1-(quar25vac0$i.num/quar0vac0$i.num), type="l", col="red",ylim=range(0:1),xlim=range(0:20),
#      lty=1,lwd=2,xlab="Time (days)",ylab="Percent reduction", main="Proportion vaccinated = 0%")
# lines(1-(quar50vac0$i.num/quar0vac0$i.num), type="l", col="blue",lty=2,lwd=2)
# lines(1-(quar75vac0$i.num/quar0vac0$i.num), type="l", col="green",lty=3,lwd=2)
# lines(1-(quar100vac0$i.num/quar0vac0$i.num), type="l", col="black",lty=4,lwd=2)
# grid(NULL, NULL,col="black")
# 
# plot(1-(quar25vac25$i.num/quar0vac25$i.num), type="l", col="red",ylim=range(0:1),xlim=range(0:20),
#      lty=1,lwd=2,xlab="Time (days)",ylab="Percent reduction", main="Proportion vaccinated = 25%")
# lines(1-(quar50vac25$i.num/quar0vac25$i.num), type="l", col="blue",lty=2,lwd=2)
# lines(1-(quar75vac25$i.num/quar0vac25$i.num), type="l", col="green",lty=3,lwd=2)
# lines(1-(quar100vac25$i.num/quar0vac25$i.num), type="l", col="black",lty=4,lwd=2)
# grid(NULL, NULL,col="black")
# 
# plot(1-(quar25vac50$i.num/quar0vac50$i.num), type="l", col="red",ylim=range(0:1),xlim=range(0:20),
#      lty=1,lwd=2,xlab="Time (days)",ylab="Percent reduction", main="Proportion vaccinated = 50%")
# lines(1-(quar50vac50$i.num/quar0vac50$i.num), type="l", col="blue",lty=2,lwd=2)
# lines(1-(quar75vac50$i.num/quar0vac50$i.num), type="l", col="green",lty=3,lwd=2)
# lines(1-(quar100vac50$i.num/quar0vac50$i.num), type="l", col="black",lty=4,lwd=2)
# grid(NULL, NULL,col="black")
# 
# plot(1-(quar25vac75$i.num/quar0vac75$i.num), type="l", col="red",ylim=range(0:1),xlim=range(0:20),
#      lty=1,lwd=2,xlab="Time (days)",ylab="Percent reduction", main="Proportion vaccinated = 75%")
# lines(1-(quar50vac75$i.num/quar0vac75$i.num), type="l", col="blue",lty=2,lwd=2)
# lines(1-(quar75vac75$i.num/quar0vac75$i.num), type="l", col="green",lty=3,lwd=2)
# lines(1-(quar100vac75$i.num/quar0vac75$i.num), type="l", col="black",lty=4,lwd=2)
# grid(NULL, NULL,col="black")
# plot(1-(quar25vac100$i.num/quar0vac100$i.num), type="l", col="red",ylim=range(0:1),xlim=range(0:20),
#      lty=1,lwd=2,xlab="Time (days)",ylab="Percent reduction", main="Proportion vaccinated = 100%")
# lines(1-(quar50vac100$i.num/quar0vac100$i.num), type="l", col="blue",lty=2,lwd=2)
# lines(1-(quar75vac100$i.num/quar0vac100$i.num), type="l", col="green",lty=3,lwd=2)
# lines(1-(quar100vac100$i.num/quar0vac100$i.num), type="l", col="black",lty=4,lwd=2)
# grid(NULL, NULL,col="black")
# 
# plot(0,type='n',axes=FALSE,ann=FALSE)
# legend("topleft",legend=c("quar=25","quar=50","quar=75","quar=100"),col=c("red","blue","green","black"),lty=c(1,2,3,4),lwd=c(2,2,2,2))
# 
# dev.off()

#regular prevalence plots
tiff("plot constant quarantine.tiff", height=18, width= 18, units= "cm", compression = "lzw", res=300)
par(mfrow=c(3,2))
plot(quar0vac0$i.num, type="l", col="red",ylim=range(0:15),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion quarantined =0%")
lines(quar0vac25$i.num, type="l", col="pink",lty=2,lwd=2)
lines(quar0vac50$i.num, type="l", col="blue",lty=3,lwd=2)
lines(quar0vac75$i.num, type="l", col="green",lty=4,lwd=2)
lines(quar0vac100$i.num, type="l", col="black",lty=5,lwd=2)

plot(quar25vac0$i.num, type="l", col="red",ylim=range(0:15),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion quarantined =25%")
lines(quar25vac25$i.num, type="l", col="pink",lty=2,lwd=2)
lines(quar25vac50$i.num, type="l", col="blue",lty=3,lwd=2)
lines(quar25vac75$i.num, type="l", col="green",lty=4,lwd=2)
lines(quar25vac100$i.num, type="l", col="black",lty=5,lwd=2)

plot(quar50vac0$i.num, type="l", col="red",ylim=range(0:15),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion quarantined =50%")
lines(quar50vac25$i.num, type="l", col="pink",lty=2,lwd=2)
lines(quar50vac50$i.num, type="l", col="blue",lty=3,lwd=2)
lines(quar50vac75$i.num, type="l", col="green",lty=4,lwd=2)
lines(quar50vac100$i.num, type="l", col="black",lty=5,lwd=2)

plot(quar75vac0$i.num, type="l", col="red",ylim=range(0:15),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion quarantined =75%")
lines(quar75vac25$i.num, type="l", col="pink",lty=2,lwd=2)
lines(quar75vac50$i.num, type="l", col="blue",lty=3,lwd=2)
lines(quar75vac75$i.num, type="l", col="green",lty=4,lwd=2)
lines(quar75vac100$i.num, type="l", col="black",lty=5,lwd=2)

plot(quar100vac0$i.num, type="l", col="red",ylim=range(0:15),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion quarantined =100%")
lines(quar100vac25$i.num, type="l", col="pink",lty=2,lwd=2)
lines(quar100vac50$i.num, type="l", col="blue",lty=3,lwd=2)
lines(quar100vac75$i.num, type="l", col="green",lty=4,lwd=2)
lines(quar100vac100$i.num, type="l", col="black",lty=5,lwd=2)

plot(0,type='n',axes=FALSE,ann=FALSE)
legend("topleft",legend=c("vac=0","vac=25","vac=50","vac=75","vac=100"),col=c("red","pink","blue","green","black"),lty=c(1,2,3,4,5),lwd=c(2,2,2,2,2))

dev.off()

tiff("plot constant vac.tiff", height=18, width= 18, units= "cm", compression = "lzw", res=300)
par(mfrow=c(3,2))


plot(quar0vac0$i.num, type="l", col="red",ylim=range(0:15),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion vaccinated = 0%")
lines(quar25vac0$i.num, type="l", col="pink",lty=2,lwd=2)
lines(quar50vac0$i.num, type="l", col="blue",lty=3,lwd=2)
lines(quar75vac0$i.num, type="l", col="green",lty=4,lwd=2)
lines(quar100vac0$i.num, type="l", col="black",lty=5,lwd=2)

plot(quar0vac25$i.num, type="l", col="red",ylim=range(0:15),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion vaccinated = 25%")
lines(quar25vac25$i.num, type="l", col="pink",lty=2,lwd=2)
lines(quar50vac25$i.num, type="l", col="blue",lty=3,lwd=2)
lines(quar75vac25$i.num, type="l", col="green",lty=4,lwd=2)
lines(quar100vac25$i.num, type="l", col="black",lty=5,lwd=2)

plot(quar0vac50$i.num, type="l", col="red",ylim=range(0:15),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion vaccinated = 50%")
lines(quar25vac50$i.num, type="l", col="pink",lty=2,lwd=2)
lines(quar50vac50$i.num, type="l", col="blue",lty=3,lwd=2)
lines(quar75vac50$i.num, type="l", col="green",lty=4,lwd=2)
lines(quar100vac50$i.num, type="l", col="black",lty=5,lwd=2)

plot(quar0vac75$i.num, type="l", col="red",ylim=range(0:15),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion vaccinated = 75%")
lines(quar25vac75$i.num, type="l", col="pink",lty=2,lwd=2)
lines(quar50vac75$i.num, type="l", col="blue",lty=3,lwd=2)
lines(quar75vac75$i.num, type="l", col="green",lty=4,lwd=2)
lines(quar100vac75$i.num, type="l", col="black",lty=5,lwd=2)

plot(quar0vac100$i.num, type="l", col="red",ylim=range(0:15),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion vaccinated = 100%")
lines(quar25vac100$i.num, type="l", col="pink",lty=2,lwd=2)
lines(quar50vac100$i.num, type="l", col="blue",lty=3,lwd=2)
lines(quar75vac100$i.num, type="l", col="green",lty=4,lwd=2)
lines(quar100vac100$i.num, type="l", col="black",lty=5,lwd=2)

plot(0,type='n',axes=FALSE,ann=FALSE)
legend("topleft",legend=c("quar=0","quar=25","quar=50","quar=75","quar=100"),col=c("red","pink","blue","green","black"),lty=c(1,2,3,4,5),lwd=c(2,2,2,2,2))
dev.off()

####################3
#regular incidence plots


tiff("plot incidence constant quarantine.tiff", height=18, width= 18, units= "cm", compression = "lzw", res=300)
par(mfrow=c(3,2))
plot(quar0vac0$se.flow, type="l", col="red",ylim=range(0:13),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion quarantined =0%")
lines(quar0vac25$se.flow, type="l", col="pink",lty=2,lwd=2)
lines(quar0vac50$se.flow, type="l", col="blue",lty=3,lwd=2)
lines(quar0vac75$se.flow, type="l", col="green",lty=4,lwd=2)
lines(quar0vac100$se.flow, type="l", col="black",lty=5,lwd=2)

plot(quar25vac0$se.flow, type="l", col="red",ylim=range(0:13),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion quarantined =25%")
lines(quar25vac25$se.flow, type="l", col="pink",lty=2,lwd=2)
lines(quar25vac50$se.flow, type="l", col="blue",lty=3,lwd=2)
lines(quar25vac75$se.flow, type="l", col="green",lty=4,lwd=2)
lines(quar25vac100$se.flow, type="l", col="black",lty=5,lwd=2)

plot(quar50vac0$se.flow, type="l", col="red",ylim=range(0:13),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion quarantined =50%")
lines(quar50vac25$se.flow, type="l", col="pink",lty=2,lwd=2)
lines(quar50vac50$se.flow, type="l", col="blue",lty=3,lwd=2)
lines(quar50vac75$se.flow, type="l", col="green",lty=4,lwd=2)
lines(quar50vac100$se.flow, type="l", col="black",lty=5,lwd=2)

plot(quar75vac0$se.flow, type="l", col="red",ylim=range(0:13),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion quarantined =75%")
lines(quar75vac25$se.flow, type="l", col="pink",lty=2,lwd=2)
lines(quar75vac50$se.flow, type="l", col="blue",lty=3,lwd=2)
lines(quar75vac75$se.flow, type="l", col="green",lty=4,lwd=2)
lines(quar75vac100$se.flow, type="l", col="black",lty=5,lwd=2)

plot(quar100vac0$se.flow, type="l", col="red",ylim=range(0:13),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion quarantined =100%")
lines(quar100vac25$se.flow, type="l", col="pink",lty=2,lwd=2)
lines(quar100vac50$se.flow, type="l", col="blue",lty=3,lwd=2)
lines(quar100vac75$se.flow, type="l", col="green",lty=4,lwd=2)
lines(quar100vac100$se.flow, type="l", col="black",lty=5,lwd=2)

plot(0,type='n',axes=FALSE,ann=FALSE)
legend("topleft",legend=c("vac=0","vac=25","vac=50","vac=75","vac=100"),col=c("red","pink","blue","green","black"),lty=c(1,2,3,4,5),lwd=c(2,2,2,2,2))

dev.off()


tiff("plot incidence constant vac.tiff", height=18, width= 18, units= "cm", compression = "lzw", res=300)
par(mfrow=c(3,2))


plot(quar0vac0$se.flow, type="l", col="red",ylim=range(0:13),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion vaccinated = 0%")
lines(quar25vac0$se.flow, type="l", col="pink",lty=2,lwd=2)
lines(quar50vac0$se.flow, type="l", col="blue",lty=3,lwd=2)
lines(quar75vac0$se.flow, type="l", col="green",lty=4,lwd=2)
lines(quar100vac0$se.flow, type="l", col="black",lty=5,lwd=2)

plot(quar0vac25$se.flow, type="l", col="red",ylim=range(0:13),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion vaccinated = 25%")
lines(quar25vac25$se.flow, type="l", col="pink",lty=2,lwd=2)
lines(quar50vac25$se.flow, type="l", col="blue",lty=3,lwd=2)
lines(quar75vac25$se.flow, type="l", col="green",lty=4,lwd=2)
lines(quar100vac25$se.flow, type="l", col="black",lty=5,lwd=2)

plot(quar0vac50$se.flow, type="l", col="red",ylim=range(0:13),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion vaccinated = 50%")
lines(quar25vac50$se.flow, type="l", col="pink",lty=2,lwd=2)
lines(quar50vac50$se.flow, type="l", col="blue",lty=3,lwd=2)
lines(quar75vac50$se.flow, type="l", col="green",lty=4,lwd=2)
lines(quar100vac50$se.flow, type="l", col="black",lty=5,lwd=2)

plot(quar0vac75$se.flow, type="l", col="red",ylim=range(0:13),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion vaccinated = 75%")
lines(quar25vac75$se.flow, type="l", col="pink",lty=2,lwd=2)
lines(quar50vac75$se.flow, type="l", col="blue",lty=3,lwd=2)
lines(quar75vac75$se.flow, type="l", col="green",lty=4,lwd=2)
lines(quar100vac75$se.flow, type="l", col="black",lty=5,lwd=2)

plot(quar0vac100$se.flow, type="l", col="red",ylim=range(0:13),xlim=range(0:17),
     lty=1,lwd=2,xlab="Time (days)",ylab="Number infectious", main="Proportion vaccinated = 100%")
lines(quar25vac100$se.flow, type="l", col="pink",lty=2,lwd=2)
lines(quar50vac100$se.flow, type="l", col="blue",lty=3,lwd=2)
lines(quar75vac100$se.flow, type="l", col="green",lty=4,lwd=2)
lines(quar100vac100$se.flow, type="l", col="black",lty=5,lwd=2)

plot(0,type='n',axes=FALSE,ann=FALSE)
legend("topleft",legend=c("quar=0","quar=25","quar=50","quar=75","quar=100"),col=c("red","pink","blue","green","black"),lty=c(1,2,3,4,5),lwd=c(2,2,2,2,2))
dev.off()
