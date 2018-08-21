library(xlsx)
library(ggplot2)
library(extrafont)
font_import()
loadfonts(device="win")       #Register fonts for Windows bitmap output
fonts()
AAprev<-read.csv("C:/Users/a/Documents/PhD/thesis/thesis doc/Thesis articles/comparing mutiple farms influenza/AA/mean prevalence values.csv")
foxcroftprev<-read.csv("C:/Users/a/Documents/PhD/thesis/thesis doc/Thesis articles/comparing mutiple farms influenza/Foxcroft/quar_0_vac_0_mean_values.csv")
JVprev<-read.csv("C:/Users/a/Documents/PhD/thesis/thesis doc/Thesis articles/comparing mutiple farms influenza/Jewel view/mean prevalence values.csv")
SHprev<-read.csv("C:/Users/a/Documents/PhD/thesis/thesis doc/Thesis articles/comparing mutiple farms influenza/Sonnenhof/quar_0_vac_0_mean_values.csv")
homog<-read.csv("C:/Users/a/Documents/PhD/thesis/thesis doc/Thesis articles/comparing mutiple farms influenza/random mixing network/quar_0_vac_0_homogenous matrix.csv")
complete<-read.csv("C:/Users/a/Documents/PhD/thesis/thesis doc/Thesis articles/comparing mutiple farms influenza/complete network/quar_0_vac_0_mean_values.csv")
minimal<-read.csv("C:/Users/a/Documents/PhD/thesis/thesis doc/Thesis articles/comparing mutiple farms influenza/network with degree 4/quar_0_vac_0_mean_values.csv")

tiff("C:/Users/a/Documents/PhD/thesis/thesis doc/Thesis articles/comparing mutiple farms influenza/comparison of prevalence with no intervention.tiff", height=180, width=180, units="mm", compression = "lzw", res=300)
###with ggplot
ggplot()+
  geom_line(data=AAprev, aes(x=time, y=i.num/num, colour="Empirical Network 1", linetype="Empirical Network 1"), size=1.5)+
  geom_line(data=foxcroftprev, aes(x=time, y=i.num/num,colour="Empirical Network 2",linetype="Empirical Network 2"), size=1.5)+
  geom_line(data=JVprev, aes(x=time, y=i.num/num,colour="Empirical Network 3", linetype="Empirical Network 3"), size=1.5)+
  geom_line(data=SHprev, aes(x=time, y=i.num/num,colour="Empirical Network 4", linetype="Empirical Network 4"), size=1.5)+
  geom_line(data=homog, aes(x=time, y=i.num/num, colour="Random Network",linetype="Random Network"), size=1.5)+
  geom_line(data=complete, aes(x=time, y=i.num/num,colour="Complete Network",linetype="Complete Network"), size=1.5)+
  geom_line(data=minimal, aes(x=time, y=i.num/num, colour="Minimal Network",linetype="Minimal Network"), size=1.5)+
  
  scale_linetype_manual("", values=c("Empirical Network 1"="solid","Empirical Network 2"="solid","Empirical Network 3"="solid","Empirical Network 4"="solid",
                                     "Random Network"="longdash","Complete Network"="dotdash","Minimal Network"="dotted"),
                        labels=c("Empirical Network 1"="Empirical Network 1","Empirical Network 2"="Empirical Network 2",
                                 "Empirical Network 3"="Empirical Network 3","Empirical Network 4"="Empirical Network 4",
                                 "Random Network"="Random Network","Complete Network"="Complete Network","Minimal Network"="Minimal Network"))+
  
  scale_colour_manual("",values = c("Empirical Network 1"="red","Empirical Network 2"="blue","Empirical Network 3"="green","Empirical Network 4"="orange",
                                    "Random Network"="black","Complete Network"="black","Minimal Network"="black"),
                      labels=c("Empirical Network 1"="Empirical Network 1","Empirical Network 2"="Empirical Network 2",
                               "Empirical Network 3"="Empirical Network 3","Empirical Network 4"="Empirical Network 4",
                               "Random Network"="Random Network","Complete Network"="Complete Network","Minimal Network"="Minimal Network"))+
  theme_classic()%+replace% 
  theme(legend.key.width = unit(2,"cm"),
        axis.text = element_text(size = 18, family="Times New Roman"),
        axis.title=element_text(size=18, family="Times New Roman"),
        legend.text=element_text(size=18, family="Times New Roman"),
        # legend.position = c(.8,.5))
        legend.position="none")+
  ylab("Proportion of  infectious horses")+
  xlab("Time (days)")+
  # scale_x_continuous(limits = c(0, 30))+
  scale_x_continuous(limits=c(0,30))+
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6))
  dev.off()







tiff("C:/Users/a/Documents/PhD/thesis/thesis doc/Thesis articles/comparing mutiple farms influenza/comparison of incidence with no intervention.tiff", height=180, width=180, units="mm", compression = "lzw", res=300)
###with ggplot
ggplot()+
  geom_line(data=AAprev,       aes(x=time, y=se.flow/num, colour="Empirical Network 1", linetype="Empirical Network 1"), size=1.5)+
  geom_line(data=foxcroftprev, aes(x=time, y=se.flow/num,colour="Empirical Network 2",linetype="Empirical Network 2"), size=1.5)+
  geom_line(data=JVprev,       aes(x=time, y=se.flow/num,colour="Empirical Network 3", linetype="Empirical Network 3"), size=1.5)+
  geom_line(data=SHprev,       aes(x=time, y=se.flow/num,colour="Empirical Network 4", linetype="Empirical Network 4"), size=1.5)+
  geom_line(data=homog,        aes(x=time, y=se.flow/num, colour="Random Network",linetype="Random Network"), size=1.5)+
  geom_line(data=complete,     aes(x=time, y=se.flow/num,colour="Complete Network",linetype="Complete Network"), size=1.5)+
  geom_line(data=minimal,      aes(x=time, y=se.flow/num, colour="Minimal Network",linetype="Minimal Network"), size=1.5)+
  
  scale_linetype_manual("", values=c("Empirical Network 1"="solid","Empirical Network 2"="solid","Empirical Network 3"="solid","Empirical Network 4"="solid",
                                 "Random Network"="longdash","Complete Network"="dotdash","Minimal Network"="dotted"),
                        labels=c("Empirical Network 1"="Empirical Network 1","Empirical Network 2"="Empirical Network 2",
                                 "Empirical Network 3"="Empirical Network 3","Empirical Network 4"="Empirical Network 4",
                                "Random Network"="Random Network","Complete Network"="Complete Network","Minimal Network"="Minimal Network"))+
  
  scale_colour_manual("",values = c("Empirical Network 1"="red","Empirical Network 2"="blue","Empirical Network 3"="green","Empirical Network 4"="orange",
                                 "Random Network"="black","Complete Network"="black","Minimal Network"="black"),
                      labels=c("Empirical Network 1"="Empirical Network 1","Empirical Network 2"="Empirical Network 2",
                               "Empirical Network 3"="Empirical Network 3","Empirical Network 4"="Empirical Network 4",
                               "Random Network"="Random Network","Complete Network"="Complete Network","Minimal Network"="Minimal Network"))+
  theme_classic()%+replace% 
  theme(legend.key.width = unit(2,"cm"),
        axis.title=element_text(size=18, family="Times New Roman"),
        axis.text = element_text(size = 18, family="Times New Roman"),
        legend.text=element_text(size=18, family="Times New Roman"),
        legend.position = c(.8,.5)
        )+
  ylab("Proportion of new infectious horses")+
  xlab("Time (days)")+
  scale_x_continuous(limits=c(0,13))+
  # scale_x_continuous(breaks=c(0,5,10,15,20))
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5, 0.6,0.7,0.8,0.9,1))
 
  # scale_y_continuous(limits = c(0, .6))
dev.off()

#################################################3
#with plot
# plot((i.num/num)~time,data=AAprev,type="l",col="red",lwd=2, ylim=range(0 :0.6), xlim=range(0:15), xlab="Time (days)", ylab="Proportion of infectious horses")
# lines((i.num/num)~time,data=foxcroftprev,type="l",lty="solid",col="blue",lwd=2)
# lines((i.num/num)~time,data=SHprev,type="l",lty="solid",lwd=2, col="green")
# lines((i.num/num)~time,data=JVprev,type="l",lty="solid",col="orange",lwd=2)
# lines((i.num/num)~time,data=homog,type="l",lty="longdash",col="black",lwd=2)
# lines((i.num/num)~time,data=complete,type="l",lty="dotdash",col="black",lwd=2)
# lines((i.num/num)~time,data=minimal,type="l",lty="dotted",col="black",lwd=2)
# 
# legend("topright",legend=c("Empirical network 1","Empirical network 2","Empirical network 3","Empirical network 4","Random mixing network",
#                            "Complete network","Minimal network"),
#        lty=c("solid","solid","solid","solid","longdash","dotdash","dotted"),
#        col=c("red","blue","green","orange","black","black","black"),
#        lwd=c(2,2,2,2,2,2,2))
# 
# dev.off()

























