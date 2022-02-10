##Pre-requisite 1: all the scripts should be in the directory Desktop/ScriptsR_modeles
##Pre-requisite 2: all the data should be in the directory Desktop/ScriptsR_modeles/Data
##Advide: Pull the directory bet-hedging and move ScriptsR_modeles to your desktop

library(stringr)
#install.packages("deming")
library(deming)

par(mar=c(4,4,3,1))
#Analysis of Bar-Even11 data
##I Importing the global data
setwd(dir="~")
setwd(dir="/Users/florian/enzyme-evolution/Data")
enzyme_data<-read.csv(file="BarEven11.txt",sep="\t",header=TRUE)#,header=TRUE)
enzyme_data$kcat<-as.numeric(enzyme_data$kcat)
head(enzyme_data)

#1.Testing for a trade-off (non log) between kcat and KM
toff01<-lm(enzyme_data$kcat~enzyme_data$KM)
summary(toff01)#Conclusion: No linear relationship R^2~0, although significant light increase
toff02<-lm(enzyme_data$KM~enzyme_data$kcat)
summary(toff02)
plot(enzyme_data$KM,enzyme_data$kcat)

#Plottng with proxies of pseudo-mechanistic parameters
kcat_KM<-enzyme_data$kcat/enzyme_data$KM
plot(kcat_KM,enzyme_data$kcat)

#par(mfrow=c(1,1))#Redefining plot prameter from multiple plots used previously

#2.Testing trade-offs with log-scaling
#Transforming the data
enzyme_data$log10kcat<-log10(enzyme_data$kcat)
enzyme_data$log10KM<-log10(enzyme_data$KM)
enzyme_data$log10kcat_KM<-log10(enzyme_data$kcat/enzyme_data$KM*10^6)#*10^6 to visualize in MicroMolar

##a.kcat vs KM trade-off on log-scale
###y-linear regression
par(mfrow=c(1,3),mai=c(0.75,0.75,0.75,0.75))
toff1<-lm(enzyme_data$log10kcat~enzyme_data$log10KM)
toff2<-lm(enzyme_data$log10KM~enzyme_data$log10kcat)
summary(toff1)#Conclusion: Light apprent trade-off with linear relationship, R^2=0.09
addtxt<-list(l=0.05,h=0.95,txt=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"),srt = 0,font=2,col="black")
plot(enzyme_data$log10KM~enzyme_data$log10kcat,pch=20,col="tomato",
     xlab="log10(kcat)",
     ylab="log10(KM)",cex=0.75,cex.lab=1.5,cex.axis=1.5)
abline(coef=c(toff2$coefficients[1],toff2$coefficients[2]),lwd=2,col="black")
#Conclusion: significant relationship, R^2=0.09, but low explained variance
#especially accounting for the fact that there should be a linear relatinoship as KM=(kcat+kr)/kf=kr/kf+(1/kf)*kcat
#This means that ln(KM)=ln(1/kf)+ln(kr+kcat)=f(ln(kr+kcat)); quality of the trade-off depends on the weight of kcat compared to that of kr and the codependance of kf and kr with kcat 

###ortho-linear regression
#toff_ortho1 <- deming(enzyme_data$log10KM~enzyme_data$log10kcat)#with package deming
abline(coef=c(toff_ortho1$coefficients[1],toff_ortho1$coefficients[2]),lwd=2,col="black",lty=2)
text(-5,6.5,"A",srt=addtxt$srt,font=2,col="black",cex=1.5)
legend("bottomright",title="Regression",legend=c("Simple linear","Deming"),lty=c(1,2),cex=1.25)
##a.kcat vs KM trade-off on log-scale

toff3<-lm(enzyme_data$log10KM~enzyme_data$log10kcat_KM)
plot(enzyme_data$log10KM~enzyme_data$log10kcat_KM,pch=20,col="tomato",
     xlab="log10(kcat/KM)",
     ylab="log10(KM)",cex=0.75,cex.lab=1.5,cex.axis=1.5)
#toff_ortho3 <- deming(enzyme_data$log10KM~enzyme_data$log10kcat_KM)#with package deming
abline(coef=c(toff3$coefficients[1],toff3$coefficients[2]),lwd=2,col="black")
abline(coef=c(toff_ortho3$coefficients[1],toff_ortho3$coefficients[2]),lwd=2,col="black",lty=2)
text(-1,6.5,"B",srt=addtxt$srt,font=2,col="black",cex=1.5)
legend("topright",title="Regression",legend=c("Simple linear","Deming"),lty=c(1,2),cex=1.25)

toff2<-lm(enzyme_data$log10kcat~enzyme_data$log10kcat_KM)
plot(enzyme_data$log10kcat~enzyme_data$log10kcat_KM,pch=20,col="tomato",
     xlab="log10(kcat/KM)",
     ylab="log10(kcat)",cex=0.75,cex.lab=1.5,cex.axis=1.5)
#toff_ortho2 <- deming(enzyme_data$log10kcat~enzyme_data$log10kcat_KM)#with package deming
abline(coef=c(toff2$coefficients[1],toff2$coefficients[2]),lwd=2,col="black")
abline(coef=c(toff_ortho2$coefficients[1],toff_ortho2$coefficients[2]),lwd=2,col="black",lty=2)
text(-1,5.6,"C",srt=addtxt$srt,font=2,col="black",cex=1.5)
legend("bottomright",title="Regression",legend=c("Simple linear","Deming"),lty=c(1,2),cex=1.25)

setwd(dir="~")
setwd(dir="Desktop/Ongoing-projects/")
dev.print(device = jpeg, file = "Bar-Even-tot-data.jpeg", , width = 1350*3,height=450*3,res=300,type="cairo")
dev.off()  
#Conclusion: significant positive relationship R^2=0.408, more significant than the previous one
#which is a strong evidence for the compensation of the increase in kcat in order to
#limit the decrease in affinity. There is no apparent trade-off between kf and kcat, although there might
#exist a limitation to its optimization masked by the many non optimized enzymes (a phenomenon which remains to be explained)
#Interpretation: need to increase kcat while containing the increase of KM. Possible influence of enzyme concentration for the "height" in the graph

##b.kcat vs kcat/KM on log-scale
#Note: this should represent an approximation of log10(kcat) versus log10(kf)
#as kcat/KM=(kcat/kcat+kr)*kf=(1-kr/kcat)*kf
toff_ortho2 <- deming(enzyme_data$log10kcat~enzyme_data$log10kcat_KM)
print(toff_ortho2)
abline(coef=c(toff_ortho2$coefficients[1],toff_ortho2$coefficients[2]),lwd=2,col="red")
#Same pattern as observed for kcat versus KM: effect seems to be higher, reality surely lying in between

#IIAnalysis of segmented data
#0.Opening the multiple data files #Need to be put in the same file #setwd(dir="Desktop/ScriptsR_modeles/Data")
setwd(dir="~")
setwd(dir="/Users/florian/enzyme-evolution/Data")
module_data<-read.csv(file="ID_modules.csv",sep=";",header=TRUE)
reaction_data<-read.csv(file="ID_reaction.csv",sep=";",header=TRUE)
enzyme_data<-read.csv(file="KM_kcat.csv",sep=";",header=TRUE)
organism_data<-read.csv(file="Organism_ID.csv",sep=";",header=TRUE)
#Renaming the columns
colnames(enzyme_data)<-c("Reaction_ID","Organism_ID","KM(muM)","kcat(1/s)")
colnames(organism_data)<-c("Organism_ID","Organism_name")
#Visualizing tables
head(module_data)
head(reaction_data)
head(enzyme_data)
head(organism_data)

#a.Attributing module type to reaction
module_Type<-c()
for (r in 1:length(reaction_data$Reaction_ID)){
        module_Type<-c(module_Type,
                       toString(module_data$Type[module_data$ID_MOD==reaction_data$Module_ID[r]]))
}
reaction_data<-cbind(reaction_data,module_Type)

#b.Attributing module type to enzyme
module_Type<-c()
for (r in 1:length(enzyme_data$Reaction_ID)){
        module_Type<-c(module_Type,
                       toString(reaction_data$module_Type[reaction_data$Reaction_ID==enzyme_data$Reaction_ID[r]]))
}

enzyme_data<-cbind(enzyme_data,module_Type)

#c.Removing data if module type unknown
enzyme_data<-enzyme_data[enzyme_data$module_Type!="",]

##c.Attributing organism name to enzyme 
organism_ID<-c()
organism_Name<-c()

for (s in 1:length(enzyme_data$Organism_ID)){
        organism_ID<-c(organism_ID,enzyme_data$Organism_ID[s])
        organism_Name<-c(organism_Name,
                         toString(organism_data[organism_data$Organism_ID==toString(enzyme_data$Organism_ID[s]),]$Organism_name[1]))
}

enzyme_data<-cbind(enzyme_data,organism_Name)

##c.For enzymes involved in many modules, picking for the putatively most demanding (1.CEM, 2.AAFAN, 3.INTER,4.SEC) 
for (s in 1:length(enzyme_data$Organism_ID)){
        if(str_detect(enzyme_data$module_Type[s],"CEM")){
                enzyme_data$module_Type[s]="CEM"
        }
        else if(str_detect(enzyme_data$module_Type[s],"AA,FA,N")){
                enzyme_data$module_Type[s]="AA,FA,N"
        }
        else if(str_detect(enzyme_data$module_Type[s],"INTER")){
                enzyme_data$module_Type[s]="INTER"
        }
        else if(str_detect(enzyme_data$module_Type[s],"SEC")){
                enzyme_data$module_Type[s]="SEC"
        }
}


#1.kcat vs KM trade-off
##a.Analysing data with regards to categories
###With regards to Metabolic modules

#log transforming the data
log10kcat<-c()
log10KM<-c()
log10kcat_KM<-c()
for (s in 1:length(enzyme_data$Organism_ID)){
        log10kcat<-c(log10kcat,log(as.numeric(as.character(enzyme_data$`kcat(1/s)`[s])),10))
        log10KM<-c(log10KM,log(as.numeric(as.character(enzyme_data$`KM(muM)`[s])),10))
        log10kcat_KM<-c(log10kcat_KM,6+log(as.numeric(as.character(enzyme_data$`kcat(1/s)`[s]))/as.numeric(as.character(enzyme_data$`KM(muM)`[s])),10))
}
enzyme_data<-cbind(enzyme_data,log10kcat,log10KM,log10kcat_KM)

#Picking the most significant trade-off
toff_SEC<-lm(enzyme_data$log10kcat[enzyme_data$module_Type=="SEC"]~enzyme_data$log10KM[enzyme_data$module_Type=="SEC"])
summary(toff_SEC)#R^2=0.28, quite explainatory for this category
##Most significant trade-off seems for the lowest demanding family of enzymes
##It seems perfectly in line with the fact the there is an obvious trade-off between the two quantities
##but one that can be opposed by a concomitant increase in both kcat and kf

#Opening a plot with the appropriate scaling
plot(enzyme_data$log10kcat~enzyme_data$log10KM,pch=20,col="white",
     xlab="log10(KM)",
     ylab="log10(kcat)",
     main="Trade-off between kcat and KM")

col_graph<-c()
leg_graph<-c()
j=1
for (i in levels(factor(enzyme_data$module_Type))){
        points(enzyme_data$log10kcat[enzyme_data$module_Type==i]~enzyme_data$log10KM[enzyme_data$module_Type==i],pch=j-1,col=j,cex=0.75)
        col_graph=c(col_graph,j)
        leg_graph<-c(leg_graph,i)
        j=j+1
}
legend("topleft", legend=leg_graph,
       col=col_graph, pch=20, cex=0.8)

abline(coef=c(toff_SEC$coefficients[1],toff_SEC$coefficients[2]),lwd=2,col=4)
#Conclusion: the difference exists between categories although it does not seem very significant
#Does not seem to need a more detailed analysis

##b.With regards to Species

#Extracting most represented species, i.e. those with at least 8 enzymes represented in the data
species_table<-c()
for (t in 1:length(table(enzyme_data$organism_Name))){
        if(table(enzyme_data$organism_Name)[t]>7)
        species_table<-c(species_table,names(table(enzyme_data$organism_Name)[t]))
}

#Opening a plot with the appropriate scaling
plot(enzyme_data$log10kcat~enzyme_data$log10KM,pch=20,col="white",
     xlab="log10(KM)",
     ylab="log10(kcat)",
     main="Trade-off between kcat and KM")

col_graph<-c()
leg_graph<-c()

for (i in species_table){
        points(enzyme_data$log10kcat[enzyme_data$organism_Name==i]~enzyme_data$log10KM[enzyme_data$organism_Name==i],pch=20,col=j)
        col_graph=c(col_graph,j)
        leg_graph<-c(leg_graph,i)
        j=j+1
}
legend("bottomright", legend=leg_graph,
       col=col_graph, pch=20, cex=0.5)
#conclusion: Coli and S.cerevisiae seem more eficient but bias very probable because of the
#enzyme module type repsrentation for each of them (see below); there may exist but mycobacterium and thermotoga would yet to be explained

#2.kcat vs KM/kcat
##a.Metabolic modules

#Two most significant categories opposing the trade-off between rate and affinity
toff_AA<-lm(enzyme_data$log10kcat[enzyme_data$module_Type=="AA,FA,N"]~enzyme_data$log10kcat_KM[enzyme_data$module_Type=="AA,FA,N"])
summary(toff_AA) #R^2=0.4854
toff_CEM<-lm(enzyme_data$log10kcat[enzyme_data$module_Type=="CEM"]~enzyme_data$log10kcat_KM[enzyme_data$module_Type=="CEM"])
summary(toff_CEM) #R^2=0.4208
toff_INTER<-lm(enzyme_data$log10kcat[enzyme_data$module_Type=="INTER"]~enzyme_data$log10kcat_KM[enzyme_data$module_Type=="INTER"])
toff_SEC2<-lm(enzyme_data$log10kcat[enzyme_data$module_Type=="SEC"]~enzyme_data$log10kcat_KM[enzyme_data$module_Type=="SEC"])
#In line with results for kcat versus KM: it is possible to overcome the rate-affinity trade-off

#Opening a plot with the appropriate scaling
par(mfrow=c(1,2),mai=c(0.95,0.95,0.95,0.95))
plot(enzyme_data$log10kcat~enzyme_data$log10kcat_KM,pch=20,col="white",
     xlab="log10(kcat/KM)",
     ylab="log10(kcat)",
     main="Simple linear regression",cex.main=1)

col_graph<-c()
leg_graph<-c()
j=1
for (i in levels(factor(enzyme_data$module_Type))){
        points(enzyme_data$log10kcat[enzyme_data$module_Type==i]~enzyme_data$log10kcat_KM[enzyme_data$module_Type==i],pch=20,col=j,cex=0.75)
        col_graph=c(col_graph,j)
        leg_graph<-c(leg_graph,i)
        j=j+1
}
legend("topleft", legend=leg_graph,
       col=col_graph, pch=20, cex=0.8)

abline(coef=c(toff_AA$coefficients[1],toff_AA$coefficients[2]),lwd=1.25,col=1)
abline(coef=c(toff_CEM$coefficients[1],toff_CEM$coefficients[2]),lwd=1.25,col=2)
abline(coef=c(toff_INTER$coefficients[1],toff_INTER$coefficients[2]),lwd=1.25,col=3)
abline(coef=c(toff_SEC2$coefficients[1],toff_SEC2$coefficients[2]),lwd=1.25,col=4)
text(0.75,-4,"A",srt=addtxt$srt,font=2,col="black",cex=1)
#Conclusion: confirmed the previously observed trend: modules seem on different isoclines, although
#less significantly that one might expect. Secondary module enzymes more prone to accept the rate-affinity trade-off

#Othogoanl regression
plot(enzyme_data$log10kcat~enzyme_data$log10kcat_KM,pch=20,col="white",
     xlab="log10(kcat/KM)",
     ylab="log10(kcat)",
     main="Orthogonal regression",cex.main=1)

col_graph<-c()
leg_graph<-c()
j=1
for (i in levels(factor(enzyme_data$module_Type))){
  points(enzyme_data$log10kcat[enzyme_data$module_Type==i]~enzyme_data$log10kcat_KM[enzyme_data$module_Type==i],pch=20,col=j,cex=0.75)
  col_graph=c(col_graph,j)
  leg_graph<-c(leg_graph,i)
  j=j+1
}
legend("topleft", legend=leg_graph,
       col=col_graph, pch=20, cex=0.8)
toff_ortho2AA <- deming(enzyme_data$log10kcat[enzyme_data$module_Type=="AA,FA,N"]~enzyme_data$log10kcat_KM[enzyme_data$module_Type=="AA,FA,N"])
#print(toff_ortho2AA)
abline(coef=c(toff_ortho2AA$coefficients[1],toff_ortho2AA$coefficients[2]),lty=2,lwd=2,col=1)
toff_ortho2CEM <- deming(enzyme_data$log10kcat[enzyme_data$module_Type=="CEM"]~enzyme_data$log10kcat_KM[enzyme_data$module_Type=="CEM"])
#print(toff_ortho2CEM)
abline(coef=c(toff_ortho2CEM$coefficients[1],toff_ortho2CEM$coefficients[2]),lty=2,lwd=2,col=2)
toff_ortho2CEM <- deming(enzyme_data$log10kcat[enzyme_data$module_Type=="INTER"]~enzyme_data$log10kcat_KM[enzyme_data$module_Type=="INTER"])
#print(toff_ortho2CEM)
abline(coef=c(toff_ortho2CEM$coefficients[1],toff_ortho2CEM$coefficients[2]),lty=2,lwd=2,col=3)
toff_ortho2SEC <- deming(enzyme_data$log10kcat[enzyme_data$module_Type=="SEC"]~enzyme_data$log10kcat_KM[enzyme_data$module_Type=="SEC"])
#print(toff_ortho2AA)
abline(coef=c(toff_ortho2SEC$coefficients[1],toff_ortho2SEC$coefficients[2]),lty=2,lwd=2,col=4)
#As previously observed, increases the trends, reality somewhere in between
legend("topleft", legend=leg_graph,
       col=col_graph, pch=20, cex=0.8)

text(0.75,-4,"B",srt=addtxt$srt,font=2,col="black",cex=1)
#legend("bottomright",title="Regression",legend=c("Simple linear","Deming"),lty=c(1,2),cex=1.25)

setwd(dir="~")
setwd(dir="Desktop/Ongoing-projects/")
dev.print(device = jpeg, file = "Bar-Even-spec-data.jpeg", , width = 1000*3,height=450*3,res=300,type="cairo")
dev.off()  

par(mfrow=c(1,2),mai=c(0.95,0.95,0.95,0.95))
plot(enzyme_data$log10kcat~enzyme_data$log10kcat_KM,pch=20,col="white",
     xlab="log10(kcat/KM)",
     ylab="log10(kcat)",
     main="",cex.main=1)

col_graph<-c()
leg_graph<-c()
j=1
for (i in c("CEM","AA,FA,N")){
  points(enzyme_data$log10kcat[enzyme_data$module_Type==i]~enzyme_data$log10kcat_KM[enzyme_data$module_Type==i],pch=20,col=j,cex=0.75)
  col_graph=c(col_graph,j)
  leg_graph<-c(leg_graph,i)
  j=j+1
}
legend("topleft", legend=leg_graph,
       col=col_graph, pch=20, cex=0.8)
text(0.75,-4,"C",srt=addtxt$srt,font=2,col="black",cex=1)
plot(enzyme_data$log10kcat~enzyme_data$log10kcat_KM,pch=20,col="white",
     xlab="log10(kcat/KM)",
     ylab="log10(kcat)",
     main="",cex.main=1)

col_graph<-c()
leg_graph<-c()
j=3
for (i in c("SEC","INTER")){
  points(enzyme_data$log10kcat[enzyme_data$module_Type==i]~enzyme_data$log10kcat_KM[enzyme_data$module_Type==i],pch=20,col=j,cex=0.75)
  col_graph=c(col_graph,j)
  leg_graph<-c(leg_graph,i)
  j=j+1
}
legend("topleft", legend=leg_graph,
       col=col_graph, pch=20, cex=0.8)

text(0.75,-4,"D",srt=addtxt$srt,font=2,col="black",cex=1)
#legend("bottomright",title="Regression",legend=c("Simple linear","Deming"),lty=c(1,2),cex=1.25)

setwd(dir="~")
setwd(dir="Desktop/Ongoing-projects/")
dev.print(device = jpeg, file = "Bar-Even-spec2-data.jpeg", , width = 1000*3,height=450*3,res=300,type="cairo")
dev.off()  








##b.Species

#Opening a plot with the appropriate scaling
plot(enzyme_data$log10kcat~enzyme_data$log10kcat_KM,pch=20,col="white",
     xlab="log10(kcat/KM)",
     ylab="log10(kcat)",
     main="Trade-off between kcat and kcat/KM")

j=1
col_graph<-c()
leg_graph<-c()
pch_graph<-c()
for (i in species_table){
        points(enzyme_data$log10kcat[enzyme_data$organism_Name==i]~enzyme_data$log10kcat_KM[enzyme_data$organism_Name==i],pch=j,col=j,cex=0.5)
        leg_graph<-c(leg_graph,i)
        col_graph=c(col_graph,j)
        pch_graph=c(pch_graph,j)
        j=j+1
}

legend("topleft", legend=leg_graph,
       col=col_graph, pch=pch_graph, cex=0.5)
##Inconclusive graph: no trend seems to be overwhelming

##Reducing the amount of species to cross species x module
species_table2<-c()

for (t in 1:length(table(enzyme_data$organism_Name))){
        if(table(enzyme_data$organism_Name)[t]>6)
                species_table2<-c(species_table2,names(table(enzyme_data$organism_Name)[t]))
}

#Opening a plot with the appropriate scaling
plot(enzyme_data$log10kcat~enzyme_data$log10kcat_KM,pch=20,col="white",
     xlab="log10(kcat/KM)",
     ylab="log10(kcat)",
     main="Trade-off between kcat and kcat/KM")
axis(tick = c())
pch_graph<-c()
col_graph<-c()
leg_graph<-c()
j=1
k=1
for (i in species_table2){
        k=1
        for (l in levels(factor(enzyme_data$module_Type))){
                print(paste(i,l))
                points(enzyme_data$log10kcat
                       [enzyme_data$organism_Name==i 
                               & enzyme_data$module_Type==l]
                       ~enzyme_data$log10kcat_KM[enzyme_data$organism_Name==i 
                                                 & enzyme_data$module_Type==l],
                       pch=k,col=j,cex=0.5)
                pch_graph=c(pch_graph,k)
                k=k+1
                leg_graph<-c(leg_graph,paste(i,l))
        }
        col_graph=c(col_graph,rep(j,length(levels(factor(enzyme_data$module_Type)))))
        j=j+1
}
legend("topleft",ncol=round(j/2), legend=leg_graph,
       col=col_graph, pch=pch_graph, cex=0.3)
#Conclusion: enzymes involved in very different modules, absolutely not representative
#No possibility to go further with more complex linear models, data too patchy and unevenly distributed

#Crossing species x "CEM"
#Opening a plot with the appropriate scaling
species_table2_bis<-c()

for (t in 1:length(table(enzyme_data$organism_Name))){
  if(table(enzyme_data$organism_Name)[t]>2)
    species_table2_bis<-c(species_table2_bis,names(table(enzyme_data$organism_Name)[t]))
  #sub_enz_data<-cbind(sub_enz_data,subset(enzyme_data,enzyme_data$organism_Name=="Homo sapiens"))
}

plot(enzyme_data$log10kcat~enzyme_data$log10kcat_KM,pch=20,col="white",
     xlab="log10(kcat/KM)",
     ylab="log10(kcat)",
     main="Trade-off between kcat and kcat/KM")
pch_graph<-c()
col_graph<-c()
leg_graph<-c()
j=1
k=1
col.pal<-rainbow(16)
for (i in species_table2_bis){
  #k=19
  print(j)
  print(i)
  for (l in c("CEM")){#levels(factor(enzyme_data$module_Type))){
    if(length(enzyme_data$log10kcat
              [enzyme_data$organism_Name==i 
                & enzyme_data$module_Type==l])>0){
                  points(enzyme_data$log10kcat
                         [enzyme_data$organism_Name==i 
                           & enzyme_data$module_Type==l]
                         ~enzyme_data$log10kcat_KM[enzyme_data$organism_Name==i 
                                                   & enzyme_data$module_Type==l],
                         pch=k,col=col.pal[j],cex=0.5)
                  print(enzyme_data$log10kcat_KM[enzyme_data$organism_Name==i 
                                                 & enzyme_data$module_Type==l])
                  pch_graph=c(pch_graph,k)
                  k=k+1
                  leg_graph<-c(leg_graph,paste(i,l))  
                  col_graph=c(col_graph,col.pal[j])
                  j=j+1
        }
    
  }
  
}
legend("topleft",ncol=round(j/4), legend=leg_graph,
       col=col_graph, pch=pch_graph, cex=0.5)


##Crossing eukaryotes/prokaryotes x module
species_table3<-c()

for (t in 1:length(table(enzyme_data$organism_Name))){
        if(table(enzyme_data$organism_Name)[t]>2)
                species_table3<-c(species_table3,names(table(enzyme_data$organism_Name)[t]))
                #sub_enz_data<-cbind(sub_enz_data,subset(enzyme_data,enzyme_data$organism_Name=="Homo sapiens"))
}

empty_col_str<-c()

for (i in 1:nrow(enzyme_data)){
        empty_col_str=c(empty_col_str,"")
}

enzyme_data<-cbind(enzyme_data,empty_col_str)
colnames(enzyme_data)[colnames(enzyme_data) == 'empty_col_str'] <- 'Organism_category'

#Manual assigning for >4 enzymes
#species_cat<-c("Prok","Euk","Prok","Euk","Prok",
               #"Euk","Prok","Euk","Euk","Prok",
               #"Euk","Euk","Prok","Prok","Prok")

#Manual assigning for >3 enzymes
#species_cat<-c("Prok","Euk","Prok","Euk","Prok",
 #              "Euk","Prok","Euk","Euk","Prok","Prok",
  #             "Prok","Euk","Euk","Prok","Prok","Prok","Euk","Prok")

#Manual assigning for >2 enzymes
species_cat<-c("Prok","Prok","Euk","Euk","Prok","Euk","Prok","Prok",
               "Euk","Prok","Euk","Euk","Prok","Euk","Prok","Prok",
               "Euk","Euk","Prok","Prok","Prok","Euk","Euk","Prok",
               "Prok","Prok","Euk","Prok")

species_table3<-data.frame(cbind(species_table3,species_cat))
enzyme_data$Organism_category<-as.character(enzyme_data$Organism_category)

for (i in 1:nrow(enzyme_data)){
        for (j in 1:nrow(species_table3)){
                if (enzyme_data[i,]$organism_Name==species_table3[j,1]){
                        print(TRUE)
                        enzyme_data[i,]$Organism_category=as.character(species_table3[j,2])
                }
        }
        
}

#Opening a plot with the appropriate scaling
plot(enzyme_data$log10kcat~enzyme_data$log10kcat_KM,pch=20,col="white",
     xlab="log10(kcat/KM)",
     ylab="log10(kcat)",
     main="Trade-off between kcat and kcat/KM")
pch_graph<-c()
col_graph<-c()
leg_graph<-c()
j=1
k=1
for (i in levels(species_table3$species_cat)){
        k=1
        for (l in levels(factor(enzyme_data$module_Type))){
                print(paste(i,l))
                points(enzyme_data$log10kcat
                       [enzyme_data$Organism_category==i 
                               & enzyme_data$module_Type==l]
                       ~enzyme_data$log10kcat_KM[enzyme_data$Organism_category==i 
                                                 & enzyme_data$module_Type==l],
                       pch=k,col=j,cex=0.5)
                pch_graph=c(pch_graph,k)
                k=k+1
                leg_graph<-c(leg_graph,paste(i,l))
        }
        col_graph=c(col_graph,rep(j,length(levels(factor(enzyme_data$module_Type)))))
        j=j+1
}
legend("topleft", legend=leg_graph,
       col=col_graph, pch=pch_graph, cex=0.5)

#Two most significant categories opposing the trade-off between rate and affinity
toff_AA_Euk<-lm(enzyme_data$log10kcat[enzyme_data$module_Type=="AA,FA,N" & enzyme_data$Organism_category=="Euk"]~enzyme_data$log10kcat_KM[enzyme_data$module_Type=="AA,FA,N" & enzyme_data$Organism_category=="Euk"])
summary(toff_AA_Euk) #R^2=0.4854
abline(coef=c(toff_AA_Euk$coefficients[1],toff_AA_Euk$coefficients[2]),lwd=2,col=1)
toff_CEM_Euk<-lm(enzyme_data$log10kcat[enzyme_data$module_Type=="CEM" & enzyme_data$Organism_category=="Euk"]~enzyme_data$log10kcat_KM[enzyme_data$module_Type=="CEM" & enzyme_data$Organism_category=="Euk"])
abline(coef=c(toff_CEM_Euk$coefficients[1],toff_CEM_Euk$coefficients[2]),lwd=2,col=2)
summary(toff_CEM_Euk) #R^2=0.4208

toff_AA_Euk_ortho<- deming(enzyme_data$log10kcat[enzyme_data$module_Type=="AA,FA,N" & enzyme_data$Organism_category=="Euk"]~enzyme_data$log10kcat_KM[enzyme_data$module_Type=="AA,FA,N" & enzyme_data$Organism_category=="Euk"])
abline(coef=c(toff_AA_Euk_ortho$coefficients[1],toff_AA_Euk_ortho$coefficients[2]),lwd=1,lty=2,col=1)
toff_CEM_Euk_ortho<- deming(enzyme_data$log10kcat[enzyme_data$module_Type=="CEM" & enzyme_data$Organism_category=="Euk"]~enzyme_data$log10kcat_KM[enzyme_data$module_Type=="CEM" & enzyme_data$Organism_category=="Euk"])
abline(coef=c(toff_CEM_Euk_ortho$coefficients[1],toff_CEM_Euk_ortho$coefficients[2]),lwd=1,lty=2,col=2)

toff_AA_Prok<-lm(enzyme_data$log10kcat[enzyme_data$module_Type=="AA,FA,N" & enzyme_data$Organism_category=="Prok"]~enzyme_data$log10kcat_KM[enzyme_data$module_Type=="AA,FA,N" & enzyme_data$Organism_category=="Prok"])
summary(toff_AA_Prok) #R^2=0.4854
abline(coef=c(toff_AA_Prok$coefficients[1],toff_AA_Prok$coefficients[2]),lwd=2,col=3)
toff_CEM_Prok<-lm(enzyme_data$log10kcat[enzyme_data$module_Type=="CEM" & enzyme_data$Organism_category=="Prok"]~enzyme_data$log10kcat_KM[enzyme_data$module_Type=="CEM" & enzyme_data$Organism_category=="Prok"])
summary(toff_CEM_Prok)
abline(coef=c(toff_CEM_Prok$coefficients[1],toff_CEM_Prok$coefficients[2]),lwd=2,col=4)

toff_AA_Prok_ortho<- deming(enzyme_data$log10kcat[enzyme_data$module_Type=="AA,FA,N" & enzyme_data$Organism_category=="Prok"]~enzyme_data$log10kcat_KM[enzyme_data$module_Type=="AA,FA,N" & enzyme_data$Organism_category=="Prok"])
abline(coef=c(toff_AA_Prok_ortho$coefficients[1],toff_AA_Prok_ortho$coefficients[2]),lwd=1,lty=2,col=3)
toff_CEM_Prok_ortho<- deming(enzyme_data$log10kcat[enzyme_data$module_Type=="CEM" & enzyme_data$Organism_category=="Prok"]~enzyme_data$log10kcat_KM[enzyme_data$module_Type=="CEM" & enzyme_data$Organism_category=="Prok"])
abline(coef=c(toff_CEM_Prok_ortho$coefficients[1],toff_CEM_Prok_ortho$coefficients[2]),lwd=1,lty=2,col=4)
