N_reso=250#250
KT_set=c(10^-5,10^-1) #Molar, Michaelis-like saturation constant for facilitated diffusion
Vm=c(10^-6,10^-3) #unit: M.s^-1, with Vm=5.10^-5pmol.s^-1.cell^-1 and Vcell=50 micrometers^3

kr_set=10^(2.8)
E_tot_conc=10^-3
##Defining variables and set of parameters
log10kcatKM_var<-c(0,10,2)
kcatKM<-10^seq(log10kcatKM_var[1],log10kcatKM_var[2],length.out=N_reso)
log10kcat_var<-c(-4,6,2)
kcat<-10^seq(log10kcat_var[1],log10kcat_var[2],length.out=N_reso)

tab_P_Phi_eq <- matrix(nrow=N_reso , ncol=N_reso)
P_Phi_eq<-list()

for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){
        KT=KT_set[s]
        Se=10*KT
        kr<-kr_set[l]
        kf=kcatKM[i]*(1+kr/kcat[j])
        if(kf/10^10>1){
          tab_P_Phi_eq[i,j]=-10^-20
        }
        else{
          a<-mpfr((Vm[s]*kf+kf*E_tot_conc*kcat[j]*(1+Se/KT)),80)
          b<- mpfr((Vm[s]*(kr+kcat[j]-Se*kf)+kf*kcat[j]*E_tot_conc*(KT+Se)),80)
          c<- mpfr(-Vm[s]*Se*(kr+kcat[j]),80)
          delta<-mpfr(b^2-4*a*c,80)
          S_conc_eq=mpfr((-b+delta^(1/2))/(2*a),80)
          ES_conc_eq=mpfr(E_tot_conc*(kf*S_conc_eq)/(kr+kcat[j]+kf*S_conc_eq),80)
          tab_P_Phi_eq[i,j]=as.numeric(kcat[j]*ES_conc_eq)
        }
      }
    }
    P_Phi_eq[[s+length(Vm)*(l-1)]]<-tab_P_Phi_eq
  }
}

##Calculating fitness accounting for non-accessibility of some values (above diffusion limit)
w_kcatKM_tab<- matrix(nrow=N_reso , ncol=N_reso)
w_kcatKM<-list()
for (l in 1:(length(kr_set))){
  for (s in 1:(length(Vm))){
    for (i in 1:N_reso){
      for (j in 1:N_reso){ 
        if(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]>0){
          w_kcatKM_tab[i,j]<-round(P_Phi_eq[[s+length(Vm)*(l-1)]][i,j]/max(P_Phi_eq[[s+length(Vm)*(l-1)]]),digits=10)
        }
        else{
          w_kcatKM_tab[i,j]<-1.1 #For values above diffusion limit
        }
      }
    }
    w_kcatKM[[s+length(Vm)*(l-1)]]<-w_kcatKM_tab
  }
}
#Defining features of plots
ncol=128
jet.colors <- colorRampPalette(c("white")) #palette for fitness
palet<-jet.colors(ncol)
pal<-list(palet)

#Low flux
addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)
addtxt2<-list(l=0.05,h=0.95,txt=c("A","B"),srt = 0,font=2,col="black",cex=1)


##Plot with nomutationl bias
f="W"
b=mut_bias_set[2]
i=1

for (n in Ne_set){
  Evo_data[[i]]<-as.data.frame(cbind(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))))
  i=i+1
}
Evo_data_full_low_b<-rbind(Evo_data[[1]],Evo_data[[2]],Evo_data[[3]],Evo_data[[4]])
##Correlation test for simulation - between kcat and kcat_KM
cor.test(Evo_data_full_low_b$V1,Evo_data_full_low_b$V2)

addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)
legtext1=c(2,3,4,5)

multiplePlot("kr=","","Vm=","M/s",KT,Vm[1],ncol=128,log10kcatKM_var,log10kcat_var,w_kcatKM,
             abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),
             xyData=Evo_data,TEXT_to_Add =addtxt1,palette=pal,image=TRUE,ncat=4,pcex=0.7,subcex=1,labcex=1.8,axcex=0.9,cextext=2,legtitle="log10(Ne)",col_data = c(1,2,3,4),pch_dat=c(rep(4,4)),colorkey=FALSE,globcex=0.5,legcex=1)

levels=c()
for(cont in 1:4){
  levels=c(levels,1-10^(-cont-1))
}
contour(w_kcatKM[[1]],xaxt="n",yaxt="n",levels=levels,labcex=1,col=c(1,2,3,4),drawlabels = FALSE,lty=6,lwd=0.5,xlim=c(0,1),add=TRUE)
abline(a=-0.32,b=1,col="lightgray",lty=1,lwd=5)
text(0.375,0.959,"0.99",srt=-90,font=1,col=1,cex=1.65,family="serif")
#text(0.375,0.2,"0.99",srt=-90,font=1,col=1,cex=1.75,family="serif")
text(0.475,0.948,"0.999",srt=-90,font=1,col=2,cex=1.65,family="serif")
text(0.58,0.937,"0.9999",srt=-90,font=1,col=3,cex=1.65,family="serif")
text(0.78,0.39,"0.99999",srt=-0,font=1,col=4,cex=1.65,family="serif")
text(addtxt2$l,addtxt2$h,addtxt2$txt[1],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=2)
setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
dev.print(device = jpeg, file = "2DFitLandscape_Evo_Results_lowF_nobias.jpeg", width = 575*3,height=460*3,res=100*3,type="cairo")

##Plot with low mutationl bias
f="L"
b=mut_bias_set[2]
i=1

for (n in Ne_set){
  Evo_data[[i]]<-as.data.frame(cbind(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))))
  i=i+1
}
Evo_data_full_low_b<-rbind(Evo_data[[1]],Evo_data[[2]],Evo_data[[3]],Evo_data[[4]])
##Correlation test for simulation - between kcat and kcat_KM
cor.test(Evo_data_full_low_b$V1,Evo_data_full_low_b$V2)

addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)
legtext1=c(2,3,4,5)

multiplePlot("kr=","","Vm=","M/s",KT,Vm[1],ncol=128,log10kcatKM_var,log10kcat_var,w_kcatKM,
             abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),
             xyData=Evo_data,TEXT_to_Add =addtxt1,palette=pal,image=TRUE,ncat=4,pcex=0.6,subcex=1,labcex=1.8,axcex=0.9,legtext = legtext1,legpos=c(1.125,0.5),cextext=2,legtitle="log10(Ne)",col_data = c(1,2,3,4),pch_dat=c(rep(19,4)),colorkey=FALSE,legcex=3,globcex=0.5)

levels=c()
for(cont in 1:4){
  levels=c(levels,1-10^(-cont-1))
}
contour(w_kcatKM[[2]],xaxt="n",yaxt="n",levels=levels,labcex=1,col=c(1,2,3,4),drawlabels = FALSE,lty=6,lwd=0.5,xlim=c(0,1),add=TRUE)
abline(a=-0.32,b=1,col="lightgray",lty=1,lwd=5)
text(0.325,0.959,"0.99",srt=-90,font=1,col=1,cex=1.75,family="serif")
text(0.425,0.948,"0.999",srt=-90,font=1,col=2,cex=1.75,family="serif")
text(0.525,0.937,"0.9999",srt=-90,font=1,col=3,cex=1.75,family="serif")
text(0.625,0.926,"0.99999",srt=-90,font=1,col=4,cex=1.75,family="serif")
legend(x=1.05,y=0.9,title="Mutational bias",legend=c("b = 0","b = - 0.1","b = - 0.2"),pch=c(4,19,17),col=1,xpd=NA,cex=1.8)

#Low mutational bias
##a.Importing simulation results
setwd(dir="~")
setwd(dir="Desktop/enzyme-evolution/Evo_Results")
flux_set=c("W")#,"H")
Ne_set=c(10^2,10^3,10^4,10^5)
mut_bias_set<-c(-0.1,-0.2)
for (f in flux_set){
  for(n in Ne_set){
    for(b in mut_bias_set){
      mu<-10^-1/n
      kr<-10^3
      print(paste(f,n,mu,b))
      load(paste("/Users/florian/Desktop/enzyme-evolution/Evo_Results/","eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,".rda",sep=""))
      if(b==0){
        assign(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi00",sep=""),evo_results)
      }
      else{
        assign(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""),evo_results)
      }
      tempo_kcat<-c()
      tempo_kf<-c()
      tempo_Phi<-c()
      #tempo_time_Fit<-list()
      for (e in 1:30){
        long=length(get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]])
        tempo_kcat<-c(tempo_kcat,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][1])
        tempo_kf<-c(tempo_kf,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][2+(long/3-1)])
        tempo_Phi<-c(tempo_Phi,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][3+2*(long/3-1)])
      }
      tempo_time_Fit<-c()
      g_set<-c()
      for (g in 1:10){
        for (e in 1:30){
          long=length(get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[g]][[e]])
          tempo_time_Fit[(g-1)*30+e]<-get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[g]][[e]][3+2*(long/3-1)]
          g_set<-c(g_set,g)
        }
      }
      assign(paste("P_Phi_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_Phi)
      assign(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_kcat)
      assign(paste("kf_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_kf)
      assign(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_kf*tempo_kcat/(tempo_kcat+kr))
      assign(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""),as.data.frame(cbind(g_set,tempo_time_Fit)))
    }
  }
}

N_rep=30
param_set<-data.frame("Ne_set"=rep(Ne_set,length(mut_bias_set)),"bias"=rep(mut_bias_set,each=length(Ne_set)))

#Getting the fitness as as feature of simulations
flux_set=c("W")#,"H")
fit_evo_eq_summary<-list()
Ne_set_data<-list()
bias_set_data<-list()
for (f in flux_set){
  fit_evo_eq_summary[[f]]<-c()
  Ne_set_data[[f]]<-c()
  bias_set_data[[f]]<-c()
  for (s in 1:8){
    n=param_set$Ne_set[s]
    b=param_set$bias[s]
    Ne_set_data[[f]]<-c(Ne_set_data[[f]],rep(n,30))
    bias_set_data[[f]]<-c(bias_set_data[[f]],rep(b,30))
    fit_evo_eq_summary[[f]]<-c(fit_evo_eq_summary[[f]],get(paste("P_Phi_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))
  }
}

data_fit_eq<-list()
col_headings <- c("Ne","bias","fitness")
for (f in flux_set){
  data_fit_eq[[f]]<-data.frame(cbind(Ne_set_data[[f]],bias_set_data[[f]],fit_evo_eq_summary[[f]]))
  names(data_fit_eq[[f]]) <- col_headings
}
b=mut_bias_set[2]
i=1
for (n in Ne_set){
  Evo_data[[i]]<-as.data.frame(cbind(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))))
  i=i+1
}

##Correlation test for simulation - between kcat and kcat_KM
Evo_data_full_high_b<-rbind(Evo_data[[1]],Evo_data[[2]],Evo_data[[3]],Evo_data[[4]])
cor.test(Evo_data_full_high_b$V1,Evo_data_full_high_b$V2)
Evo_data_full_lowf<-rbind(Evo_data_full_low_b,Evo_data_full_high_b)
cor.test(Evo_data_full_lowf$V1,Evo_data_full_lowf$V2)

xyData=Evo_data
xlim=c(0,10)
ylim=c(-4,6)
for (i in 1:length(xyData)){
  points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
         (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
         col=i,pch=17,cex=0.8)
}

text(addtxt2$l,addtxt2$h,addtxt2$txt[2],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=2)

dev.print(device = jpeg, file = "2DFitLandscape_Evo_Results_highF_withbias.jpeg", width = 575*3,height=460*3,res=100*3,type="cairo")




##Plot with low mutationl bias
mut_bias_set<-c(-0.1,-0.2)
f="W"
b=mut_bias_set[1]
i=1

for (n in Ne_set){
  Evo_data[[i]]<-as.data.frame(cbind(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))))
  i=i+1
}
Evo_data_full_low_b<-rbind(Evo_data[[1]],Evo_data[[2]],Evo_data[[3]],Evo_data[[4]])
##Correlation test for simulation - between kcat and kcat_KM
cor.test(Evo_data_full_low_b$V1,Evo_data_full_low_b$V2)

addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)
legtext1=c(2,3,4,5)

multiplePlot("","","","","","",ncol=128,log10kcatKM_var,log10kcat_var,w_kcatKM,
             abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),
             xyData=Evo_data,TEXT_to_Add =addtxt1,palette=pal,image=TRUE,ncat=4,pcex=0.7,subcex=1,labcex=1.8,axcex=0.9,legtext = legtext1,legpos=c(1.125,0.5),
             cextext=2,legtitle="log10(Ne)",col_data = c(1,2,3,4),pch_dat=c(rep(19,4)),colorkey=FALSE,legcex=2.5,globcex=0.5)

levels=c()
for(cont in 1:4){
  levels=c(levels,1-10^(-cont-1))
}
contour(w_kcatKM[[1]],xaxt="n",yaxt="n",levels=levels,labcex=1,col=c(1,2,3,4),drawlabels = FALSE,lty=6,lwd=0.5,xlim=c(0,1),add=TRUE)
abline(a=-0.32,b=1,col="lightgray",lty=1,lwd=5)
text(0.375,0.959,"0.99",srt=-90,font=1,col=1,cex=1.75,family="serif")
text(0.475,0.948,"0.999",srt=-90,font=1,col=2,cex=1.75,family="serif")
text(0.575,0.937,"0.9999",srt=-90,font=1,col=3,cex=1.75,family="serif")
text(0.675,0.926,"0.99999",srt=-90,font=1,col=4,cex=1.75,family="serif")
legend(x=1.08,y=0.9,title="Mutational bias",legend=c("b = - 0.1","b = - 0.2"),pch=c(19,17),col=1,xpd=NA,cex=1.6)
#legend(x=1.115,y=0.6,title="Protein cost",legend=c(expression(paste(10^{-3})),expression(paste(10^{-4}))),col=1,pch=c(6,2),xpd=TRUE,cex=1.6)
text(addtxt2$l,addtxt2$h,addtxt2$txt[1],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
axis(side=1,labels=FALSE,tck=FALSE)
axis(side=4,labels=FALSE,tck=FALSE)
#High mutational bias
b=mut_bias_set[2]
i=1
for (n in Ne_set){
  Evo_data[[i]]<-as.data.frame(cbind(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))))
  i=i+1
}

##Correlation test for simulation - between kcat and kcat_KM
Evo_data_full_high_b<-rbind(Evo_data[[1]],Evo_data[[2]],Evo_data[[3]],Evo_data[[4]])
#cor.test(Evo_data_full_high_b$V1,Evo_data_full_high_b$V2)
Evo_data_full_lowf<-rbind(Evo_data_full_low_b,Evo_data_full_high_b)
#cor.test(Evo_data_full_lowf$V1,Evo_data_full_lowf$V2)
xyData=Evo_data
xlim=c(0,10)
ylim=c(-4,6)
for (i in 1:length(xyData)){
  points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
         (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
         col=i,pch=17,cex=0.9)
}

#text(addtxt2$l,addtxt2$h,addtxt2$txt[2],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=2)
setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
dev.print(device = jpeg, file = "2DFitLandscape_Evo_Results_lowF_withbias_Def.jpeg", width = 575*3,height=460*3,res=100*3,type="cairo")




#Low mutational bias
##a.Importing simulation results
setwd(dir="~")
setwd(dir="Desktop/enzyme-evolution/Evo_Results")
flux_set=c("H")#,"H")
Ne_set=c(10^2,10^3,10^4,10^5)
mut_bias_set<-c(-0.1,-0.2)
for (f in flux_set){
  for(n in Ne_set){
    for(b in mut_bias_set){
      mu<-10^-1/n
      kr<-10^3
      print(paste(f,n,mu,b))
      load(paste("/Users/florian/Desktop/enzyme-evolution/Evo_Results/","eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,".rda",sep=""))
      if(b==0){
        assign(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi00",sep=""),evo_results)
      }
      else{
        assign(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""),evo_results)
      }
      tempo_kcat<-c()
      tempo_kf<-c()
      tempo_Phi<-c()
      #tempo_time_Fit<-list()
      for (e in 1:30){
        long=length(get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]])
        tempo_kcat<-c(tempo_kcat,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][1])
        tempo_kf<-c(tempo_kf,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][2+(long/3-1)])
        tempo_Phi<-c(tempo_Phi,get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[10]][[e]][3+2*(long/3-1)])
      }
      tempo_time_Fit<-c()
      g_set<-c()
      for (g in 1:10){
        for (e in 1:30){
          long=length(get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[g]][[e]])
          tempo_time_Fit[(g-1)*30+e]<-get(paste("eq_",f,"_Ne",n,"_mu",-log10(mu),"_bi0",-b*10,sep=""))[[g]][[e]][3+2*(long/3-1)]
          g_set<-c(g_set,g)
        }
      }
      assign(paste("P_Phi_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_Phi)
      assign(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_kcat)
      assign(paste("kf_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_kf)
      assign(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""),tempo_kf*tempo_kcat/(tempo_kcat+kr))
      assign(paste("Fit_evo_time",f,"_Ne",n,"_bi0",-b*10,sep=""),as.data.frame(cbind(g_set,tempo_time_Fit)))
    }
  }
}

N_rep=30
param_set<-data.frame("Ne_set"=rep(Ne_set,length(mut_bias_set)),"bias"=rep(mut_bias_set,each=length(Ne_set)))

#Getting the fitness as as feature of simulations
flux_set=c("W")#,"H")
fit_evo_eq_summary<-list()
Ne_set_data<-list()
bias_set_data<-list()
for (f in flux_set){
  fit_evo_eq_summary[[f]]<-c()
  Ne_set_data[[f]]<-c()
  bias_set_data[[f]]<-c()
  for (s in 1:8){
    n=param_set$Ne_set[s]
    b=param_set$bias[s]
    Ne_set_data[[f]]<-c(Ne_set_data[[f]],rep(n,30))
    bias_set_data[[f]]<-c(bias_set_data[[f]],rep(b,30))
    fit_evo_eq_summary[[f]]<-c(fit_evo_eq_summary[[f]],get(paste("P_Phi_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))
  }
}



##Plot with low mutationl bias
mut_bias_set<-c(-0.1,-0.2)
f="H"
b=mut_bias_set[1]
i=1

for (n in Ne_set){
  Evo_data[[i]]<-as.data.frame(cbind(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))))
  i=i+1
}
Evo_data_full_low_b<-rbind(Evo_data[[1]],Evo_data[[2]],Evo_data[[3]],Evo_data[[4]])
##Correlation test for simulation - between kcat and kcat_KM
cor.test(Evo_data_full_low_b$V1,Evo_data_full_low_b$V2)

addtxt1<-list(l=0.75,h=0.2,txt="Diffusion limit",srt = 45,font=1,col="black",cex=2)
addtxt12<-list(l=0.8,h=0.18,txt=expression(paste("(kr=",10^3,s^{-1},")")),srt = 45,font=3,col="black",cex=1)
legtext1=c(2,3,4,5)

multiplePlot("","","","","","",ncol=128,log10kcatKM_var,log10kcat_var,w_kcatKM,
             abs="log10 kcat/KM",ord="log10 kcat",scale=c(0,1),
             xyData=Evo_data,TEXT_to_Add =addtxt1,palette=pal,image=TRUE,ncat=4,pcex=0.7,subcex=1,labcex=1.8,axcex=0.9,legtext = legtext1,legpos=c(1.125,0.5),
             cextext=2,legtitle="log10(Ne)",col_data = c(1,2,3,4),pch_dat=c(rep(19,4)),colorkey=FALSE,legcex=2.5,globcex=0.5)

levels=c()
for(cont in 1:4){
  levels=c(levels,1-10^(-cont-1))
}
contour(w_kcatKM[[2]],xaxt="n",yaxt="n",levels=levels,labcex=1,col=c(1,2,3,4),drawlabels = FALSE,lty=6,lwd=0.5,xlim=c(0,1),add=TRUE)
abline(a=-0.32,b=1,col="lightgray",lty=1,lwd=5)
text(0.275,0.959,"0.99",srt=-90,font=1,col=1,cex=1.75,family="serif")
text(0.375,0.948,"0.999",srt=-90,font=1,col=2,cex=1.75,family="serif")
text(0.475,0.937,"0.9999",srt=-90,font=1,col=3,cex=1.75,family="serif")
text(0.575,0.926,"0.99999",srt=-90,font=1,col=4,cex=1.75,family="serif")
legend(x=1.08,y=0.9,title="Mutational bias",legend=c("b = - 0.1","b = - 0.2"),pch=c(19,17),col=1,xpd=NA,cex=1.6)
#legend(x=1.115,y=0.6,title="Protein cost",legend=c(expression(paste(10^{-3})),expression(paste(10^{-4}))),col=1,pch=c(6,2),xpd=TRUE,cex=1.6)
text(addtxt2$l,addtxt2$h,addtxt2$txt[2],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
axis(side=1,labels=FALSE,tck=FALSE)
axis(side=4,labels=FALSE,tck=FALSE)
#High mutational bias
b=mut_bias_set[2]
i=1
for (n in Ne_set){
  Evo_data[[i]]<-as.data.frame(cbind(log10(get(paste("kcat_KM_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep=""))),log10(get(paste("kcat_ev_eq_fl",f,"_Ne",n,"_bi0",-b*10,sep="")))))
  i=i+1
}

##Correlation test for simulation - between kcat and kcat_KM
Evo_data_full_high_b<-rbind(Evo_data[[1]],Evo_data[[2]],Evo_data[[3]],Evo_data[[4]])
#cor.test(Evo_data_full_high_b$V1,Evo_data_full_high_b$V2)
Evo_data_full_lowf<-rbind(Evo_data_full_low_b,Evo_data_full_high_b)
#cor.test(Evo_data_full_lowf$V1,Evo_data_full_lowf$V2)
xyData=Evo_data
xlim=c(0,10)
ylim=c(-4,6)
for (i in 1:length(xyData)){
  points((xyData[[i]][[1]]-(xlim[1]))/((xlim[2]-xlim[1])),
         (xyData[[i]][[2]]-(ylim[1]))/(ylim[2]-ylim[1]),
         col=i,pch=17,cex=0.9)
}

#text(addtxt2$l,addtxt2$h,addtxt2$txt[2],srt=addtxt2$srt,font=addtxt2$font,col=addtxt2$col,cex=2)
text(addtxt12$l,addtxt12$h,addtxt12$txt[1],srt=addtxt12$srt,font=addtxt12$font,col=addtxt12$col,cex=2)
setwd(dir="~")
setwd(dir="enzyme-evolution/Paper/Figures")
dev.print(device = jpeg, file = "2DFitLandscape_Evo_Results_highF_withbias_Def.jpeg", width = 575*3,height=460*3,res=100*3,type="cairo")
