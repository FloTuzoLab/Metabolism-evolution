simul=15
flux="H"
deterministic=TRUE

setwd(dir="~")
#setwd(dir="Desktop/ScriptsR_modeles/")
#source("AP1.Plotting_multiple_images.R") #To plot multiple gradient and color plots
#source("AP2.RaphsonNewton.R") #To find equilibrium for 2 variables problems
options(digits=15)
setwd(dir="Desktop/enzyme-evolution/Evo_Results")

Ne_set=c(10^2,10^3,10^4,10^5)
N_rep=30
mut_bias_set<-c(-0.2)#<-c(-0.1)
cost_set<-c(10^-5,10^-4,10^-3,10^-2)
param_set<-data.frame("Ne_set"=rep(Ne_set,length(cost_set)),"cost_set"=rep(cost_set,each=length(Ne_set)))

kcat_i<-10^-1
kf_i<-10^3
kr_i<-10^3
E_tot_conc_i=10^-5
N_enz_init=round(10^-5/(1.6*10^-9))

KT_set=c(10^-5,10^-1)
Vm_set=c(10^-6,10^-3)

if(flux=="W"){
  KT=KT_set[1]
  Vm=Vm_set[1]
}
if(flux=="H"){
  KT=KT_set[2]
  Vm=Vm_set[2]
}

Se=10*KT

K<-2.5*10^-3#corresponding to half fitness concentration for which proteins correspond to half the cost of metabolic budget
#Flux determination

  
  fit<-function(kf,kr,kcat,E_conc_m,cost_Enz,N_ind){
    ##If deterministic relationship
    if(deterministic){
      E_conc=E_conc_m
      Etot_back<-2.5*10^-3#broadly corresponds to E.coli or S.cerevisae
      #Net_Flux=0
      kf_act<-kf*10^(-(E_conc+Etot_back)/(5*10^-3))
      a<- (Vm*kf_act+kf_act*E_conc*kcat*(1+Se/KT))
      b<- (Vm*(kr+kcat-Se*kf_act)+kf_act*kcat*E_conc*(KT+Se))
      c<- -Vm*Se*(kr+kcat)
      delta<-b^2-4*a*c
      S_conc_eq=(-b+delta^(1/2))/(2*a)
      ES_conc_eq=E_conc*(kf_act*S_conc_eq)/(kr+kcat+kf_act*S_conc_eq)
      Flux=kcat*ES_conc_eq
        #print(Flux)
        #Net_Flux=Flux*K/(E_conc+2*K)
      if(Flux-E_conc*cost_Enz>0){
        Net_Flux=Flux-E_conc*cost_Enz
      }
      else{
        Net_Flux=10^-20
      }
    }
    ##If noisiness introduced
    else{
      Etot_back<-2.5*10^-3#broadly corresponds to E.coli or S.cerevisae
      Net_Flux=0
      if(N_ind!=0){
        for(n in 1:N_ind){
          E_conc<-rpois(1,E_conc_m/(100*1.6*10^-9))*(100*1.6*10^-9)
          kf_act<-kf*10^(-(E_conc+Etot_back)/(5*10^-3))
          a<- (Vm*kf_act+kf_act*E_conc*kcat*(1+Se/KT))
          b<- (Vm*(kr+kcat-Se*kf_act)+kf_act*kcat*E_conc*(KT+Se))
          c<- -Vm*Se*(kr+kcat)
          delta<-b^2-4*a*c
          S_conc_eq=(-b+delta^(1/2))/(2*a)
          ES_conc_eq=E_conc*(kf_act*S_conc_eq)/(kr+kcat+kf_act*S_conc_eq)
          Flux=kcat*ES_conc_eq
          #print(Flux)
          #Net_Flux=Flux*K/(E_conc+2*K)
          if(Flux-E_conc*cost_Enz>0){
            Net_Flux=Net_Flux+Flux-E_conc*cost_Enz
          }
          else{
            Net_Flux=Net_Flux+10^-20
          }
        }
        Net_Flux=Net_Flux/N_ind
      }
    return(Net_Flux)
  }
}

#Max_val<-fit(10^8,10^3,10^4,10^-6,10^-3)
evo_gen_results<-list()
evo_results<-list()

#G<-5*10^2*Ne_set[1]#used for mutational bias
#for(l in 1:length(mut_bias_set)){
 # evo_results[[l]]<-list()
  #mut_bias<-mut_bias_set[l]
  #for (n in 1:length(Ne_set)){
#for (g in 1:length(G_set)){
mut_bias=mut_bias_set[1]
  #if(mut_bias=="increasing"){
    #mut_bias_coeff=-0.05
  #}else{
    #mut_bias<-as.numeric(mut_bias)
  #}
  cost<-param_set$cost_set[simul]
  Ne<-param_set$Ne_set[simul]
  G<-10^3*Ne##2.5*10^2#In principle 10^3
  mu<-10^-1/Ne
  nam<-paste("eq_",flux,"_Ne",Ne,"_mu",-log10(mu),"_bi0",-mut_bias*10,sep="")
  #G<-G_set[g]
  #evo_results[[g]]<-list()
  for (g in 1:10){
    evo_results[[g]]<-list()
  }
  
    for(r in 1:N_rep){
      print(r)
      g=1
      #print(l)
    #print(r)
      
    #List of (Ne) cells
      Cell_tab<-list()
    #Cell_tab_curr<-list()
      Cell_tab[["Econc"]]<-c()
      Cell_tab[["kcat"]]<-c()
      Cell_tab[["kf"]]<-c()
      Cell_tab[["fit"]]<-c()
      Cell_tab[["w"]]<-c()
      Cell_tab[["N"]]<-c()
    #Resident phenotype
      for (i in 1:1){
        Cell_tab[["Econc"]]<-c(Cell_tab[["Econc"]],E_tot_conc_i)
        Cell_tab[["kcat"]]<-c(Cell_tab[["kcat"]],kcat_i)
        Cell_tab[["kf"]]<-c(Cell_tab[["kf"]],kf_i)
        Cell_tab[["fit"]]<-c(Cell_tab[["fit"]],1)
        Cell_tab[["w"]]<-c(Cell_tab[["w"]],1)
        Cell_tab[["N"]]<-c(Cell_tab[["N"]],Ne)
      }

      
      for (t in 1:G){
        ##Absolute fitness calculation
        for (i in 1:length(Cell_tab[["N"]])){
          #if(Cell_tab[["N"]][i]>100){
          #  Cell_tab[["fit"]][i]<-fit(Cell_tab[["kf"]][i],kr_i,Cell_tab[["kcat"]][i],Cell_tab[["Econc"]][i],cost,100)
          #}
          #else{
          #  Cell_tab[["fit"]][i]<-fit(Cell_tab[["kf"]][i],kr_i,Cell_tab[["kcat"]][i],Cell_tab[["Econc"]][i],cost,Cell_tab[["N"]][i])
          #}
          Cell_tab[["fit"]][i]<-fit(Cell_tab[["kf"]][i],kr_i,Cell_tab[["kcat"]][i],Cell_tab[["Econc"]][i],cost,1)
        }
        if(Cell_tab[["fit"]][i]>10^-3){
          print(paste(Cell_tab[["fit"]][i],Cell_tab[["kf"]][i],Cell_tab[["kcat"]][i],Cell_tab[["Econc"]][i]))
        }
        #Mean absolute fitness
        w_m<-weighted.mean(Cell_tab[["fit"]],Cell_tab[["N"]])
      ##Relative fitness calculation
        for (i in 1:length(Cell_tab[["N"]])){
          Cell_tab[["w"]][i]<-Cell_tab[["fit"]][i]/w_m
        }
        #print(length(Cell_tab[["w"]]))
        ##Wright-Fisher sampling
        Cell_tab[["N"]]<-as.vector(rmultinom(1, Ne, prob=Cell_tab[["w"]]*Cell_tab[["N"]]))
        i=1
        while (i < length(Cell_tab[["N"]])+1){
          if(Cell_tab[["N"]][i]==0){
            Cell_tab[["Econc"]]<-Cell_tab[["Econc"]][-i]
            Cell_tab[["kcat"]]<-Cell_tab[["kcat"]][-i]
            Cell_tab[["kf"]]<-Cell_tab[["kf"]][-i]
            Cell_tab[["fit"]]<-Cell_tab[["fit"]][-i]
            Cell_tab[["w"]]<-Cell_tab[["w"]][-i]
            Cell_tab[["N"]]<-Cell_tab[["N"]][-i]
          }
          else{
            i=i+1
          }
        }
        mut_table<-c()
        for (i in 1:length(Cell_tab[["N"]])){
          mut_table<-c(mut_table,rbinom(1,Cell_tab[["N"]][i],mu))
        }
        pop_Size_prev<-length(Cell_tab[["N"]])
        for (i in 1:pop_Size_prev){
          if(mut_table[i]>0){
            Cell_tab[["N"]][i]=Cell_tab[["N"]][i]-mut_table[i]
            for (m in 1:mut_table[i]){
              Par_mut<-rbinom(1,1,1/2)##Parameters affected by the mutation; 0: concentration, 1: efficiency
              if(mut_bias=="increasing"){
                if(Par_mut==1){
                  mut_delta_kcat<-rnorm(1,mean=-0.2+mut_bias_coeff*log10(Cell_tab[["kcat"]][i]/10^-4),sd=0.30)
                  mut_delta_kf<-rnorm(1,mean=-0.2+mut_bias_coeff*log10(Cell_tab[["kf"]][i]),sd=0.30)
                  Cell_tab[["kcat"]][length(Cell_tab[["N"]])+1]<-10^(log10(Cell_tab[["kcat"]][i])+mut_delta_kcat)
                  Cell_tab[["kf"]][length(Cell_tab[["N"]])+1]<-10^(log10(Cell_tab[["kf"]][i])+mut_delta_kf)
                  Cell_tab[["Econc"]][length(Cell_tab[["N"]])+1]<-Cell_tab[["Econc"]][i]
                }
                else{
                  mut_delta_enz<-rnorm(1,mean=0,sd=0.5)
                  Cell_tab[["Econc"]][length(Cell_tab[["N"]])+1]<-10^(log10(Cell_tab[["Econc"]][i])+mut_delta_enz)
                  Cell_tab[["kcat"]][length(Cell_tab[["N"]])+1]<-Cell_tab[["kcat"]][i]
                  Cell_tab[["kf"]][length(Cell_tab[["N"]])+1]<-Cell_tab[["kf"]][i]
                }
              }
              else{
                if(Par_mut==1){
                mut_delta<-rnorm(2,mean=mut_bias,sd=0.30)
                Cell_tab[["kcat"]][length(Cell_tab[["N"]])+1]<-10^(log10(Cell_tab[["kcat"]][i])+mut_delta[1])
                Cell_tab[["kf"]][length(Cell_tab[["N"]])+1]<-10^(log10(Cell_tab[["kf"]][i])+mut_delta[2])
                Cell_tab[["Econc"]][length(Cell_tab[["N"]])+1]<-Cell_tab[["Econc"]][i]
                }
                else{
                  mut_delta_enz<-rnorm(1,mean=0,sd=0.5)
                  Cell_tab[["Econc"]][length(Cell_tab[["N"]])+1]<-10^(log10(Cell_tab[["Econc"]][i])+mut_delta_enz)
                  Cell_tab[["kcat"]][length(Cell_tab[["N"]])+1]<-Cell_tab[["kcat"]][i]
                  Cell_tab[["kf"]][length(Cell_tab[["N"]])+1]<-Cell_tab[["kf"]][i]
                }
              }
              
              #if(Cell_tab[["kcat"]][length(Cell_tab[["N"]])]>10^6){
               # Cell_tab[["kcat"]][length(Cell_tab[["N"]])]<-10^6
              #}
              if(Cell_tab[["kf"]][length(Cell_tab[["N"]])]>10^10){
                Cell_tab[["kf"]][length(Cell_tab[["N"]])]<-10^10
              }
              if(Cell_tab[["kcat"]][length(Cell_tab[["N"]])]>10^6){
                Cell_tab[["kcat"]][length(Cell_tab[["N"]])]<-10^6
              }
              Cell_tab[["fit"]][length(Cell_tab[["N"]])+1]<-1
              Cell_tab[["w"]][length(Cell_tab[["N"]])+1]<-1
              Cell_tab[["N"]][length(Cell_tab[["N"]])+1]<-1
            }
          }
        }
        #if(t%%10^5==0){
         # print(paste(Cell_tab[["kcat"]],Cell_tab[["kf"]]))
          #print(t)
        #}
        if(t==G/10*g){
          print(paste(Cell_tab[["Econc"]][which(Cell_tab[["N"]]==max(Cell_tab[["N"]]))],Cell_tab[["kcat"]][which(Cell_tab[["N"]]==max(Cell_tab[["N"]]))],Cell_tab[["kf"]][which(Cell_tab[["N"]]==max(Cell_tab[["N"]]))],Cell_tab[["fit"]][which(Cell_tab[["N"]]==max(Cell_tab[["N"]]))]))
          evo_results[[g]][[r]]=c(Cell_tab[["Econc"]][which(Cell_tab[["N"]]==max(Cell_tab[["N"]]))],
                                  Cell_tab[["kcat"]][which(Cell_tab[["N"]]==max(Cell_tab[["N"]]))],
                                  Cell_tab[["kf"]][which(Cell_tab[["N"]]==max(Cell_tab[["N"]]))],
                                  Cell_tab[["fit"]][which(Cell_tab[["N"]]==max(Cell_tab[["N"]]))])
          print(g)
          g=g+1
        }
      }
      #evo_results[[l]][[r]]
    #print(m)
    }
#}
  #}
#}
  if(mut_bias=="increasing"){
    save(simul,evo_results,file=paste("eq_",flux,"_Ne",Ne,"_mu",-log10(mu),"_bi0",3,"EnzVarCostCrow.rda",sep=""))
  }
  if(mut_bias!="increasing"){
    save(simul,evo_results,file=paste("eq_",flux,"_Ne",Ne,"_mu",-log10(mu),"_co1_",cost^-1,"EnzVarCostCrow.rda",sep=""))
  }

  
  
  
  
#max_fit<-0
#for(i in 1:1000){
 # cur<-fit(10^8,10^3,10^4,10^-6,10^-3,1)
 # if(cur>max_fit){
 #   max_fit<-cur
 # }
#}
 #  max_fit

#fitnz<-c()

#fitnz<-fit(10^8,10^3,10^4,10^-6,10^-3,100)

#mean(fitnz/max_fit)
