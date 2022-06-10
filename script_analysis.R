library(here)
library(ggplot2)
library(RColorBrewer)
library(plsgenomics)
library(animation)
library(tidyverse)
library(stringr)
library(plyr)
library(dplyr)

#####Scenario 1#####

#Creation of the mean dispersal data table over time


files_names <- list.files()
nb_files <- length(files_names)
data_names <- vector("list",length=nb_files)
for (i in 1 : nb_files){
  data_names[i] <- strsplit(files_names[i], split=".txt")
}

dim=50
ntime=(length(read.table(files_names[1],h=F)[,1]))/dim
nb_replicat=20
nstock=10


calcul_disp_moy = function(data,ntime,dim){
  dispmoy=c()
  for (l in 1:ntime){
    vecmat=as.vector(as.matrix(data[((l-1)*dim+1):(((l-1)*dim+1)+(dim-1)),]))
    vecmoy=c()
    for (m in 1 : length(vecmat)){
      if (vecmat[m]!=0){
        if(vecmat[m]!=-1){
          vecmoy=c(vecmoy,vecmat[m])
        }}
    }
    dispmoy=c(dispmoy,mean(vecmoy))
  }
  dispmoy
}


dispmoytot=c()
for (i in 1:length(data_names)){
  dispmoytot=c(dispmoytot,calcul_disp_moy(read.table(files_names[i]),ntime,dim))
  print((i*100)/length(data_names))
}

time=as.numeric(c(rep(c(1,(1:(ntime-1))*nstock),length(files_names))))


datadispmoy=data.frame(dispmoy=as.numeric(dispmoytot),fragmentation=rep(0,length(files_names)*ntime),aggregation=rep(0,length(files_names)*ntime),replicat=rep(0,length(files_names)*ntime),time=as.numeric(time))

for (q in 1:length(files_names)){
  
  dfrag=str_locate(files_names[q],"frag")
  ffrag=str_locate(files_names[q],"_paggr")
  daggr=str_locate(files_names[q],"_paggr")
  faggr=str_locate(files_names[q],"_rep")
  drep=str_locate(files_names[q],"_rep")
  frep=str_locate(files_names[q],".txt")
  
  datadispmoy[(q*ntime-ntime+1):(q*ntime),2]=rep(as.numeric(str_sub(files_names[q],(dfrag[1,2]+1),(ffrag[1,1]-1))),ntime)
  datadispmoy[(q*ntime-ntime+1):(q*ntime),3]=rep(as.numeric(str_sub(files_names[q],(daggr[1,2]+1),(faggr[1,1]-1))),ntime)
  datadispmoy[(q*ntime-ntime+1):(q*ntime),4]=rep(str_sub(files_names[q],(drep[1,2]+1),(frep[1,1]-1)),ntime)
  
  
}

#Summary of the mean dispersal for each combination of parameters

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


dispmoysumm=data_summary(datadispmoy, varname="dispmoy", groupnames=c("fragmentation","aggregation","time"))

write.table(x = dispmoysumm, "dispmoysumm.txt")

#Occurrence of each dispersal strategy over time

calcul_occ_disp = function(data, ntime, dim, pasdetemps){
  maxdisp=max(data)
  time=c(1:ntime)[c(1:ntime)%%pasdetemps==0]
  occdisp=c(rep(0,length(time)*maxdisp))
  temps=c()
  dispersion=c()
  for (l in time){
    for (i in 1:maxdisp){
      occdisp[maxdisp*((l/pasdetemps)-1)+i]=sum(data[((l-1)*dim+1):(((l-1)*dim+1)+(dim-1)),]==i)
      temps=c(temps,l)
      dispersion=c(dispersion,i)
    }}
  table_occ_disp=cbind(temps,dispersion,occdisp)
  table_occ_disp
}


pmut = c(0,1,10,100)
pfrag = c(1,10,100)
paggr = c(0,10,20,40,60,80)
repli = c(1:40)

dim = 50
pasdetemps = 10

data_occ=c()

for (i in 1:length(pmut)){
  for (j in 1:length(pfrag)){
    for (k in 1:length(paggr)){
      
      conditions = paste0("pmut",pmut[i],"pfrag",pfrag[j],"_paggr",paggr[k],"_rep")
      
      for (g in 1:length(repli)){
        
        ntime=length(read.table(paste0(conditions,repli[g],".txt"))[,1])/50
        data_occ_trans=data.frame(calcul_occ_disp(read.table(paste0(conditions,repli[g],".txt")),ntime,dim,pasdetemps))
        data_occ_trans2=cbind(data_occ_trans,replicats=rep(g,length(data_occ_trans[,1])), pmut=rep(pmut[i],length(data_occ_trans[,1]))
                              ,pfragr=rep(pfrag[j],length(data_occ_trans[,1])),paggr=rep(paggr[k],length(data_occ_trans[,1])))
        data_occ=rbind(data_occ,data_occ_trans2)
      }}}}

write.table(data_occ, "data_occ.txt")


#####Scenario 2#####

dim=50
nb_replicat=40
nstock=20

files_names <- list.files()
nb_files <- length(files_names)
data_names <- vector("list",length=nb_files)
for (i in 1 : nb_files){
  data_names[i] <- strsplit(files_names[i], split=".txt")
}


#Remove of STOP

for (j in 1:length(data_names)){
  if (str_detect(files_names[j],"STOP")==TRUE){
    nom=files_names[j]
    dSTOP=str_locate(nom,"STOP")
    fSTOP=str_locate(nom,".txt")
    newname=paste0(str_sub(nom,start = 1, end = dSTOP[1]-2),".txt")
    file.rename(nom,newname)
  }
}

files_names <- list.files()
nb_files <- length(files_names)
data_names <- vector("list",length=nb_files)
for (i in 1 : nb_files){
  data_names[i] <- strsplit(files_names[i], split=".txt")
}

#Creation of the mean dispersal data table over time

calcul_disp_moy = function(data,ntime,dim){
  dispmoy=c()
  for (l in 1:ntime){
    vecmat=as.vector(as.matrix(data[((l-1)*dim+1):(((l-1)*dim+1)+(dim-1)),]))
    vecmoy=c()
    for (m in 1 : length(vecmat)){
      if (vecmat[m]!=0){
        if(vecmat[m]!=-1){
          vecmoy=c(vecmoy,vecmat[m])
        }}
    }
    dispmoy=c(dispmoy,mean(vecmoy))
  }
  dispmoy
}


dispmoytot=c()
time=c()
ntimetot=c()
for (i in 1:length(data_names)){
  ntime=(length(read.table(files_names[i],h=F)[,1]))/dim
  ntimetot=c(ntimetot,ntime)
  dispmoytot=c(dispmoytot,calcul_disp_moy(read.table(files_names[i]),ntime,dim))
  time=c(time,c((1:(ntime))*nstock))
  print((i*100)/length(data_names))
}

datadispmoy=data.frame(dispmoy=as.numeric(dispmoytot),pmut=rep(0,length(dispmoytot)),pfrag=rep(0,length(dispmoytot)),aggregation=rep(0,length(dispmoytot)),replicat=rep(0,length(dispmoytot)),time=as.numeric(time))

for (q in 1:length(files_names)){
  
  dpmut=str_locate(files_names[q],"pmut")
  fpmut=str_locate(files_names[q],"pfrag")
  dpfrag=str_locate(files_names[q],"pfrag")
  fpfrag=str_locate(files_names[q],"_paggr")
  daggr=str_locate(files_names[q],"_paggr")
  faggr=str_locate(files_names[q],"_rep")
  drep=str_locate(files_names[q],"_rep")
  frep=str_locate(files_names[q],".txt")
  
  datadispmoy[(sum(ntimetot[1:q])-ntimetot[q]+1):(sum(ntimetot[1:q])),2]=rep(as.numeric(str_sub(files_names[q],(dpmut[1,2]+1),(fpmut[1,1]-1))),ntimetot[q])
  datadispmoy[(sum(ntimetot[1:q])-ntimetot[q]+1):(sum(ntimetot[1:q])),3]=rep(as.numeric(str_sub(files_names[q],(dpfrag[1,2]+1),(fpfrag[1,1]-1))),ntimetot[q])
  datadispmoy[(sum(ntimetot[1:q])-ntimetot[q]+1):(sum(ntimetot[1:q])),4]=rep(as.numeric(str_sub(files_names[q],(daggr[1,2]+1),(faggr[1,1]-1))),ntimetot[q])
  datadispmoy[(sum(ntimetot[1:q])-ntimetot[q]+1):(sum(ntimetot[1:q])),5]=rep(str_sub(files_names[q],(drep[1,2]+1),(frep[1,1]-1)),ntimetot[q])
  
  
}

#Number of individuals over time

calcul_nb_ind = function(data,ntime,dim){
  nb_ind=c()
  for (l in 1:ntime){
    vecmat=as.vector(as.matrix(data[((l-1)*dim+1):(((l-1)*dim+1)+(dim-1)),]))
    vecmatind=vecmat[vecmat != 0 & vecmat != -1]
    nb_ind=c(nb_ind,length(vecmatind))
  }
  nb_ind
}


nb_ind_tot=c()
time=c()
ntimetot=c()
for (i in 1:length(data_names)){
  ntime=(length(read.table(files_names[i],h=F)[,1]))/dim
  ntimetot=c(ntimetot,ntime)
  nb_ind_tot=c(nb_ind_tot,calcul_nb_ind(read.table(files_names[i]),ntime,dim))
  time=c(time,c((1:(ntime))*nstock))
  print((i*100)/length(data_names))
}

datatot=cbind(datadispmoy,nb_ind=nb_ind_tot)

#Number of empty patch over time

calcul_nb_vide = function(data,ntime,dim){
  nb_vide=c()
  for (l in 1:ntime){
    vecmat=as.vector(as.matrix(data[((l-1)*dim+1):(((l-1)*dim+1)+(dim-1)),]))
    vecmatvide=vecmat[vecmat == 0]
    nb_vide=c(nb_vide,length(vecmatvide))
  }
  nb_vide
}

nb_vide_tot=c()
time=c()
ntimetot=c()
for (i in 1:length(data_names)){
  ntime=(length(read.table(files_names[i],h=F)[,1]))/dim
  ntimetot=c(ntimetot,ntime)
  nb_vide_tot=c(nb_vide_tot,calcul_nb_vide(read.table(files_names[i]),ntime,dim))
  time=c(time,c((1:(ntime))*nstock))
  print((i*100)/length(data_names))
}

datatot=cbind(datatot,nb_vide=nb_vide_tot)


poccupation=datatot$nb_ind/(datatot$nb_ind+datatot$nb_vide)

datatot=cbind(datatot,poccupation=poccupation)

#Time and percentage of fragmentation at extinction

levfrag=levels(as.factor(datatot$pfrag))
levmut=levels(as.factor(datatot$pmut))
levaggr=levels(as.factor(datatot$aggregation))
levrep=levels(as.factor(datatot$replicat))

tempsext=c()
pfrag=c()
pmut=c()
aggregation=c()
replicat=c()
fragext=c()

for ( i in 1:length(levfrag)){
  for ( j in 1:length(levmut)){
    for (k in 1:length(levaggr)){
      for (l in 1:length(levrep)){
        tempsext=c(tempsext,max(datatot$time[datatot$pmut==levmut[j] & datatot$pfrag==levfrag[i] & datatot$aggregation ==levaggr[k] & datatot$replicat==levrep[l]]))
        tempsexti=max(datatot$time[datatot$pmut==levmut[j] & datatot$pfrag==levfrag[i] & datatot$aggregation ==levaggr[k] & datatot$replicat==levrep[l]])
        nb_frag=2500-(datatot$nb_vide[datatot$pmut==levmut[j] & datatot$pfrag==levfrag[i] & datatot$aggregation ==levaggr[k] & datatot$replicat==levrep[l] & datatot$time==tempsexti] + datatot$nb_ind[datatot$pmut==levmut[j] & datatot$pfrag==levfrag[i] & datatot$aggregation ==levaggr[k] & datatot$replicat==levrep[l] & datatot$time==tempsexti] )
        fragext = c(fragext,((nb_frag*100)/2500))
        pfrag=c(pfrag,levfrag[i])
        pmut=c(pmut,levmut[j])
        aggregation=c(aggregation,levaggr[k])
        replicat=c(replicat,levrep[l])
      }
    }
  }
}

datatempsext=data.frame(pfrag,pmut,aggregation,replicat,tempsext,fragext)

#Summary each combination of parameters

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


dispmoysumm=data_summary(datatot, varname="dispmoy", groupnames=c("pmut","pfrag","aggregation","time"))
nbindsumm=data_summary(datatot, varname="nb_ind", groupnames=c("pmut","pfrag","aggregation","time"))
poccupasumm=data_summary(datatot, varname="poccupation", groupnames=c("pmut","pfrag","aggregation","time"))
tempsextsumm=data_summary(datatempsext, varname="tempsext", groupnames=c("pmut","pfrag","aggregation"))
fragextsumm=data_summary(datatempsext, varname="fragext", groupnames=c("pmut","pfrag","aggregation"))

write.table(dispmoysumm, "dispmoysumm.txt")
write.table(nbindsumm, "nbindsumm.txt")
write.table(poccupasumm, "poccupasumm.txt")
write.table(tempsextsumm, "tempsextsumm.txt")
write.table(fragextsumm, "fragextsumm.txt")