 ###### Script of the model of dispersal evolution in fragmented landscape ######
 ###### Basile Finand
 
 #Packages

library(doParallel)
library(NLMR) 
library(raster)
library(lattice)
library(landscapetools) 
library(plsgenomics)


detectCores()             
cl<-makeCluster(20)       
registerDoParallel(cl)

#Parameters of the model

pext=0.05         #Extinction probability
dim=50            #Dimension of the landscape
nstartcol=10      #Number of populations at the start of the simulation
dispstartcol=12   #colonisation capacity of the first populations
ntime=50000       #Number of time steps of the simulation
pmut=0.10         #Mutation rate
nb_replicas=20    #Number of replicates
nstock=10         #Number of time steps for data storage

frag_change=0     #0 for the Scenario 1: Evolution of dispersal in fixed fragmented landscapes, 1 for the Scenario 2: Evolutionary rescue under progressively increasing fragmentation

pfrag=0.01        #If scenario 2, fragmentation rate
frag_stop=1       #If scenario 2, maximal fragmentation to stop the simulation
temps_stab=200    #f scenario 2, time to stabilize before increase of fragmentation

mod_frag=c(0,0.20,0.40,0.60,0.80,0.90,0.95,0.99)      #Percentage of fragmentation to be tested
mod_aggr=c(0,0.20,0.40,0.60,0.80)                     #Percentage of aggregation to be tested  



#Function to extend the landscape to do the Torus

expansion=function(x){
  maxcol=max(x);
  extension=matrix(0,nrow=dim+maxcol*2,ncol=dim+maxcol*2);
  extension[maxcol+1:dim,maxcol+1:dim]=x;
  extension[1:maxcol,maxcol+1:dim]=x[(dim-maxcol+1):dim,];
  extension[(dim+maxcol+1):(dim+maxcol*2),maxcol+1:dim]=x[1:maxcol,];
  extension[maxcol+1:dim,1:maxcol]=x[,(dim-maxcol+1):dim];
  extension[maxcol+1:dim,(dim+maxcol+1):(dim+maxcol*2)]=x[,1:maxcol];
  extension[1:maxcol,1:maxcol]=x[(dim-maxcol+1):dim,(dim-maxcol+1):dim];
  extension[(dim+maxcol+1):(dim+maxcol*2),(dim+maxcol+1):(dim+maxcol*2)]=x[1:maxcol,1:maxcol];
  extension[1:maxcol,(dim+maxcol+1):(dim+maxcol*2)]=x[(dim-maxcol+1):dim,1:maxcol];
  extension[(dim+maxcol+1):(dim+maxcol*2),1:maxcol]=x[1:maxcol,(dim-maxcol+1):dim];
  extension
}

#Mutation function

mutation<-function(x){
  if (runif(1,min=0,max=1)<pmut) {
    if(x==1){
      x=x+1
    }else{
      x=x+sample(c(-1,1),1)
    }}                                 
  x
}

#Function to find all patches defined at a x distance of another one

patch=function(x){
  coorpatch=(-x+1):(x-1)
  nb=length((-x+1):(x-1))
  patchdisp1=data.frame(x=c(rep(0,nb)),y=c(rep(0,nb)))
  for (i in 1:nb){
    patchdisp1[i,1]=x
    patchdisp1[i,2]=coorpatch[i]
  }
  patchdisp2=data.frame(x=c(rep(0,nb)),y=c(rep(0,nb)))
  for (j in 1:nb){
    patchdisp2[j,1]=-x
    patchdisp2[j,2]=coorpatch[j]
  }
  patchdisp3=data.frame(x=c(rep(0,nb)),y=c(rep(0,nb)))
  for (k in 1:nb){
    patchdisp3[k,1]=coorpatch[k]
    patchdisp3[k,2]=x
  }
  patchdisp4=data.frame(x=c(rep(0,nb)),y=c(rep(0,nb)))
  for (l in 1:nb){
    patchdisp4[l,1]=coorpatch[l]
    patchdisp4[l,2]=-x
  }
  patchdisp5=data.frame(x=c(rep(0,4)),y=c(rep(0,4)))
  patchdisp5[1,1]=x
  patchdisp5[1,2]=x
  patchdisp5[2,1]=-x
  patchdisp5[2,2]=-x
  patchdisp5[3,1]=-x
  patchdisp5[3,2]=x
  patchdisp5[4,1]=x
  patchdisp5[4,2]=-x
  patchdisp=rbind(patchdisp1,patchdisp2,patchdisp3,patchdisp4,patchdisp5)
  patchdisp
}

#Function to calculate the mean dispersal on a grid

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


#########Start of the simulation #######

nouvfragvec=c()
fragreelvec=c()

for (b in 1:length(mod_frag)){      
  for (s in 1:length(mod_aggr)){    
    
    
        foreach (c = 1:nb_replicas)  %dopar% {  

      library(doParallel)
      library(NLMR)
      library(raster)
      library(lattice)
      library(landscapetools)
      library(plsgenomics)
      
      frag=mod_frag[b]          #Definition of the percentage of fragmentation
      paggr=mod_aggr[s]         #Definition of the percentage of the aggregaion    
      if (paggr == 0 ) {presence_aggr=0} else {presence_aggr=1} 
      
      #Landscape creation with aggregation
      
      if(presence_aggr==1){
        if((paggr*2)>=1.6){
          aggr=nlm_fbm(ncol = dim, nrow=dim,fract_dim = paggr*2, modus_operandi="sloppy")}
        else{
          aggr=nlm_fbm(ncol = dim, nrow=dim,fract_dim = paggr*2)}  
        
        landbyn=util_classify(x=aggr,weighting=c(1-frag,frag))  
        landaggr=matrix(landbyn,dim,dim)                        
        landaggr[landaggr==1]=0                                 
        landaggr[landaggr==2]=-1                                
        
        for (v in 1:nstartcol){               
          nrow=sample(1:dim,1)                
          ncol=sample(1:dim,1)                
          while (landaggr[nrow,ncol]!=0){     
            nrow=sample(1:dim,1)
            ncol=sample(1:dim,1)
          }
          landaggr[nrow,ncol]=dispstartcol}   
      }
      
      #Landscape creation without aggregation
      
      if (presence_aggr==0){
        distr=c(rep(-1,frag*(dim*dim)), rep(dispstartcol,nstartcol), rep(0, dim*dim-frag*(dim*dim)-nstartcol))  
        distrrand=sample(distr,length(distr), replace=FALSE)                                                   
        mdistr=matrix(data=distrrand,nrow=dim, ncol=dim)                                                       
      }
      
      
      if(presence_aggr==0){mdistr=mdistr}
      if(presence_aggr==1){mdistr=landaggr}  
      
      
      
      #Start of the time if scenario 1
      
      if (frag_change==0){
        
        file_name=paste0("frag", frag*100,"_paggr",paggr*100,"_rep",c,".txt")
        write.table(x = mdistr,file = file_name, row.names = FALSE, col.names=FALSE)
        print(file_name)
        
        
        for (time in 1:ntime){
          
          mdistrexp=expansion(mdistr)   
          mdistrexpbis= mdistrexp       
          maxcol=max(mdistr)            
          
          #Competition/colonization process
          
          for (i in (maxcol+1):(dim+maxcol)){                                                 
            for (j in (maxcol+1):(dim+maxcol)){                                               
              if (mdistrexp[i,j]==0){                                                         
                colonizer=c()
                for (k in 1:maxcol){                                                          
                  potpatch=patch(k)                                                           
                  potcol=c()
                  for (l in 1:length(potpatch[,1])){                                          
                    if (mdistrexp[i+potpatch[l,1],j+potpatch[l,2]]>=k){                       
                      potcol=c(potcol,mdistrexp[i+potpatch[l,1],j+potpatch[l,2]])}            
                  }
                  if (is.element(k,potcol)==TRUE){                                            
                    colonizer=k                                                               
                    break}                                                                    
                  colonizer=c(colonizer,potcol)                                               
                }
                if (is.null(colonizer)== FALSE){mdistrexpbis[i,j]=mutation(min(colonizer))}   #Mutation process
              }}}
          mdistr=mdistrexpbis[(maxcol+1):(dim+maxcol),(maxcol+1):(dim+maxcol)]                
          
          #Extinction process
          
          for (i in 1:dim){
            for (j in 1:dim){                                     
              if (mdistr[i,j]!= 0) {                              
                if (mdistr[i,j]!= -1){                            
                  if (runif(1,min=0,max=1)<pext) mdistr[i,j]=0}   
              }}}
          
          #For each storage time, recording the landscape
          
          if (time%%nstock==0){
            write.table(x = mdistr,file = file_name,append = TRUE, row.names = FALSE, col.names=FALSE)  
          }
          
          if(sum(sapply(c(0,-1),function(x) sum(mdistr==x)))==dim*dim){                               
            if (time%%nstock!=0){
              write.table(x = mdistr,file = file_name,append = TRUE, row.names = FALSE, col.names=FALSE)  
            }
            file_nameSTOP = paste0("pmut", pmut*1000,"pfrag",pfrag*10000,"_paggr",paggr*100,"_rep",c,"_STOP",time,".txt")  
            file.rename(file_name,file_nameSTOP)
          } 
          
          print((time/ntime)*100)
          
          
          #If all patch are empty due to extinction, end of the simulation
          
          if (sum(sapply(c(0,-1),function(x) sum(mdistr==x)))==dim*dim){break} 
        }}
      
      
      
      #Start of the time if scenario 2
      
      
      if (frag_change==1){
        
        file_name=paste0("pmut", pmut*1000,"pfrag",pfrag*10000,"_paggr",paggr*100,"_rep",c,".txt")
        write.table(x = mdistr,file = file_name, row.names = FALSE, col.names=FALSE)
        print(file_name)
        time=0
        nouvfrag=frag
        dispmoyvec=c()
        
    
        while (nouvfrag<frag_stop){
          
          time=time+1
          
          #stabilization of the landscape before the increase of fragmentation
          
          if(time <= temps_stab){
            
            mdistrexp=expansion(mdistr)   
            mdistrexpbis= mdistrexp       
            maxcol=max(mdistr)      
            
            #Competition/colonization process
            
            for (i in (maxcol+1):(dim+maxcol)){                                                 
              for (j in (maxcol+1):(dim+maxcol)){                                              
                if (mdistrexp[i,j]==0){                                                         
                  colonizer=c()
                  for (k in 1:maxcol){                                                          
                    potpatch=patch(k)                                                           
                    potcol=c()
                    for (l in 1:length(potpatch[,1])){                                         
                      if (mdistrexp[i+potpatch[l,1],j+potpatch[l,2]]>=k){                       
                        potcol=c(potcol,mdistrexp[i+potpatch[l,1],j+potpatch[l,2]])}            
                    }
                    if (is.element(k,potcol)==TRUE){                                            
                      colonizer=k                                                               
                      break}                                                                   
                    colonizer=c(colonizer,potcol)                                              
                  }
                  if (is.null(colonizer)== FALSE){mdistrexpbis[i,j]=min(colonizer)}    
                }}}
            mdistr=mdistrexpbis[(maxcol+1):(dim+maxcol),(maxcol+1):(dim+maxcol)]                
            
            #Extinction process
            
            for (i in 1:dim){
              for (j in 1:dim){                                    
                if (mdistr[i,j]!= 0) {                              
                  if (mdistr[i,j]!= -1){                            
                    if (runif(1,min=0,max=1)<pext) mdistr[i,j]=0}   
                }}}
            
            #For each storage time, recording the landscape
            
            if (time%%nstock==0){
              write.table(x = mdistr,file = file_name,append = TRUE, row.names = FALSE, col.names=FALSE)  
            }
            
            if(sum(sapply(c(0,-1),function(x) sum(mdistr==x)))==dim*dim){                           
              if (time%%nstock!=0){
                write.table(x = mdistr,file = file_name,append = TRUE, row.names = FALSE, col.names=FALSE)  
              }
              file_nameSTOP = paste0("pmut", pmut*1000,"pfrag",pfrag*10000,"_paggr",paggr*100,"_rep",c,"_STOP",time,".txt")  
              file.rename(file_name,file_nameSTOP)
            } 
            
            
            if (sum(sapply(c(0,-1),function(x) sum(mdistr==x)))==dim*dim){break}
            
            
          #start of the increase of fragmentation
              
          }else{
          
          mdistrexp=expansion(mdistr)   
          mdistrexpbis= mdistrexp       
          maxcol=max(mdistr)         
          
          #Competition/colonization process
          
          for (i in (maxcol+1):(dim+maxcol)){                                                 
            for (j in (maxcol+1):(dim+maxcol)){                                               
              if (mdistrexp[i,j]==0){                                                         
                colonizer=c()
                for (k in 1:maxcol){                                                         
                  potpatch=patch(k)                                                          
                  potcol=c()
                  for (l in 1:length(potpatch[,1])){                                         
                    if (mdistrexp[i+potpatch[l,1],j+potpatch[l,2]]>=k){                       
                      potcol=c(potcol,mdistrexp[i+potpatch[l,1],j+potpatch[l,2]])}            
                  }
                  if (is.element(k,potcol)==TRUE){                                            
                    colonizer=k                                                               
                    break}                                                                    
                  colonizer=c(colonizer,potcol)                                               
                }
                if (is.null(colonizer)== FALSE){mdistrexpbis[i,j]=mutation(min(colonizer))}   #Mutation process
                }}}
          mdistr=mdistrexpbis[(maxcol+1):(dim+maxcol),(maxcol+1):(dim+maxcol)]                
          
          #Patchs extinction
          
          for (i in 1:dim){
            for (j in 1:dim){                                     
              if (mdistr[i,j]!= 0) {                              
                if (mdistr[i,j]!= -1){                           
                  if (runif(1,min=0,max=1)<pext) mdistr[i,j]=0}   
              }}}
          
          #For each storage time, recording the landscape
          
          if (time%%nstock==0){
            write.table(x = mdistr,file = file_name,append = TRUE, row.names = FALSE, col.names=FALSE)  
          }
          
          if(sum(sapply(c(0,-1),function(x) sum(mdistr==x)))==dim*dim){                               
            if (time%%nstock!=0){
              write.table(x = mdistr,file = file_name,append = TRUE, row.names = FALSE, col.names=FALSE)  
            }
            file_nameSTOP = paste0("pmut", pmut*1000,"pfrag",pfrag*10000,"_paggr",paggr*100,"_rep",c,"_STOP",time,".txt")  
            file.rename(file_name,file_nameSTOP)
          } 
          
          
          print((time/ntime)*100)
          
          if (sum(sapply(c(0,-1),function(x) sum(mdistr==x)))==dim*dim){break}  
          
          
          #If scenario 2, increase of fragmentation process
          #When there is no aggregation
          
          if(presence_aggr==0){
            nouvfrag=(1-nouvfrag)*pfrag+nouvfrag
            for (k in 1:dim){
              for (l in 1:dim){
                if (mdistr[k,l]!=-1){
                  if (runif(1,min=0,max=1)<pfrag){ mdistr[k,l]=-1}}}}
          }
          
          
          #when there is aggregation
          
          if(presence_aggr==1){
            nouvfrag=(1-nouvfrag)*pfrag+nouvfrag
            landbyn=util_classify(x=aggr,weighting=c(1-(nouvfrag),(nouvfrag)))  
            landaggr=matrix(landbyn,dim,dim)                        
            landaggr[landaggr==1]=0                                 
            landaggr[landaggr==2]=-1 
              for (o in 1:dim){
              for (p in 1:dim){
                if (landaggr[o,p]==-1){mdistr[o,p]=landaggr[o,p]}}}
          
          nouvfragvec=c(nouvfragvec,nouvfrag)}
          fragreel=(sum(mdistr==-1)/(dim*dim))
          fragreelvec=c(fragreelvec,fragreel)
          dispmoyvec=c(dispmoyvec,calcul_disp_moy(mdistr,1,dim))
        }}}}}}

stopCluster(cl)

