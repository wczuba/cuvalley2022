# DAta LOAD , EDA, SUBSET & transform for training dataset


# Library ----

library("lubridate")
library("hexbin")
library("latticeExtra")
library("dplyr")
library("pheatmap")
library("caTools")
library("gridExtra")

# Data load ----

data00<-read.csv("D:/CUVALLEY2/CASE_3/_data/_proc_data/input_data.csv",stringsAsFactors=FALSE)
temp_mes<-read.csv2("D:/CUVALLEY2/CASE_3/_data/temp_zuz.csv",stringsAsFactors=FALSE)
ident<-read.csv2("D:/CUVALLEY2/CASE_3/_data/opis_zmiennych.csv",stringsAsFactors=FALSE)


# Column names & type----

u=0
for (p in c(1:nrow(ident)))
{
  kolumny<-colnames(data00)
  z<-grep(toupper(ident[p,2]), toupper(kolumny), value = F)
  print(z)
  colnames(data00)[z]<-paste0(ident[p,3],"_",ident[p,4])
}


# Time convert & merge data ----

data00$czas<-substr(data00$czas,1,19)
data01<-merge(data00, temp_mes, by.x = "czas", by.y = "Czas",all.x=T)
data01<-data01[substr(data01$czas,1,2)=="20",]
data01$czas<-as.POSIXct(data01$czas)

# New temp. columns ----

data02<-data01
data02$dT<-NA
data02$odst_pom<-NA
data02[is.na(data02$temp_zuz)==F,"odst_pom"]<-as.numeric(c(1,diff(data02[is.na(data02$temp_zuz)==F,"czas"])))
data02[is.na(data02$temp_zuz)==F,"dT"]<-c(0,diff(data02[is.na(data02$temp_zuz)==F,"temp_zuz"]))/as.numeric(c(1,diff(data02[is.na(data02$temp_zuz)==F,"czas"])))*3600
for (i in c(2:58))
{
  data02[,i]<-as.numeric(data02[,i])
}



# New mass balance columns ----

data02$M_N_c<-(data02[,2]+data02[,3])*data02[,56]/100/60*1000
data02$M_N_Fe<-(data02[,2]+data02[,3])*data02[,57]/100/60*1000
data02$M_N_S<-(data02[,2]+data02[,3])*data02[,58]/100/60*1000
data02$M_N_R<-(data02[,2]+data02[,3])*(100-data02[,56]-data02[,57]-data02[,58])/100/60*1000

data02$M_R_Fe<-data02[,4]*data02[,54]/100/60*1000
data02$M_R_sog<-data02[,4]*data02[,55]/100/60*1000
data02$M_R_R<-(data02[,4])*(100-data02[,54]-data02[,55])/100/60*1000

data02$M_P<-data02[,5]/60*1000

data02$M_all<-data02$M_N_c+data02$M_N_Fe+data02$M_N_S+data02$M_N_R+data02$M_R_Fe+data02$M_R_sog+data02$M_R_R+data02$M_P

data02$SR_heat<-data02[,6]*(data02[,15]-data02[,52])+
  data02[,7]*(data02[,16]-data02[,52])+
  data02[,8]*(data02[,17]-data02[,52])+
  data02[,9]*(data02[,18]-data02[,52])+
  data02[,10]*(data02[,19]-data02[,52])+
  data02[,11]*(data02[,20]-data02[,52])+
  data02[,12]*(data02[,21]-data02[,52])+
  data02[,13]*(data02[,22]-data02[,52])

# New time marks  ----

data02[is.na(data02$temp_zuz)==F,"data"]<-substr(data02[is.na(data02$temp_zuz)==F,"czas"],1,10)
data02[is.na(data02$temp_zuz)==F,"godzina"]<-hour(data02[is.na(data02$temp_zuz)==F,"czas"])+minute(data02[is.na(data02$temp_zuz)==F,"czas"])/60
data02[is.na(data02$temp_zuz)==F,"dzien_tyg"]<-wday(data02[is.na(data02$temp_zuz)==F,"czas"])+data02[is.na(data02$temp_zuz)==F,"godzina"]/24
data02[is.na(data02$temp_zuz)==F,"dzien_nar"]<-as.numeric(difftime(data02[is.na(data02$temp_zuz)==F,"czas"],strptime("01.10.2020", format = "%d.%m.%Y"),units="days"))


# Temp. analysis ----

temp_mes2<-data02[is.na(data02$temp_zuz)==F,c("czas","temp_zuz","dT","odst_pom","data","godzina","dzien_tyg","dzien_nar")]



# Temp. dynamics ----

hexbinplot(jitter(temp_zuz,factor=.5)~jitter(dT,factor=0.5),aspect=1,data=temp_mes2,xlim=c(-30,30),ylim=c(1280,1330),xbins=320,trans = sqrt, inv = function(x) x^2,ylab="Temp. zuzel [st.C]",xlab="Zmiana temperatury [st.C/1h]",
           panel=function(x,y,...)
             {
             panel.hexbinplot(x,y,...)
             panel.ellipse(x,y,...,col="red",lwd=2, level = 0.5, segments = 50)
             #panel.lmline(x,y,...,col="blue",lwd=2)
             panel.rug(x,y,...,col="blue")
           })

# Ecdf for temp.  ----


w1<-ecdfplot(~odst_pom,data=temp_mes2,xlim=c(0,7200),aspect=1,lwd=2,xlab="Odstep pomiarowy [s]")
w2<-ecdfplot(~dT,data=temp_mes2,xlim=c(-100,100),aspect=1,lwd=2,xlab="Zmiana temp. dT [st C/1h]")
w3<-ecdfplot(~temp_zuz,data=temp_mes2,xlim=c(1200,1400),aspect=1,lwd=2,xlab="Temperatura żużel [st. C]")
grid.arrange(w1,w2,w3,ncol=1)

# ACF for temp.  ----

acf(temp_mes2$temp_zuz,lag.max=300)


# DAY ACF for 1h MW ----


temp_mes2<-temp_mes2[order(temp_mes2$godzina,temp_mes2$data),]
xyplot(runmean(temp_mes2$temp_zuz,11387/24)~jitter(temp_mes2$godzina,100),xlab="Czas zegarowy [h]",ylab="Temp. żużel [st.C]",
       panel=function(x,y,...)
         {
         panel.xyplot(x,y,...,type="p")
         #panel.violin(x,y,...)
         panel.xyplot(x,y,...,type="spline",lwd=3,col="red")
       })

set.seed(42)
rows <- sample(nrow(temp_mes2))
temp_mes2<-temp_mes2[rows,]
xyplot(runmean(temp_mes2$temp_zuz,11387/24)~jitter(temp_mes2$godzina,100),
       panel=function(x,y,...)
       {
         panel.xyplot(x,y,...,type="p")
         panel.xyplot(x,y,...,type="spline",lwd=2,col="red")
       })


# Weeks ACF with 1D MW----

temp_mes2<-temp_mes2[order(temp_mes2$dzien_tyg),]
xyplot(runmean(temp_mes2$temp_zuz,11387/7)~temp_mes2$dzien_tyg,ylab="Temp. żużel [st.C]",xlab="Dzień tygodnia",
       panel=function(x,y,...)
       {
         panel.xyplot(x,y,...,type="p")
         panel.xyplot(x,y,...,type="spline",lwd=2,col="red")
       })

set.seed(42)
rows <- sample(nrow(temp_mes2))
temp_mes2<-temp_mes2[rows,]
xyplot(runmean(temp_mes2$temp_zuz,11387/7)~temp_mes2$dzien_tyg,
       panel=function(x,y,...)
       {
         panel.xyplot(x,y,...,type="p")
         panel.xyplot(x,y,...,type="spline",lwd=2,col="red")
       })


# Trends for weeks MW  ----

temp_mes2<-temp_mes2[order(temp_mes2$dzien_nar),]
xyplot(runmean(temp_mes2$temp_zuz,(11387/(length(unique(temp_mes2$data)))*7))~temp_mes2$dzien_nar,xlab="Czas eksploatacji pieca",ylab="Temp. żużel [st.C]",
       panel=function(x,y,...)
       {
         panel.xyplot(x,y,...,type="p")
         panel.xyplot(x,y,...,type="spline",lwd=2,col="red")
       })

set.seed(42)
rows <- sample(nrow(temp_mes2))
temp_mes2<-temp_mes2[rows,]
xyplot(runmean(temp_mes2$temp_zuz,(11387/(length(unique(temp_mes2$data)))*7))~temp_mes2$dzien_nar,
       panel=function(x,y,...)
       {
         panel.xyplot(x,y,...,type="p")
         panel.xyplot(x,y,...,type="spline",lwd=2,col="red")
       })





# Procesing temp events (+6h;-2h TSA window)----



kol<- c(14,62:70,2:13,15:61)
temp_mes2<-temp_mes2[order(temp_mes2$dzien_nar),]



for (col_nr in kol)
{

wind1=6
wind2=2
probs=200

smr <- data.frame(matrix(NA, ncol=probs+4, nrow=0))
nazwy<-c("czas","temp_zuz","dT","odst_pom",sprintf("T%02d", seq(1,probs)))
colnames(smr)<-nazwy

smc <- data.frame(matrix(NA, ncol=6, nrow=0))
nazwy<-c("czas","temp_zuz","dT","odst_pom","time","value")
colnames(smc)<-nazwy
max_events=nrow(temp_mes2)
max_events=5000

for (n in c(10:max_events))
{

  smc[c((1+probs*(n-1)):(probs+probs*(n-1))),1]<-temp_mes2[n,"czas"]
  smc[c((1+probs*(n-1)):(probs+probs*(n-1))),2]<-temp_mes2[n,"temp_zuz"]
  smc[c((1+probs*(n-1)):(probs+probs*(n-1))),3]<-temp_mes2[n,"dT"]
  smc[c((1+probs*(n-1)):(probs+probs*(n-1))),4]<-temp_mes2[n,"odst_pom"]
  values<-data02[data02$czas>=(temp_mes2[n,"czas"]-3600*wind1)&data02$czas<=(temp_mes2[n,"czas"]+3600*wind2),c("czas",colnames(data02)[col_nr])]
  values$time<-as.numeric(values$czas-temp_mes2[n,"czas"])
  values[,2]<-as.numeric(values[,2])
  colnames(values)[2]<-"value"
  values$time<-cut(values$time,breaks=quantile(values$time,probs=seq(0,1,by=1/(probs-1))),labels=F,dig.lab = 5)
  values<-values%>% group_by(time) %>% summarise(sr=mean(value)) 
  
  smc[c((1+probs*(n-1)):(probs+probs*(n-1))),"value"]<-values$sr
  if (exists("old_values"))
  { 
    smc[c((1+probs*(n-1)):(probs+probs*(n-1))),"value_diff"]<-old_values$sr-values$sr 
    }
  old_values<-values
  smc[c((1+probs*(n-1)):(probs+probs*(n-1))),"time"]<-as.numeric(c(1:probs)/probs*(wind1+wind2)-wind1)
  #colnames(smc)<-nazwy
  smr[n,c(1:4)]<-c(temp_mes2[n,"czas"],temp_mes2[n,"temp_zuz"],temp_mes2[n,"dT"],temp_mes2[n,"odst_pom"])
  smr[n,c(5:(probs+4))]<-values$sr
  print(n)
}


smc<-na.omit(smc)
smc<-smc[smc$odst_pom<4000,]
smc<-smc[smc$dT>=-25&smc$dT<=25,]
smc<-smc[smc$temp_zuz>=1270,]





p01<-hexbinplot(value~time|cut(temp_zuz,breaks=c(quantile(smc$temp_zuz, probs = seq(0, 1, by = 0.2),na.rm=T))),ylab=colnames(data02)[col_nr],data=smc[order(smc$time),],aspect=1,layout=c(5,1),ylim=c(quantile(smc$value,probs=0.05,na.rm = T),quantile(smc$value,probs=0.95,na.rm = T)),xlab="Time [h]",main="T model",
           panel=function(x,y,...)
             {
             panel.hexbinplot(x,y,...)
             panel.abline(v=0,col="white")
             panel.xyplot(x,runmean(y,nrow(smc)/(100)),...,type="spline",col="red",lwd=3)
             panel.xyplot(x,runquantile(y,nrow(smc)/(100),probs=0.5),...,type="spline",col="blue",lwd=3)
           })


p02<-hexbinplot(value_diff~time|cut(dT,breaks=c(quantile(smc$dT, probs = seq(0, 1, by = 0.2),na.rm=T))),data=smc[order(smc$time),],aspect=1,layout=c(5,1),ylab=colnames(data02)[col_nr],ylim=c(quantile(smc$value_diff,probs=0.05,na.rm = T),quantile(smc$value_diff,probs=0.95,na.rm = T)),xlab="Time [h]",main="dT model",
           panel=function(x,y,...)
           {
             panel.hexbinplot(x,y,...)
             panel.abline(v=0,col="white")
             panel.abline(h=0,col="white")
             panel.xyplot(x,runmean(y,nrow(smc)/(100)),...,type="spline",col="red",lwd=3)
             panel.xyplot(x,runquantile(y,nrow(smc)/(100),probs=0.5),...,type="spline",col="blue",lwd=3)
           })

png(filename = paste0("hex_",col_nr,".png"),width = 1280, height = 720, units = "px")
grid.arrange(p01,p02,ncol=1)
dev.off()


saveRDS(smc, file = paste0("smc_",col_nr,".rds"))
saveRDS(smr, file = paste0("smr_",col_nr,".rds"))


smr<-na.omit(smr)

mtscaled <- as.matrix((smr[,c(5:204)]))


annotation_row = data.frame(
  Temp_zuz = cut(smr$temp_zuz,breaks=c(quantile(smr$temp_zuz, probs = seq(0, 1, by = 0.2),na.rm=T))))

rownames(annotation_row) = rownames(smr)

#png(filename = paste0("heat_temp_",col_nr,".png"),width = 1280, height = 720, units = "px")
#print(pheatmap(mtscaled, annotation_row = annotation_row,cluster_cols=F,labels_row=""))
#dev.off()

mtscaled <- as.matrix((smr[,c(5:204)]))

for (i in c(5:204))
{
smr[,i]<-c(0,diff(smr[,i]))
}

saveRDS(smr, file = paste0("smr2_",col_nr,".rds"))

mtscaled <- as.matrix((smr[,c(5:204)]))

annotation_row = data.frame(
  dT = cut(smr$dT,breaks=c(quantile(smr$dT, probs = seq(0, 1, by = 0.2),na.rm=T))))

rownames(annotation_row) = rownames(smr)

#png(filename = paste0("heat_dT_",col_nr,".png"),width = 1280, height = 720, units = "px")
#print(pheatmap(mtscaled, annotation_row = annotation_row,cluster_cols=F,labels_row=""))
#dev.off()

}


# Learning dataset creation ----

data03<-data02

data03<-data03[,c(1,14,62:71,73:74,59:61)]
colnames(data03)[2]<-"Moc_ciep"

for (i in c(1:ncol(data03)))
{
  colnames(data03)[i]<-paste0(i,"_",colnames(data03)[i])
}

max_events=nrow(temp_mes2)
time_delays_15min=24
process_variables=c(2,11,12)
process_variables_count=length(process_variables)

smc <- data.frame(matrix(NA, ncol=process_variables_count*time_delays_15min+6, nrow=0))
columns<-c(paste0(rep(colnames(data03)[process_variables],each=time_delays_15min),"_",rep(c(1:24),process_variables_count)))
colnames(smc)[1:(time_delays_15min*process_variables_count)]<-columns
colnames(smc)[(time_delays_15min*process_variables_count+1):(time_delays_15min*process_variables_count+6)]<-colnames(data03)[c(1,13:17)]

for (n in c(10:max_events))
{
  
  smc[n,c(time_delays_15min*process_variables_count+1):(time_delays_15min*process_variables_count+6)]<-temp_mes2[n,c(1,6,7,2,3,4)]
  values<-data03[data03[,1]>=(temp_mes2[n,"czas"]-(time_delays_15min*15*60))&data03[,1]<=(temp_mes2[n,"czas"]),c(1:12)]
  values$`1_czas`<-as.numeric(values$`1_czas`-temp_mes2[n,"czas"])
  values$`1_czas`<-values$`1_czas`/(15*60)
  values$`1_czas`<-floor(values$`1_czas`)
  values<-values%>% group_by(`1_czas`) %>% summarise(across(everything(), list(mean)))
  values<-values[order(-values$`1_czas`),]
  values<-values[values$`1_czas`<0&values$`1_czas`>-25,]
  for (t in process_variables)
  {
    values[,t]<-cumsum(values[,t])
  }
  if (nrow(values)==24)
  {
  for (col_nr in c(1:process_variables_count))
  {
    smc[n,c(((col_nr-1)*time_delays_15min+1):((col_nr-1)*time_delays_15min+time_delays_15min))]<-c(t(values[,process_variables[col_nr]]))
  }
    print(n)
  }
  
}


# Save Learning dataset ----


saveRDS(smc, file = paste0("learning_dset.rds"))


