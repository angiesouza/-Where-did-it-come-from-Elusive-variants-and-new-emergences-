
library(ggplot2)
library(tidyverse)
library(viridis)
library(gridExtra)
library(dbplyr)

# Figure 2 - Probability of establishment ####

dados <-read.table("Output_prob_stablishment.txt", header=F, dec=".")
colnames(dados) <- c("Intro", "Tx_rep", "Tam_Prop", "Distancia" , "PropIntro", "Sobreviveu","Num_sobrev")

bs=5


trep= 15 
nrep=200
ps=c(1,10,100)

# Function

xi <- seq(0,4,0.1)
yi <- exp((-xi**2)/2)
yii <- rep(c(yi), 3)
z <- rep("Function", 123)

#

tint = 1

tamanho <- mutate(dados,dis_int = Distancia*10)

Tam_Prop3 <- filter(tamanho, Intro == tint)



grupo <- data.frame(TamProp=numeric()
                    ,Distancia=numeric()
                    ,Probabilidade1=numeric()
                    ,Probabilidade200=numeric()
                    ,Funcao=numeric())

nc=1

for(i in 1:3){
  filter1 <- filter(Tam_Prop3, Tam_Prop == ps[i])
  
  for (n in 0:40) {
    filter2 <- filter(filter1, dis_int == n)
    prob <- mean(filter2$Sobreviveu)
    
    grupo[nc,]=NA
    grupo$TamProp[nc] <- ps[i]
    grupo$Distancia[nc] <- n
    grupo$Probabilidade1[nc] <- prob
    
    nc=nc+1
  }}

#

tint=200

Tam_Prop3 <- filter(tamanho, Intro == tint)


nc=1

for(i in 1:3){
  filter1 <- filter(Tam_Prop3, Tam_Prop == ps[i])
  
  for (n in 0:40) {
    filter2 <- filter(filter1, dis_int == n)
    prob <- mean(filter2$Sobreviveu)
    grupo$Probabilidade200[nc] <- prob
    nc=nc+1
  }}

grupo$Funcao <- yii
grupo$gFuncao <- z

note <- expression(paste("Distance between optimal phenotypes imposed by hosts (", italic("D"), ")"))

Fig2 <- ggplot (data=grupo) +
  geom_line( aes(x=Distancia/10, y=Funcao,group=gFuncao,color=factor(gFuncao))) +
  geom_line( aes(x=Distancia/10, y=Probabilidade1, group=TamProp,color=factor(TamProp)), linetype="dashed")+
  geom_line( aes(x=Distancia/10, y=Probabilidade200, group=TamProp,color=factor(TamProp))) +
  labs(x= note, y="Probability of\n establishment" , colour="Propagule size")+
  scale_color_manual(values =  c("#481567FF", "#20A387FF", "#FDE725FF","grey")) +
  theme_set(theme_bw(base_size = bs))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

ggsave(paste0("Figure2.png"), width = 3, height = 1.5, dpi = 1000)
grid.arrange(Fig2, ncol=1)
dev.off() 

# ####

# Figure 3 - Histogram ####

df1 <-read.table("phenotypic_variability.txt", header=T, dec=".")
df2 <-read.table("propagule_id.txt", header=T, dec=".")

ind=100

df1a <- filter(df1, individuos==ind & distancia<=3.0)
df1a <- filter(df1a, host == 2)
df1a <- mutate(df1a, sobreviveu=sucesso)
df1a <- mutate(df1a, nind=individuos)
df1a <- select(df1a, nind, distancia, repeticao, tempo, id_linhagem, fenotipo, sobreviveu)

df2a <- mutate(df2, distancia = distancia/10)
df2a <- filter(df2a, nind==ind & distancia<=3.0)
df2a <- mutate(df2a, tx_rep = NULL)

df3a <- filter(df2a, sobreviveu==1  & distancia<=3.0)

df1a <- mutate(df1a, Status="Donor population")
df2a <- mutate(df2a, Status="Propagule")
df3a <- mutate(df3a, Status="Survivals")

df <- rbind(df1a,df2a, df3a)


Fig3 <- ggplot(df, aes(x = fenotipo, fill=Status)) +
  geom_histogram(binwidth=1, position="identity")  +
  scale_x_continuous(limits=c(480,540))+
  scale_fill_viridis_d()+
  facet_grid(distancia~.)+
  geom_vline(data=filter(df3a), aes(xintercept=500), colour="yellow", linetype="dotted",   lwd=0.35) +
  geom_vline(data=filter(df3a, distancia=="0"), aes(xintercept=500), colour="black", linetype="dotted", lwd=0.35) +
  geom_vline(data=filter(df3a, distancia=="0.5"), aes(xintercept=505), colour="black", linetype="dotted", lwd=0.35) + 
  geom_vline(data=filter(df3a, distancia=="1"), aes(xintercept=510), colour="black", linetype="dotted", lwd=0.35) + 
  geom_vline(data=filter(df3a, distancia=="1.5"), aes(xintercept=515), colour="black", linetype="dotted", lwd=0.35) + 
  geom_vline(data=filter(df3a, distancia=="2"), aes(xintercept=520), colour="black", linetype="dotted", lwd=0.35) + 
  geom_vline(data=filter(df3a, distancia=="2.5"), aes(xintercept=525), colour="black", linetype="dotted", lwd=0.35) + 
  geom_vline(data=filter(df3a, distancia=="3"), aes(xintercept=530), colour="black", linetype="dotted", lwd=0.35) + 
  geom_vline(data=subset(df3a, distancia=="0"), aes(xintercept=mean(df3a$fenotipo[df3a$distancia=="0"])), colour="yellow",  lwd=0.35) +
  geom_vline(data=subset(df3a, distancia=="0.5"), aes(xintercept=mean(df3a$fenotipo[df3a$distancia=="0.5"])), colour="yellow",  lwd=0.35) +
  geom_vline(data=subset(df3a, distancia=="1"), aes(xintercept=mean(df3a$fenotipo[df3a$distancia=="1"])), colour="yellow",  lwd=0.35) +
  geom_vline(data=subset(df3a, distancia=="1.5"), aes(xintercept=mean(df3a$fenotipo[df3a$distancia=="1.5"])), colour="yellow",  lwd=0.35) +
  geom_vline(data=subset(df3a, distancia=="2"), aes(xintercept=mean(df3a$fenotipo[df3a$distancia=="2"])), colour="yellow",  lwd=0.35) +
  geom_vline(data=subset(df3a, distancia=="2.5"), aes(xintercept=mean(df3a$fenotipo[df3a$distancia=="2.5"])), colour="yellow", lwd=0.35) +
  geom_vline(data=subset(df3a, distancia=="3"), aes(xintercept=mean(df3a$fenotipo[df3a$distancia=="3"])), colour="yellow",  lwd=0.35) +
  labs(x="Phenotypes", y="Frequency")+
  theme_set(theme_bw(base_size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(), legend.position="bottom") +
  guides(size = "none")

ggsave(paste0("Figure3.png"), width = 2.3, height = 3.6)
grid.arrange(Fig3, ncol=1)
dev.off() 

# ####

# Figure 4 - Boxplot + line graphic ####

# boxplot

df_fen<-read.table("phenotypic_variability.txt", header=T, dec=".")
names(df_fen)[7] <- "Host"
df_fen$Host <- as.factor(df_fen$Host)
levels( df_fen$Host)[levels( df_fen$Host)=="1"] <- "New population"
levels( df_fen$Host)[levels( df_fen$Host)=="2"] <- "Donor population"

filter1 <- filter(df_fen, tempo == 1000 & sucesso == 1 & distancia < 3.5 )
filter1$distancia <- as.factor(filter1$distancia)

#cores=c( "#414487FF","#440154FF")
cores=c( "#365D8DFF","#440154FF")


A <- ggplot(data = filter1, aes(x =  distancia , y =  fenotipo, fill = Host)) +
  geom_boxplot( size=0.18, outlier.size=0.001, alpha=0.9) +
  coord_cartesian(ylim=c(485,542))+
  scale_fill_manual(values=cores)+
  theme_set(theme_bw(base_size = 6))+
  #   geom_vline(aes(xintercept= 7.5), colour="black",linetype= "dashed", lwd=0.15)+
  labs (x= "Optimal phenotype distance between hosts", y= "Phenotype",tag = "A")+
  #scale_shape(name = "Hospedeiro", labels = c("New", "Old"))+
  theme(axis.text.x = element_text( ),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position="bottom")


# line graph

df_fen<-read.table("phenotypic_variability.txt", header=T, dec=".")
df2 <- filter(df_fen, tempo == 1000)
df2 <- mutate(df2,dis_int = distancia*10)

resumo=function(data,ic){
  distancias=unique(data$Distancia)
  data_sum=data.frame(Distancia=as.numeric()
                      ,Fenotipos_md=numeric()
                      ,Fenotipos_sup=numeric()
                      ,Fenotipos_inf=numeric()
                      ,Amplitude_md=numeric()
                      ,Amplitude_sup=numeric()
                      ,Amplitude_inf=numeric()
                      ,Linhagens_md=numeric()
                      ,Linhagens_sup=numeric()
                      ,Linhagens_inf=numeric()
  )

  
  nl=0
  for(n in distancias){
    sub=data[data$Distancia==n,]
    replicas=nrow(sub)
    nl=nl+1
    data_sum[nl,]<-NA
    data_sum$Distancia[nl]<-n
    data_sum$Fenotipos_md[nl]<-mean(sub$Fenotipos)
    data_sum$Amplitude_md[nl]<-mean(sub$Amplitude)
    data_sum$Linhagens_md[nl]<-mean(sub$Linhagens)
    
    rank_fen <- sort(sub$Fenotipos)
    rank_amp <- sort(sub$Amplitude)
    rank_lin <- sort(sub$Linhagens)
    
    inf=1+round(((1-ic)/2)*replicas,0)
    sup=replicas-round(((1-ic)/2)*replicas,0)
    
    data_sum$Fenotipos_sup[nl]<-rank_fen[sup]
    data_sum$Amplitude_sup[nl]<-rank_amp[sup]
    data_sum$Linhagens_sup[nl]<-rank_lin[sup]
    
    data_sum$Fenotipos_inf[nl]<-rank_fen[inf]
    data_sum$Amplitude_inf[nl]<-rank_amp[inf]
    data_sum$Linhagens_inf[nl]<-rank_lin[inf]
    
  }
  return(data_sum)
}


nc=1

grupo <- data.frame(Distancia=numeric()
                    ,Fenotipos=numeric()
                    ,Amplitude=numeric()
                    ,Linhagens=numeric())

for (n in 0:30) {
  filter1 <- filter(df2, dis_int == n & sucesso == 1)
  
  nrep=unique(filter1$repeticao)
  
  if(length(nrep)>9) {
    
    soma_ampl=0
    soma_fen=0
    soma_li=0
    
    for (r in nrep){
      cond <- filter1$repeticao == r
      am <- max(filter1$fenotipo[cond])-min(filter1$fenotipo[cond])
      fe <- length(unique(filter1$fenotipo[cond]))
      li <- length(unique(filter1$id_linhagem[cond]))
      
      grupo[nc,]=NA
      grupo$Distancia[nc] <- n/10
      grupo$Fenotipos[nc] <- fe
      grupo$Amplitude[nc] <- am
      grupo$Linhagens[nc] <- li
      
      nc=nc+1    
      
    }}}

grupo <- na.omit(grupo)

data <- resumo(grupo,0.95)


B <- ggplot(data, aes(Distancia)) + 
  geom_line(aes(y = Fenotipos_md, colour = "Number of phenotypes"))+
  geom_line(aes(y = Amplitude_md, colour = "Phenotypic amplitude"))+
  geom_ribbon(aes(ymin=Fenotipos_inf, ymax=Fenotipos_sup,fill = "Number of phenotypes"), alpha=0.2,colour = NA)+
  geom_ribbon(aes(ymin=Amplitude_inf, ymax=Amplitude_sup,fill = "Phenotypic amplitude"), alpha=0.2,colour = NA)+
  scale_color_viridis_d(name= " ")+
  scale_fill_viridis_d(name= " ")+
  labs(x="Optimal phenotype distance between hosts", y=" ",tag = "B")+
  theme_set(theme_bw(base_size = 6))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position="bottom")


teste <- grid.arrange(A, B, ncol=2, nrow=1)
ggsave(paste0("Figure4.png"), teste, width = 4.8, height = 2)
dev.off() 

# ####

# Figure 5 - rare events ####

df1 <-read.table("phenotypic_variability.txt", header=T, dec=".")
df2 <-read.table("propagule_id.txt", header=T, dec=".")

restricao <- unique(df2$repeticao)

df1a <- df1[df1$repeticao%in%restricao,]
df1a <- mutate(df1a, Status="Donor population")
df1a$Status[df1a$host == "1"] <-"New population"
df1a <- select(df1a,nind, distancia, repeticao, tempo, id_linhagem, fenotipo, sobreviveu, Status)

df2a <- mutate(df2, distancia = distancia/10, tempo = 200)
df3a <- filter(df2a, sobreviveu==1 )

df2a <- mutate(df2a, Status="Propagule")
df3a <- mutate(df3a, Status="Survivals")

df <- rbind(df1a,df2a,df3a)

rest <- c(200, 205, 210, 220,250, 1000)

df<- df[df$tempo%in%rest,]

cores=c( "#440154FF","#365D8DFF","#1F968BFF", "#FDE725FF")

plotraras <- ggplot(df, aes(x = fenotipo, fill=Status)) +
  geom_histogram(binwidth=1, position="identity")  +
  scale_x_continuous(limits=c(490,545))+
  scale_fill_manual(values=cores)+
  facet_grid(tempo~.)+
  geom_vline(data=filter(df), aes(xintercept=500), colour="#FDE725FF", linetype="dotted",   lwd=0.35) +
  geom_vline(data=filter(df), aes(xintercept=535), colour="black", linetype="dotted",   lwd=0.35) +
  geom_vline(data=subset(df, distancia=="3.5"), aes(xintercept=mean(df3a$fenotipo[df3a$distancia=="3.5"])), colour="#FDE725FF", lwd=0.35) +
  labs(x="Phenotypes", y="Frequency")+
  theme_set(theme_bw(base_size = 6))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.position="bottom")+
  guides(fill=guide_legend(ncol=2))


ggsave(paste0("Figure5.png"), width = 2, height = 3)
grid.arrange(plotraras, ncol=1)
dev.off() 

# ####
