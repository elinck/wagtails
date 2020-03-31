# plots for semenov et al. in prep

library(data.table)
library(magrittr)
library(dplyr)
library(plyr)
library(zoo)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(tidyr)

# set directory
setwd("~/Dropbox/wagtails/")

# summarize mummer output
files <- list.files("data/mumer_YW_coord",full.names = T)
files <- files[grepl("\\.coords",files)]
i <- 1
for(f in files){
  if(i==1){
    dat <- fread(f,sep=" ",data.table=F)
    i=i+1
  } else {
    tmp <- fread(f,sep=" ",data.table=F)
    dat <- rbind(tmp,dat)
    i=i+1
  }
}
colnames(dat) <- c("refStart","refStop","sep1","qStart","qStop","sep2","refLength","qLength","sep3",
                   "p.identity","sep4","names")
dat$refName <- strsplit(dat$names,"\t") %>% sapply(function(e) unlist(e)[1])
dat$qName <- strsplit(dat$names,"\t") %>% sapply(function(e) unlist(e)[2])
dat <- arrange(dat,refName,refStart)
sum <- ddply(dat,
             .(qName,refName),
             summarize,
             totalMatch=sum(qLength),
             refStart=min(refStart))
sum <- arrange(sum,refName,refStart,totalMatch)
sum.a <- subset(sum,totalMatch>50000)   
sum.b <- subset(sum,totalMatch<50000)
sum.b$refName <- "NA"
sum <- rbind.data.frame(sum.a, sum.b)
sum <- ddply(sum,.(qName),function(e){                                             #get contigs with 1 hit > 1000bp
  a <- subset(e,refName==e$refName[e$totalMatch==max(totalMatch)])
  if(nrow(a)==1){
    a
  }
})

# see all chr names
unique(sum$refName)

# drop "Chr" header
sum$refName <- gsub("Chr","", sum$refName)

# chromosome order for pretty plots
chr_order <- c("1","1A","2","3","4","4A",as.character(5:15),
               as.character(17:28),"LGE22","LG34","Z","NA")

# read in fst data, sympatric pops
sym.fst <- fread("data/sympatric.parental.windowed.weir.fst")
sym.fst <- sym.fst[sym.fst$N_VARIANTS>150,]
sym.fst$comparison <- rep("sympatric", nrow(sym.fst))
sym.fst$outlier <- sym.fst$WEIGHTED_FST>=quantile(sym.fst$WEIGHTED_FST,0.995, na.rm = TRUE)
quantile(sym.fst$WEIGHTED_FST,0.995, na.rm = TRUE) #0.04650244 
sym.fst$ID <- paste0(sym.fst$CHROM,"_",sym.fst$BIN_START)

# read in fst data, allopatric pops
allo.fst <- fread("data/allopatric.parental.windowed.weir.fst")
allo.fst <- allo.fst[allo.fst$N_VARIANTS>150,]
allo.fst$comparison <- rep("allopatric", nrow(allo.fst))
allo.fst$outlier <- allo.fst$WEIGHTED_FST>=quantile(allo.fst$WEIGHTED_FST,0.995, na.rm = TRUE)
quantile(allo.fst$WEIGHTED_FST,0.995, na.rm = TRUE) #0.3600221
allo.fst$ID <- paste0(allo.fst$CHROM,"_",allo.fst$BIN_START)

# select only shared windows
id.set <- intersect(sym.fst$ID, allo.fst$ID)
sym.fst <- sym.fst[which(sym.fst$ID %in% id.set),]
allo.fst <- allo.fst[which(allo.fst$ID %in% id.set),]
sym.fst <- arrange(sym.fst,CHROM,BIN_START)
allo.fst <- arrange(allo.fst,CHROM,BIN_START)

# set row numbers, rolling means
sym.fst$rollmean <- rollmean(sym.fst$WEIGHTED_FST,150,na.pad = TRUE)
allo.fst$rollmean <- rollmean(allo.fst$WEIGHTED_FST,150,na.pad = TRUE)

# merge dataframes
win.df <- rbind.data.frame(sym.fst,allo.fst)
colnames(win.df) <- c("scaffold","start","end","n_variants","weighted_fst",
                      "mean_fst","comparison","outlier","unique_ID","rollmean")

#merge mummer info with windowed stats
win.df <- merge(win.df,sum,by.x="scaffold",by.y="qName",all.x=T,all.y=F)
win.df$refName[!win.df$refName %in% chr_order] <- "NA"
win.df$chr <- factor(win.df$refName,levels=chr_order)
win.df <- arrange(win.df,chr,refStart)
sym.tmp <- win.df[win.df$comparison=="sympatric",]
sym.tmp$row <- 1:nrow(sym.tmp)
allo.tmp <- win.df[win.df$comparison=="allopatric",]
allo.tmp$row <- 1:nrow(allo.tmp)
win.df <- rbind.data.frame(sym.tmp,allo.tmp)

# chr labels dataframe
chr_labels <- ddply(win.df,.(chr),summarize,mid=median(row),start=min(row),stop=max(row))
chr_labels$chr <- as.character(chr_labels$chr)
chr_labels$chr[chr_labels$chr %in% as.character(12:19)] <- "12-19"
chr_labels$mid[chr_labels$chr=="12-19"] <- median(chr_labels$mid[chr_labels$chr=="12-19"],na.rm=T)
chr_labels$start[chr_labels$chr=="12-19"] <- min(chr_labels$start[chr_labels$chr=="12-19"],na.rm=T)
chr_labels$stop[chr_labels$chr=="12-19"] <- max(chr_labels$stop[chr_labels$chr=="12-19"],na.rm=T)
chr_labels$chr[chr_labels$chr %in% as.character(21:28)] <- "21-28"
chr_labels$mid[chr_labels$chr=="21-28"] <- median(chr_labels$mid[chr_labels$chr=="21-28"],na.rm=T)
chr_labels$start[chr_labels$chr=="21-28"] <- min(chr_labels$start[chr_labels$chr=="21-28"],na.rm=T)
chr_labels$stop[chr_labels$chr=="21-28"] <- max(chr_labels$stop[chr_labels$chr=="21-28"],na.rm=T)
chr_labels$chr[chr_labels$chr=="LGE22"] <- ""
chr_labels$mid[chr_labels$chr==""] <- median(chr_labels$mid[chr_labels$chr==""],na.rm=T)
chr_labels$start[chr_labels$chr==""] <- min(chr_labels$start[chr_labels$chr==""],na.rm=T)
chr_labels$stop[chr_labels$chr==""] <- max(chr_labels$stop[chr_labels$chr==""],na.rm=T)
chr_labels$start[chr_labels$chr=="1B"] <- 0
chr_labels$stop[chr_labels$chr=="1B"] <- 0
chr_labels <- subset(chr_labels,!is.na(chr) & !duplicated(chr))

# plot!
p1 <- ggplot(data=win.df,aes(x=row,y=weighted_fst,col=chr))+
  theme_bw()+
  facet_grid(comparison~.)+
  theme(text=element_text(size=10),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        axis.title.x=element_blank(),
        strip.text=element_text(size=10),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position="none")+
  scale_y_continuous(breaks=c(0,.25,.5,.75,1))+
  scale_color_manual(values=rep(c("#2b2a27","#a6a5a2"),length(unique(win.df$chr))/2+1))+
  geom_point(data=subset(win.df, weighted_fst>0), size=1.1,shape=21)+
  #geom_point(data=subset(win.df, outlier="TRUE"),size=1.2,shape=21,col="red")+
  #geom_line(aes(y=rollmean),lwd=0.5,col="black")+
  #geom_segment(data=chr_segments,aes(x=start+50,xend=stop-50,y=Fst,yend=Fst,col=NA),col="black")+
  #geom_point(data=subset(win.df, outlier="TRUE"),size=0.4,shape=21,stroke=0.4,col="red")+
  geom_text(data=chr_labels,aes(label=chr,x=mid,y=-0.1,col=NA),
            col="black",size=2,angle=0) +  
  geom_vline(xintercept = chr_labels[chr_labels$chr=="20",]$start, linetype="solid", 
             color = "yellow", size=3, alpha=0.3) +
  geom_vline(xintercept = 13750.0, linetype="solid", 
             color = "blue", size=3, alpha=0.3) +
  geom_hline(aes(yintercept=0)) +
  geom_hline(data=subset(win.df, comparison=="sympatric"),aes(yintercept=0.04650244),linetype="dashed") +
  geom_hline(data=subset(win.df, comparison=="allopatric"),aes(yintercept=0.3600221),linetype="dashed") +
  labs(y=expression(F[ST])) 

# read in gemma results
gemma.tmp <- fread("data/face.assoc.txt")

# merge w/ mummer data
gemma.df <- merge(gemma.tmp,sum,by.x="chr",by.y="qName",all.x=T,all.y=F)
gemma.df$chrom <- factor(gemma.df$refName,levels=chr_order)
gemma.df <- arrange(gemma.df,chrom,refStart)
gemma.df$p_wald <- log10(gemma.df$p_wald)
gemma.df$row <- 1:nrow(gemma.df)
gemma.sub <- sample_n(gemma.df, 1000000)
gemma.sub <- arrange(gemma.sub,chrom,refStart)
gemma.sub$label <- 'GEMMA'

# gemma chr labels 
chr_labels_2 <- ddply(gemma.sub,.(chrom),summarize,mid=median(row),start=min(row),stop=max(row))
chr_labels_2$chrom <- as.character(chr_labels_2$chrom)
chr_labels_2$chrom[chr_labels_2$chrom %in% as.character(12:19)] <- "12-19"
chr_labels_2$mid[chr_labels_2$chrom=="12-19"] <- median(chr_labels_2$mid[chr_labels_2$chrom=="12-19"],na.rm=T)
chr_labels_2$start[chr_labels_2$chrom=="12-19"] <- min(chr_labels_2$start[chr_labels_2$chrom=="12-19"],na.rm=T)
chr_labels_2$stop[chr_labels_2$chrom=="12-19"] <- max(chr_labels_2$stop[chr_labels_2$chrom=="12-19"],na.rm=T)
chr_labels_2$chrom[chr_labels_2$chrom %in% as.character(21:28)] <- "21-28"
chr_labels_2$mid[chr_labels_2$chrom=="21-28"] <- median(chr_labels_2$mid[chr_labels_2$chrom=="21-28"],na.rm=T)
chr_labels_2$start[chr_labels_2$chrom=="21-28"] <- min(chr_labels_2$start[chr_labels_2$chrom=="21-28"],na.rm=T)
chr_labels_2$stop[chr_labels_2$chrom=="21-28"] <- max(chr_labels_2$stop[chr_labels_2$chrom=="21-28"],na.rm=T)
chr_labels_2$chrom[chr_labels_2$chrom=="LGE22"] <- ""
chr_labels_2$mid[chr_labels_2$chrom==""] <- median(chr_labels_2$mid[chr_labels_2$chrom==""],na.rm=T)
chr_labels_2$start[chr_labels_2$chrom==""] <- min(chr_labels_2$start[chr_labels_2$chrom==""],na.rm=T)
chr_labels_2$stop[chr_labels_2$chrom==""] <- max(chr_labels_2$stop[chr_labels_2$chrom==""],na.rm=T)
chr_labels_2$start[chr_labels_2$chrom=="1B"] <- 0
chr_labels_2$stop[chr_labels_2$chrom=="1B"] <- 0
chr_labels_2 <- subset(chr_labels_2,!is.na(chrom) & !duplicated(chrom))

# get outlier cutoff
quantile(gemma.sub$p_wald,0.005, na.rm = TRUE) #-2.253437 

# plot gemma
p2 <- ggplot(data=gemma.sub,aes(x=row,y=p_wald,col=chrom))+
  facet_grid(label~.)+
  theme_bw()+
  theme(text=element_text(size=10),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        axis.title.x=element_blank(),
        strip.text=element_text(size=10),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position="none")+
  geom_point(size=1.1,shape=21) +
  scale_color_manual(values=rep(c("#2b2a27","#a6a5a2"),length(unique(gemma.sub$chrom))/2+1))+
  #geom_line(aes(y=rollmean),lwd=0.5,col="black")+
  geom_text(data=chr_labels_2,aes(label=chrom,x=mid,y=1,col=NA),
            col="black",size=2,angle=0) +  
  geom_vline(xintercept = chr_labels_2[chr_labels_2$chrom=="20",]$start, linetype="solid", 
             color = "yellow", size=3, alpha=0.3) +
  geom_vline(xintercept = 1345000, linetype="solid", 
             color = "blue", size=3, alpha=0.3) +
  geom_hline(aes(yintercept=0)) +
  geom_hline(aes(yintercept=-2.25479),linetype="dashed") +
  labs(y=expression(P[Wald])) 

png(file="figures/manhattan_gemma.png",res=300,width=8.5,height=6.5,units="in")
plot_grid(p1,p2,ncol=1,rel_heights = c(2,1))  
dev.off()

# load likelihood data
im <- read.table("~/Dropbox/wagtails/data/IM_realparams.txt") %>% as.data.frame()
sc <- read.table("~/Dropbox/wagtails/data/SC_realparams.txt") %>% as.data.frame()
si <- read.table("~/Dropbox/wagtails/data/SI_realparams.txt") %>% as.data.frame()
colnames(im) <- c("nPer","nAlb","tSplit","m12","m21","ll_model","theta")
colnames(sc) <- c("nPer","nAlb","t1","t2","m12","m21","ll_model","theta")
colnames(si) <- c("nPer","nAlb","tSplit","ll_model")
im.df <- cbind.data.frame(im$ll_model, rep("IM",nrow(im)))
colnames(im.df) <- c("log_likelihood","model")
im.df$AIC <- (-2*im.df$log_likelihood) + (2*6)
sc.df <- cbind.data.frame(sc$ll_model, rep("SC",nrow(sc)))
colnames(sc.df) <- c("log_likelihood","model")
sc.df$AIC <- (-2*sc.df$log_likelihood) + (2*7)
si.df <- cbind.data.frame(si$ll_model, rep("SI",nrow(si)))
colnames(si.df) <- c("log_likelihood","model")
si.df$AIC <- (-2*si.df$log_likelihood) + (2*4)
mod.df <- rbind.data.frame(im.df,sc.df,si.df)

# get best supported model
best_aic <- mod.df[which.min(mod.df$AIC),]
best_ll <- mod.df[which.max(mod.df$log_likelihood),]

# get optimization median
mod_dat_aic <- ddply(mod.df, "model", summarise, median_param=median(AIC))
mod_dat_ll <- ddply(mod.df, "model", summarise, median_param=median(log_likelihood))


# plot model comparison
p3 <- ggplot(mod.df, aes(x=model,y=log_likelihood)) +
  geom_boxplot() +
  geom_jitter(pch=21) +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  ylab("Log-likelihood") +
  xlab("Model")

# load param bootstraps
sc_boots <- read.table("~/Dropbox/wagtails/data/SC_realparams_boots.txt") %>% as.data.frame()
colnames(sc_boots) <- c("nPer","nAlb","t1","t2","m12","m21","ll_model","theta")
sc_boots <- sc_boots[,1:6]
colnames(sc_boots) <- c(expression("N[e~italic(~personata)]"), expression("N[e~italic(~alba)]"), 
                        "Duration~of~isolation~(years)", "Duration~of~hybridization~(years)",
                        expression("N[m~italic(~alba-personata)]"),expression("N[m~italic(~personata-alba)]"))
sc_boots.df <- gather(sc_boots)
colnames(sc_boots.df) <- c("parameter","value")

# get mean param center
pdat <- ddply(sc_boots.df, "parameter", summarise, median_param=median(value))
pdat_sd <- ddply(sc_boots.df, "parameter", summarise, sd_param=sd(value))

  
p4 <- ggplot(sc_boots.df, aes(x=value)) +
  #geom_histogram() +
  geom_density(alpha=0.2,color="black",fill="black") +
  facet_wrap(~ parameter, scales="free",labeller=label_parsed) +
  theme_classic() +
  ylab("Density") +
  xlab(element_blank()) +
  scale_x_continuous(breaks = scales::pretty_breaks(4), limits = c(0, NA)) +
  geom_vline(data=pdat, aes(xintercept=median_param),
             linetype="solid", size=1) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size=9),
        axis.text = element_text(size=9))

pdf(file="figures/demographic_inference.pdf",width=12,height=4.5)
plot_grid(p3,p4,labels="AUTO",ncol=2,rel_widths = c(1,1.75)) 
dev.off()

# read hybrid zone simulation results
cline_df <- read.csv("~/Dropbox/wagtails/data/simulation_cline_parameters.csv")

# extract final generation
cline_500 <- cline_df[cline_df$generation==500,]

# trajectory plot
p5 <- ggplot(cline_df,aes(x=generation,y=cline_center,col=cline_type,
                          group=interaction(replicate, cline_type)))+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        strip.background = element_blank(),
        strip.text=element_text(size=12))+
  scale_color_manual(values=c("red3","steelblue3"))+
  xlab("Generations")+ylab("Transect Location")+
  geom_line(alpha=0.7,size=0.75) +
  ylim(0,1) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid = element_blank())

# get mean cline center
cdat <- ddply(cline_500, "cline_type", summarise, mean_center=mean(cline_center))

# plot distribution
p6 <- ggplot(cline_500,aes(x=cline_center,fill=cline_type)) +
  theme_bw()+
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c("red3","steelblue3")) +
  geom_vline(data=cdat, aes(xintercept=mean_center,  color=cline_type),
             linetype="dashed", size=1) +
  scale_color_manual(values=c("red3","steelblue3")) +
  xlim(0,1) +
  coord_flip() +
  ylab("Density") +
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text = element_text(size=12),
        axis.title.x = element_text(size=12),
        legend.title = element_blank(),
        panel.grid = element_blank())

# load plotting function from hzam
plot_freqs <-  function(format_output, fit_neutral=FALSE, fit_trait=FALSE){
  
  # base plot
  p <- ggplot(format_output, aes(x=location, y=freqs, color=vars)) +
    geom_point(size=4,alpha=0.5,pch=21) +
    theme_bw() +
    scale_color_manual(values=c("steelblue3","red3"),
                       breaks=c("neutral","traits"),
                       labels=c("Neutral loci","Mating trait loci")) +
    ylab("Hybrid Index") +
    xlab("Transect Location") +
    #coord_flip() +
    theme(axis.title = element_text(size=12),
          axis.text = element_text(size=12),
          legend.position = "bottom",
          legend.title = element_blank(),
          panel.grid = element_blank())
  
  if(fit_neutral==T){
    model_n <- gam(freqs ~ s(location), data=subset(format_output,vars=="neutral"), quasibinomial(link = "logit"),
                   method = "P-ML")
    model_np <- predict_gam(model_n, length_out = 1000)
    model_np$fit <- inv.logit(model_np$fit)
    colnames(model_np) <- c("location", "freqs", "se.fit")
    p <- p + geom_smooth(data = model_np, inherit.aes = FALSE, 
                         mapping=aes(location, freqs),linetype="solid",col="black",se=FALSE) +
      labs(linetype="Cline")
    
  }
  if(fit_trait==T){
    model_t <- gam(freqs ~ s(location), data=subset(format_output,vars=="traits"), quasibinomial(link = "logit"),
                   method = "P-ML")
    model_tp <- predict_gam(model_t, length_out = 10000)
    model_tp$fit <- inv.logit(model_tp$fit)
    colnames(model_tp) <- c("location", "freqs", "se.fit")
    p <- p + geom_smooth(data = model_tp, inherit.aes = FALSE, 
                         mapping=aes(location, freqs),linetype="dashed",col="black",se=FALSE) +
      labs(linetype="Cline")
  }
  return(p)
}

pop_data <- read.csv("~/Dropbox/wagtails/data/displaced_cline_example.csv")


# plot hybrid zone
p7 <- plot_freqs(pop_data, fit_neutral = TRUE, fit_trait = TRUE)

pdf(file="figures/simulation_results.pdf",width=15,height=6)
plot_grid(p5,p6,p7, labels="AUTO",ncol=3)
dev.off()



