library(data.table)
library(magrittr)
library(dplyr)
library(plyr)
library(zoo)
library(ggplot2)
library(ggrepel)
library(cowplot)

# set directory
setwd("~/Dropbox/wagtails/")

# summarize mummer output
files <- list.files("mumer_YW_coord",full.names = T)
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
sum.a <- subset(sum,totalMatch>20000)   
sum.b <- subset(sum,totalMatch<20000)
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
sym.fst <- fread("sympatric.parental.windowed.weir.fst")
sym.fst <- sym.fst[sym.fst$N_VARIANTS>150,]
sym.fst$comparison <- rep("sympatric", nrow(sym.fst))
sym.fst$outlier <- sym.fst$WEIGHTED_FST>=quantile(sym.fst$WEIGHTED_FST,0.995, na.rm = TRUE)
quantile(sym.fst$WEIGHTED_FST,0.995, na.rm = TRUE) #0.04650244 
sym.fst$ID <- paste0(sym.fst$CHROM,"_",sym.fst$BIN_START)

# read in fst data, allopatric pops
allo.fst <- fread("allopatric.parental.windowed.weir.fst")
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
  geom_point(data=subset(win.df, weighted_fst>0), size=1,shape=21,alpha=0.8)+
  geom_line(aes(y=rollmean),lwd=0.5,col="black")+
  geom_text(data=chr_labels,aes(label=chr,x=mid,y=-0.1,col=NA),
            col="black",size=2,angle=0) +  
  geom_vline(xintercept = chr_labels[chr_labels$chr=="20",]$start, linetype="solid", 
             color = "yellow", size=3, alpha=0.3) +
  geom_hline(aes(yintercept=0)) +
  geom_hline(data=subset(win.df, comparison=="sympatric"),aes(yintercept=0.04650244),linetype="dashed") +
  geom_hline(data=subset(win.df, comparison=="allopatric"),aes(yintercept=0.3600221),linetype="dashed") +
  labs(y=expression(F[ST])) 

png("fst_test.png",width = 9, res = 300, height = 4, units = "in")
p1
dev.off()


# read in gemma results
gemma.tmp <- fread("face.assoc.txt")

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
quantile(gemma.sub$p_wald,0.005, na.rm = TRUE) #0.3600221

# calculate rolling mean
gemma.sub$rollmean <- rollmean(gemma.sub$p_wald,150,na.pad = TRUE)

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
  geom_point(size=1,shape=21,alpha=0.8) +
  scale_color_manual(values=rep(c("#2b2a27","#a6a5a2"),length(unique(gemma.sub$chrom))/2+1))+
  geom_line(aes(y=rollmean),lwd=0.5,col="black")+
  geom_text(data=chr_labels_2,aes(label=chrom,x=mid,y=1,col=NA),
            col="black",size=2,angle=0) +  
  geom_vline(xintercept = chr_labels_2[chr_labels_2$chrom=="20",]$start, linetype="solid", 
             color = "yellow", size=3, alpha=0.3) +
  geom_hline(aes(yintercept=0)) +
  geom_hline(aes(yintercept=-2.25479),linetype="dashed") +
  labs(y=expression(P[Wald])) 

# write to png together
png(file="manhattan_gemma.png",res=300,width=8,height=6,units="in")
plot_grid(p1,p2,ncol=1,rel_heights = c(2,1))  
dev.off()


