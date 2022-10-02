
library(vegan)
#step1 install.packages("vegan")
#step2 library(vegan)
library(psych)A
library(ggplot2)
library(ggrepel)
library(reshape)
library(reshape2)
library(dplyr)
library(multcompView)
library (adespatial)
library(extrafont)

font_import()
loadfonts(device = "win", quiet = TRUE)#change font to times new roman
par(tck = 0.02)
#data<- read.csv("G:/LJ/csv/Cs0723_9.csv",row.names = 1)#Storage path

A<- read.csv("G:/RDA_T.csv",row.names = 1) 
pore<- read.csv("G:/porewater.csv",row.names = 1) 

a<- A[,c(3:13)]#target variables >=2
a.log <- log1p (a)#returns the natural logarithm of (1 + number)

#env = pore[,c(3:13)]# environmental variables 
env<- pore[,c(4:6,9,11,13)]#delete vif.cca>10


env$ORP<- log10(env$ORP+1-min(env$ORP))

rda_tb <- rda(a~., env, scale = FALSE)#redundancy analysis

rda_adj <- RsquareAdj (rda_tb)$adj.r.squared#rda adjustment

screeplot(rda_tb)

#checking Variance Inflation Factors,delete factors larger than 10
vif.cca(rda_tb)#delete "EC","NA"

#forward selection
rda_tb_fs <- ordiR2step(rda(a.log~1, env, scale = FALSE),
                        scope = formula(rda_tb), 
                        R2scope = rda_adj,
                        direction = 'forward', permutations = 999)
rda_tb_fs$anova#forward selection results

rda_tb_fs.scaling2 <- summary(rda_tb_fs, scaling = 2)
rda_tb_fs.site <- data.frame(rda_tb_fs.scaling2$sites)[1:2]
rda_tb_fs.env <- data.frame(rda_tb_fs.scaling2$biplot)[1:2]
rda_tb_fs.sp <- data.frame(rda_tb_fs.scaling2$species)[1:2]

group <-  A[,c(1:2)]
group$sample <- rownames(group)

rda_tb_fs.site$sample <- rownames(rda_tb_fs.site)
rda_tb_fs.site <- merge(rda_tb_fs.site, group, by = "sample")

rda_tb_fs.env$sample <- NA
rda_tb_fs.env$group <- rownames(rda_tb_fs.env)
rda_tb_fs.sp$group <- rownames(rda_tb_fs.sp)
# constrained&unconstrained explained proportion of each RDA aixs
rda_eig <-rda_tb_fs$CCA$eig
RD1 <- round(rda_eig[1] / rda_tb_fs$tot.chi,3) * 100
RD2 <- round(rda_eig[2] / rda_tb_fs$tot.chi,3) * 100

#ggplot
library(ggplot2)
p <- ggplot(rda_tb_fs.site, aes(RDA1, RDA2)) +
  geom_point(aes(pch=pos),size=3,stroke=1.5) +
  scale_shape_manual(values = c(2,16,1,0,4,3))+
  labs(x = paste("RDA 1 (",RD1,"%)"),
       y = paste("RDA 2 (",format(RD2,nsmall =1),"%)") )+
  scale_y_continuous(breaks = seq(-1,1,0.5),limits =c(-1,1) ) +
  scale_x_continuous(breaks = seq(-1.4,1.4,0.4),limits = c(-1.41,1.41)) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title = element_blank(), 
        legend.position = c(.1, 0.82),
        legend.background  = element_rect(color=1),#border color
        legend.key = element_rect(fill = 'white'),
        legend.text = element_text(size=14),
        plot.margin = margin(t = 0.6,  # Top margin
                             r = 1,  # Right margin
                             b = 0.8,  # Bottom margin
                             l = 0.5,  # Left margin
                             unit = "cm"),
        axis.ticks = element_line(),
        axis.ticks.length=unit(-0.25, "cm"),
        axis.title.x = element_text(size = 17,color = "black",margin = margin(t = 15)),
        axis.text.x = element_text( size = 16,vjust = -0.5,color="black"),
        axis.title.y = element_text(size = 17,color = "black",margin = margin(r = 15)),
        axis.text.y = element_text( size = 16,hjust = -0.5,color = "black"),
        text=element_text(size=12,family="Times New Roman"))+
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  geom_segment(data = rda_tb_fs.env,
               aes(x = 0,y = 0, xend = RDA1,yend = RDA2), 
               arrow = arrow(length = unit(0.4, 'cm'),
                             angle=20,
                             type = "open"), 
               size = 1,
               lineend="round",
               linejoin="mitre",
               color='blue'
  ) +
  geom_text_repel(data = rda_tb_fs.env, 
            aes(RDA1*1.06, RDA2 *1.05,
                label=c("K","pH","ORP",
                        "HCO3","WTO")
            ),
            size =5,
            fontface = "bold",
            family="Times New Roman",
            position=position_jitter(width=0.14,height=0.165),
            parse = T,
            color='blue')+#environment varibles(black arrows)
  geom_segment(data = rda_tb_fs.sp,
               aes(x = 0,y = 0, xend = RDA1,yend = RDA2), 
               arrow = arrow(length = unit(0.5, 'cm'),
                             angle=20,
                             type = "closed"),
               lineend="round",
               linejoin="mitre",
               size = 1, 
               linetype=5,
               color ='red'
  ) +
  geom_text_repel(data = rda_tb_fs.sp, 
            aes(RDA1, RDA2,
                label=c("Cd","Cr","Zn","Ni",
                        "Cu","Pb","Co","Fe",
                        "Mn","As","Ti"
                        )
                ),
            size =5,
            fontface = "bold",
            family="Times New Roman",
            color='red',
            parse = T)#target variables(red arrows)
p
#display plot

