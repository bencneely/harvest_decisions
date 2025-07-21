#--------------------------------------------------------------
#Ben Neely
#07/10/2025
#Summarize decision tree data and make donut charts for each node
#--------------------------------------------------------------

## Clear R
cat("\014")  
rm(list=ls())

## Install and load packages
## Checks if package is installed, installs if not, activates for current session
if("FSA" %in% rownames(installed.packages()) == FALSE) {install.packages("FSA")}
library(FSA)

if("rio" %in% rownames(installed.packages()) == FALSE) {install.packages("rio")}
library(rio)

if("ggforce" %in% rownames(installed.packages()) == FALSE) {install.packages("ggforce")}
library(ggforce)

if("grid" %in% rownames(installed.packages()) == FALSE) {install.packages("grid")}
library(grid)

if("patchwork" %in% rownames(installed.packages()) == FALSE) {install.packages("patchwork")}
library(patchwork)

if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse")}
library(tidyverse)

## Set ggplot theme
pubtheme=theme_classic()+
  theme(panel.grid=element_blank(), 
        panel.background=element_blank(),
        plot.background=element_blank(),
        panel.border=element_rect(fill="transparent"),
        axis.title=element_text(size=22,color="black",face="bold"),
        axis.text=element_text(size=18,color="black"),
        legend.position="inside",
        legend.position.inside=c(0.01,0.99),
        legend.justification=c("left","top"),
        legend.title=element_blank())
options(scipen=999)

## Set working directory
setwd("C:/Users/Ben.Neely/OneDrive - State of Kansas, OITS/Desktop/decision trees")

## Read in data
load("mods.Rdata")

################################################################################
## Create function to extract rpart information from the models
extract_rpart_info=function(mod,model_name="unknown",id_start=1) {
  frame=mod$frame
  node_ids=as.numeric(rownames(frame))  # actual node IDs
  
  if (length(node_ids)==0 || is.null(frame$yval2)) return(NULL)
  
  total=frame$n
  yval2_all=frame$yval2
  
  harvest_counts=yval2_all[,2]
  release_counts=yval2_all[,3]
  
  pct_harvest=round(yval2_all[,4]*100,0)
  pct_release=round(yval2_all[,5]*100,0)
  
  # Get paths for each node by correct node IDs
  node_paths_list=rpart::path.rpart(mod,nodes=node_ids,print.it=FALSE)
  
  # Collapse each path into one string
  split_rules=sapply(node_paths_list, function(path_vec) {
    if (length(path_vec)==1 && path_vec=="root") {
      return("root")
    }
    path_vec=path_vec[path_vec!="root"]
    paste(path_vec,collapse=" & ")
  })
  
  data.frame(split_rule=split_rules,
             id=seq(id_start,length.out=length(node_ids)),
             total=total,
             harvest=harvest_counts,
             release=release_counts,
             p_harvest=pct_harvest,
             p_release=pct_release,
             model_name=model_name,
             stringsAsFactors=FALSE)
}

## Run function to make data ready for plotting
summary_df=do.call(rbind,lapply(seq_along(mods),function(i) {
  extract_rpart_info(mods[[i]],model_name=names(mods)[i],id_start=i)
}))

################################################################################
## Build function to make donut charts for global and each spp group
makedonut=function(harvest,release,total,title=NULL,hole_size=0.75,max_text_size=14) {
  df=data.frame(outcome=c("Harvest","Release"),
                count=c(harvest,release))
  df=df%>%
    mutate(prop=count/sum(count),
           start=c(0,head(cumsum(prop),-1))*2*pi,
           end=cumsum(prop)*2*pi)
  
  center_label=paste0(scales::percent(df$prop[1],accuracy=1)," Harvest\n",
                      scales::percent(df$prop[2],accuracy=1)," Release\n",
                      "N = ",scales::comma(total))
  
  ggplot(df)+
    geom_arc_bar(aes(x0=0,y0=0,r0=hole_size,r=1,start=start,end=end,fill=outcome),
                 color="white")+
    coord_fixed()+
    scale_fill_viridis_d(option="turbo",begin=0.1,end=0.9)+
    annotate("text",x=0,y=0,label=center_label,size=max_text_size,
             lineheight=1,hjust=0.5,vjust=0.5)+
    theme_void()+
    theme(legend.position="none",
          legend.title=element_blank())
  }

################################################################################
## Create and export global plots
## Use function to create plots
global=lapply(seq_len(nrow(subset(summary_df,model_name=="global"))), function(i) {
  row=summary_df[i, ]
  makedonut(harvest=row$harvest,
            release=row$release,
            total=row$total,
            title=paste0("Node ",row$id),
            hole_size=0.75,
            max_text_size=14)
})

## Export plots
for (i in seq_along(global)) {
  filename=paste0("global/global_", subset(summary_df,model_name=="global")$id[i],".png")
  ggsave(filename,plot=global[[i]],width=6,height=6,dpi=300,bg="transparent")
}

################################################################################
## Create and export bcf plots
## Use function to create plots
bcfdat=filter(summary_df,model_name=="bcf")
bcf=lapply(seq_len(nrow(bcfdat)), function(i) {
  row=bcfdat[i, ]
  makedonut(harvest=row$harvest,
            release=row$release,
            total=row$total,
            title=paste0("Node ",row$id),
            hole_size=0.75,
            max_text_size=14)
})

## Export plots
for (i in seq_along(bcf)) {
  filename=paste0("bcf/bcf_", subset(summary_df,model_name=="bcf")$id[i],".png")
  ggsave(filename,plot=bcf[[i]],width=6,height=6,dpi=300,bg="transparent")
}

################################################################################
## Create and export ccarp plots
## Use function to create plots
ccarpdat=filter(summary_df,model_name=="ccarp")
ccarp=lapply(seq_len(nrow(ccarpdat)), function(i) {
  row=ccarpdat[i, ]
  makedonut(harvest=row$harvest,
            release=row$release,
            total=row$total,
            title=paste0("Node ",row$id),
            hole_size=0.75,
            max_text_size=14)
})

## Export plots
for (i in seq_along(ccarp)) {
  filename=paste0("ccarp/ccarp_", subset(summary_df,model_name=="ccarp")$id[i],".png")
  ggsave(filename,plot=ccarp[[i]],width=6,height=6,dpi=300,bg="transparent")
}

################################################################################
## Create and export ccf plots
## Use function to create plots
ccfdat=filter(summary_df,model_name=="ccf")
ccf=lapply(seq_len(nrow(ccfdat)), function(i) {
  row=ccfdat[i, ]
  makedonut(harvest=row$harvest,
            release=row$release,
            total=row$total,
            title=paste0("Node ",row$id),
            hole_size=0.75,
            max_text_size=14)
})

## Export plots
for (i in seq_along(ccf)) {
  filename=paste0("ccf/ccf_", subset(summary_df,model_name=="ccf")$id[i],".png")
  ggsave(filename,plot=ccf[[i]],width=6,height=6,dpi=300,bg="transparent")
}

################################################################################
## Create and export crappie plots
## Use function to create plots
crappiedat=filter(summary_df,model_name=="crappie")
crappie=lapply(seq_len(nrow(crappiedat)), function(i) {
  row=crappiedat[i, ]
  makedonut(harvest=row$harvest,
            release=row$release,
            total=row$total,
            title=paste0("Node ",row$id),
            hole_size=0.75,
            max_text_size=14)
})

## Export plots
for (i in seq_along(crappie)) {
  filename=paste0("crappie/crappie_", subset(summary_df,model_name=="crappie")$id[i],".png")
  ggsave(filename,plot=crappie[[i]],width=6,height=6,dpi=300,bg="transparent")
}

################################################################################
## Create and export fhc plots
## Use function to create plots
fhcdat=filter(summary_df,model_name=="fhc")
fhc=lapply(seq_len(nrow(fhcdat)), function(i) {
  row=fhcdat[i, ]
  makedonut(harvest=row$harvest,
            release=row$release,
            total=row$total,
            title=paste0("Node ",row$id),
            hole_size=0.75,
            max_text_size=14)
})

## Export plots
for (i in seq_along(fhc)) {
  filename=paste0("fhc/fhc_", subset(summary_df,model_name=="fhc")$id[i],".png")
  ggsave(filename,plot=fhc[[i]],width=6,height=6,dpi=300,bg="transparent")
}

################################################################################
## Create and export fwd plots
## Use function to create plots
fwddat=filter(summary_df,model_name=="fwd")
fwd=lapply(seq_len(nrow(fwddat)), function(i) {
  row=fwddat[i, ]
  makedonut(harvest=row$harvest,
            release=row$release,
            total=row$total,
            title=paste0("Node ",row$id),
            hole_size=0.75,
            max_text_size=14)
})

## Export plots
for (i in seq_along(fwd)) {
  filename=paste0("fwd/fwd_", subset(summary_df,model_name=="fwd")$id[i],".png")
  ggsave(filename,plot=fwd[[i]],width=6,height=6,dpi=300,bg="transparent")
}

################################################################################
## Create and export lmb plots
## Use function to create plots
lmbdat=filter(summary_df,model_name=="lmb")
lmb=lapply(seq_len(nrow(lmbdat)), function(i) {
  row=lmbdat[i, ]
  makedonut(harvest=row$harvest,
            release=row$release,
            total=row$total,
            title=paste0("Node ",row$id),
            hole_size=0.75,
            max_text_size=14)
})

## Export plots
for (i in seq_along(lmb)) {
  filename=paste0("lmb/lmb_", subset(summary_df,model_name=="lmb")$id[i],".png")
  ggsave(filename,plot=lmb[[i]],width=6,height=6,dpi=300,bg="transparent")
}

################################################################################
## Create and export percid plots
## Use function to create plots
perciddat=filter(summary_df,model_name=="percid")
percid=lapply(seq_len(nrow(perciddat)), function(i) {
  row=perciddat[i, ]
  makedonut(harvest=row$harvest,
            release=row$release,
            total=row$total,
            title=paste0("Node ",row$id),
            hole_size=0.75,
            max_text_size=14)
})

## Export plots
for (i in seq_along(percid)) {
  filename=paste0("percid/percid_", subset(summary_df,model_name=="percid")$id[i],".png")
  ggsave(filename,plot=percid[[i]],width=6,height=6,dpi=300,bg="transparent")
}

################################################################################
## Create and export smb plots
## Use function to create plots
smbdat=filter(summary_df,model_name=="smb")
smb=lapply(seq_len(nrow(smbdat)), function(i) {
  row=smbdat[i, ]
  makedonut(harvest=row$harvest,
            release=row$release,
            total=row$total,
            title=paste0("Node ",row$id),
            hole_size=0.75,
            max_text_size=14)
})

## Export plots
for (i in seq_along(smb)) {
  filename=paste0("smb/smb_", subset(summary_df,model_name=="smb")$id[i],".png")
  ggsave(filename,plot=smb[[i]],width=6,height=6,dpi=300,bg="transparent")
}

################################################################################
## Create and export sunfish plots
## Use function to create plots
sunfishdat=filter(summary_df,model_name=="sunfish")
sunfish=lapply(seq_len(nrow(sunfishdat)), function(i) {
  row=sunfishdat[i, ]
  makedonut(harvest=row$harvest,
            release=row$release,
            total=row$total,
            title=paste0("Node ",row$id),
            hole_size=0.75,
            max_text_size=14)
})

## Export plots
for (i in seq_along(sunfish)) {
  filename=paste0("sunfish/sunfish_", subset(summary_df,model_name=="sunfish")$id[i],".png")
  ggsave(filename,plot=sunfish[[i]],width=6,height=6,dpi=300,bg="transparent")
}

################################################################################
## Create and export whb plots
## Use function to create plots
whbdat=filter(summary_df,model_name=="whb")
whb=lapply(seq_len(nrow(whbdat)), function(i) {
  row=whbdat[i, ]
  makedonut(harvest=row$harvest,
            release=row$release,
            total=row$total,
            title=paste0("Node ",row$id),
            hole_size=0.75,
            max_text_size=14)
})

## Export plots
for (i in seq_along(whb)) {
  filename=paste0("whb/whb_", subset(summary_df,model_name=="whb")$id[i],".png")
  ggsave(filename,plot=whb[[i]],width=6,height=6,dpi=300,bg="transparent")
}

################################################################################
## Create and export wiper plots
## Use function to create plots
wiperdat=filter(summary_df,model_name=="wiper")
wiper=lapply(seq_len(nrow(wiperdat)), function(i) {
  row=wiperdat[i, ]
  makedonut(harvest=row$harvest,
            release=row$release,
            total=row$total,
            title=paste0("Node ",row$id),
            hole_size=0.75,
            max_text_size=14)
})

## Export plots
for (i in seq_along(wiper)) {
  filename=paste0("wiper/wiper_", subset(summary_df,model_name=="wiper")$id[i],".png")
  ggsave(filename,plot=wiper[[i]],width=6,height=6,dpi=300,bg="transparent")
}