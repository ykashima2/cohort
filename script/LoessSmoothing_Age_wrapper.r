library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggridges)
library(tibble)
library(patchwork)
library(ggpmisc)
library(ggrepel)
set.seed(1234)

########################################################################
# Parameter
########################################################################

Type_list <- c("Bcell","NKcell", "Mono" ,"Tcell" )

param <- c()
param$ATAC_Filter1 = 4        
param$ATAC_Open_th = 0.1      
param$RNA_Filter1  = 4        
param$RNA_filt_cellrate = 0.1  

param$Ass_ATAC_delt = 0.1
param$Ass_RNAmed    = 1.5
param$Ass_RNAfc     = 2

#--------
INPUTDIR <- "data/"
outdir = "output/"
DATAdir = "Input/celltype/"
SampleListfile <- paste0( INPUTDIR , "Sample_inf.txt")

# OUTPUT
Pearson_dir <- paste0( outdir , "Figure3/")
OUTdir1 <- paste0( Pearson_dir , "SummaryB1/")
OUTdir2 <- paste0( Pearson_dir , "SummaryB2/")

## create Dir
for(dir_path in c( Pearson_dir , OUTdir1 , OUTdir2)){
        if (!dir.exists(dir_path)) { #ディレクトリが存在するか確認
                dir.create(dir_path) #ディレクトリが存在しない場合、作成
                print(paste("ディレクトリ", dir_path, "を作成しました"))
        }
}


# OUTPUT
metaInf <- read.table(file = SampleListfile, header=T,stringsAsFactors=F,sep ="\t")
metaInf <- metaInf  %>% filter( type != "2nd")
metaInf <- metaInf %>%
                mutate(groupe=if_else(type == "1st" ,"old1",
                                        if_else(age >=55 , "old2","young"  )  ) )

source( "script/1_LoessSmoothing.data_util.r" )
metaInfUse <- metaInf %>% filter( groupe == "old1" )
sample_set <-metaInfUse %>% pull(Sample)
forPlot_GeneStats_res <- doall(  metaInf = metaInfUse , sample_set ,DATAdir ,Type_list , OUTdir1 ,param )

param2 <- param
param2$Ass_RNAfc <- 1.5
param2$Ass_ATAC_delt <- 0.05

metaInfUse <- metaInf %>% filter( groupe == "old1"| groupe == "old2" )
sample_set <-metaInfUse %>% pull(Sample)
forPlot_GeneStats_res_B2 <-doall(  metaInf = metaInfUse , sample_set ,DATAdir ,Type_list , OUTdir2 ,param2 )

chcekAssoOverlap( OUTdir1,OUTdir2,param, param2 )

Out1 <- OUTdir1
Out2 <- OUTdir2

pGene="KDM6A"
p1 <- GenePlot( pGene ,  forPlot_GeneStats_res[[3]]  )
p2 <- GenePlot( pGene ,  forPlot_GeneStats_res_B2[[3]]  )
p1/p2

pGene="KDM6A"
p2 <-GenePlot( pGene ,  forPlot_GeneStats_res_B2[[3]]  )
ggsave( file = paste0( OUTdir2,"/Plot_",pGene,".pdf" ), p2 , width= 13 , height=4 )

pGene="MAML2"
p2 <-GenePlot( pGene ,  forPlot_GeneStats_res_B2[[3]]  )
ggsave( file = paste0( OUTdir2,"/Plot_",pGene,".pdf" ), p2 , width= 13 , height=4 )

pGene="IRF1"
p2 <-GenePlot( pGene ,  forPlot_GeneStats_res_B2[[3]]  )
ggsave( file = paste0( OUTdir2,"/Plot_",pGene,".pdf" ), p2 , width= 13 , height=4 )


