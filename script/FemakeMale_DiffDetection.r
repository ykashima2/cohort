#setwd("work/utokyo_250327/")
#install.packages(c("dplyr","tibble","tidyr","ggplot2","purrr","cowplot","patchwork"))

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(purrr)
library(cowplot)
library(patchwork)
library(ggpubr)
library(rlang)
library(ggpmisc)

set.seed(1234)
options(width = 500)


#---------------------------------------------------------------------------
# Setup input file
#---------------------------------------------------------------------------
INPUTDIR <- "data/"
SampleListfile <- paste0( INPUTDIR , "Sample_inf.txt")
AnnotationFILE <- paste0( INPUTDIR , "Cell_annotaion.txt")
# INPUT exp data
DATAdir <- "Input/"
RNAdir <- paste0( DATAdir , "bulk/")
ATACdir <- paste0( DATAdir , "bulk/")

# OUTPUT
Pearson_dir <-  "output/Figure2/"

# OUTPUT
Pearson_dir <- paste0( DATAdir , "Figure2_sex_diff")

for(dir_path in c( Pearson_dir)){
  if (!dir.exists(dir_path)) { # ディレクトリが存在するか確認
    dir.create(dir_path) # ディレクトリが存在しない場合、作成
    print(paste("ディレクトリ", dir_path, "を作成しました"))
  }
}

# Read Meta
metaInf <- read.table(file = SampleListfile, header=T,stringsAsFactors=F,sep ="\t")
metaInf <- metaInf  %>% filter( type != "2nd")
metaInf <- metaInf %>%
                mutate(groupe=if_else(type == "1st" ,"old1",
                                        if_else(age >=55 , "old2","young"  )  ) )

Type_list <- list( "Mono" ,"Tcell","Bcell" ,"NKcell" , "Bulk" )

#---------------------------------------------------------------------------
# Read input file
#---------------------------------------------------------------------------
age_data.rna  <- metaInf  %>% filter( RNA == 1 ) %>% dplyr::select(age) %>% unlist
age_data.atac <- metaInf  %>% filter( ATAC == 1 ) %>% dplyr::select(age) %>% unlist
names(age_data.rna)  <- metaInf  %>% filter( RNA == 1 ) %>% dplyr::select(Sample) %>% unlist
names(age_data.atac) <- metaInf  %>% filter( ATAC == 1 ) %>% dplyr::select(Sample) %>% unlist

##################
## RNA Plot
################
celltype <- "Bulk"
RNAfile = paste0(RNAdir ,"RNA_", celltype ,"_mat.txt" )
Exp_mat <- read.table(file = RNAfile, header=T,stringsAsFactors=F,sep ="\t" , row.names=1)
# selcet col of  exprresion values
colnames(Exp_mat) <- gsub( paste0("_",celltype,"_Mean") ,"" , colnames(Exp_mat) )
colnames(Exp_mat) <- gsub( paste0("_Mean") ,"" , colnames(Exp_mat) )

Box_col = c( Male_Young="#7A81FF" , Male_Old="#011993" , Female_Young = "#FF85FF" , Female_Old= "#FF2F92"  )
Dot_col = c( M = "#4F53F7" , F = "#FCB4CF" )

plot_Expl <- Exp_mat %>%
  mutate( Gene = rownames(Exp_mat)) %>%
  pivot_longer(-Gene , names_to = "sample" ,values_to="exp" ) %>%
  left_join( . ,metaInf , by=c("sample"= "Sample") ) %>%
  mutate( cate  = if_else( groupe== "old1" | groupe =="old2" , "old" , "young"))  %>% 
  mutate( cate2 = if_else( sex == "M" ,
			if_else( groupe== "old1" | groupe =="old2" , "Male_Old" , "Male_Young") , 
			if_else( groupe== "old1" | groupe =="old2" , "Female_Old" , "Female_Young") 
							 )) %>%
  mutate( cate2 = factor( cate2, levels = names(Box_col)) ,
          sex = factor( sex ,levels = c("M","F"))
        ) %>%
  filter( RNA == 1 ) %>% 
  select( Gene, sample , exp , age,sex , groupe , cate ,cate2 ) 

all_Gene_stat <-  plot_Expl %>%
   group_by(Gene) %>%
   summarise(
   pWil_yo  = wilcox.test(exp[cate == "young"], exp[cate == "old"])$p.value,
   pWil_Myo = wilcox.test(exp[cate2 == "Male_Young"], exp[cate2 == "Male_Old"])$p.value,
   pWil_Fyo = wilcox.test(exp[cate2 == "Female_Young"], exp[cate2 == "Female_Old"])$p.value,
   pWil_FM = wilcox.test(exp[sex == "F"], exp[sex == "M"])$p.value,
   cor_y  = cor(exp[cate == "young"] , age[cate == "young"] ),
   cor_o  = cor(exp[cate == "old"]   , age[cate == "old"] ),
   cor_My = cor(exp[cate2 == "Male_Young" ]  , age[cate2 == "Male_Young" ] ),
   cor_Mo = cor(exp[cate2 == "Male_Old"  ]   , age[cate2 == "Male_Old"  ] ),
   cor_Fy = cor(exp[cate2 == "Female_Young"] , age[cate2 == "Female_Young"] ),
   cor_Fo = cor(exp[cate2 == "Female_Old" ]  , age[cate2 == "Female_Old" ] ),
   corM = cor(exp[ sex=="M"] , age[sex=="M"] ),
   corF = cor(exp[ sex=="F"] , age[sex=="F"] ),
   mean_y  = mean(exp[cate == "young"]),
   mean_o  = mean(exp[cate == "old"]),
   mean_My = mean(exp[cate2 == "Male_Young" ]),
   mean_Mo = mean(exp[cate2 == "Male_Old"  ]),
   mean_Fy = mean(exp[cate2 == "Female_Young" ]),
   mean_Fo = mean(exp[cate2 == "Female_Old"  ]),
   mean_M  = mean(exp[sex == "M"]),
   mean_F  = mean(exp[sex == "F"])
  )   %>% 
   mutate(
    pWil_yo_adj  = p.adjust(pWil_yo, method = "BH"), # BH法による多重検定補正
    pWil_Myo_adj = p.adjust(pWil_Myo, method = "BH"), # BH法による多重検定補正
    pWil_Fyo_adj = p.adjust(pWil_Fyo, method = "BH"), # BH法による多重検定補正
    pWil_FM_adj  = p.adjust(pWil_FM, method = "BH"), # BH法による多重検定補正
    fc_yo   = log( mean_y /mean_o,2),
    fc_M_yo = log( mean_My /mean_Mo,2),    
    fc_F_yo = log( mean_Fy /mean_Fo,2),
    fc_FM   = log( mean_F /mean_M,2),
  ) %>% ungroup()

all_Gene_stat.use <- all_Gene_stat %>% 
	select(Gene , pWil_yo_adj,pWil_Myo_adj,pWil_Fyo_adj,pWil_FM_adj,corM, corF,
	        mean_y,mean_o,mean_My,mean_Mo,mean_Fy,mean_Fo,mean_M,mean_F,fc_yo,fc_M_yo,fc_F_yo,fc_FM) %>% 
    mutate( diffCor = abs( corM - corF ))

write.table( all_Gene_stat.use , file= paste0( Pearson_dir , "/allGene_Mean_Cor_wilcox_stats.slim.txt")  , quote=F, col.names=NA , sep="\t")


testPlotS <-function(gene,plot_Expl ,stats , Pearson_dir){
  stats <- stats %>% filter(Gene == gene)%>% select( Gene , pWil_yo_adj,pWil_Myo_adj,pWil_Fyo_adj,pWil_FM_adj,corM, corF)

  max_val <- plot_Expl %>%
	filter(Gene == gene) %>%
	select(exp) %>%
 	unlist%>%
    max()

  p1<- plot_Expl %>% 
  filter(Gene == gene) %>%
  ggplot() +
  geom_boxplot(aes(x = cate2, y = exp, fill = cate2)) +
  scale_fill_manual(values = Box_col) +
  theme_classic() +
  ggtitle( paste0("Boxplot : ", gene)) + 
  theme(text = element_text(size = 20), axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotate("text", x = 1.5, y = max_val * 1.05, label = paste("p =", stats$pWil_Myo_adj %>% round(.,2)), size = 5) +
  annotate("text", x = 3.5, y = max_val * 1.05, label = paste("p =", stats$pWil_Fyo_adj %>% round(.,2) ), size = 5)

  p2<-plot_Expl %>% filter( Gene == gene) %>%
    ggplot( aes(x = age ,y=exp , col=sex)  )+geom_point()+
    scale_colour_manual(values = Dot_col ) +
    geom_smooth(method="lm", aes(color=sex) ,se=F ) +
	theme_classic() +
    ggtitle( paste0("Correlation : ", gene , "\nMale:R = " , stats$corM %>% round(.,2) ,  "\nFemale:R = " , stats$corF %>% round(.,2) )) +
    theme(text = element_text(size = 20) , axis.title.x = element_blank()) +
      stat_poly_eq(formula = y ~ x, aes(label = paste(stat(eq.label) )),  parse = TRUE)
  print(p1+p2)
  outfig <- paste0( Pearson_dir ,"/plot_" ,gene, "PlotS.pdf" )
  ggsave( p1+p2 , file = outfig , height= 7,width= 12 )
  return(p1+p2)
}

testPlotS( "FRMD3" , plot_Expl, all_Gene_stat.use ,Pearson_dir)
testPlotS( "IL7R" , plot_Expl, all_Gene_stat.use , Pearson_dir)
testPlotS( "EFNA5" , plot_Expl, all_Gene_stat.use ,Pearson_dir)
testPlotS( "SPTB" , plot_Expl, all_Gene_stat.use ,Pearson_dir)
testPlotS( "FAM13A" , plot_Expl, all_Gene_stat.use ,Pearson_dir)
testPlotS( "CDKN2A" , plot_Expl, all_Gene_stat.use ,Pearson_dir)
testPlotS( "CDKN1A" , plot_Expl, all_Gene_stat.use ,Pearson_dir)

#女性で変化
filterparm =  " abs(corF) >= 0.3 & ( corM > -0.2 &  corM < 0.2)  & (  mean_F>0.005 )"
selectCol = c( "Gene","pWil_Myo_adj","pWil_Fyo_adj","corM","corF","mean_F","mean_M","fc_F_yo","fc_M_yo","fc_FM"  )
all_Gene_stat.use %>% filter(!!parse_expr(filterparm)) %>% select( all_of(selectCol))  %>% as.data.frame %>%
 write.table( . , file= paste0( Pearson_dir , "/Female_diffGene.txt")  , quote=F, col.names=NA , sep="\t")

#男性で変化
filterparm =  " abs(corM) >= 0.3 & ( corF > -0.2 &  corF < 0.2)  & (  mean_M>0.005 )"
selectCol = c( "Gene","pWil_Myo_adj","pWil_Fyo_adj","corM","corF","mean_F","mean_M","fc_F_yo","fc_M_yo","fc_FM"  )
all_Gene_stat.use %>% filter(!!parse_expr(filterparm)) %>% select( all_of(selectCol))  %>% as.data.frame %>%
 write.table( . , file= paste0( Pearson_dir , "/Male_diffGene.txt")  , quote=F, col.names=NA , sep="\t")

