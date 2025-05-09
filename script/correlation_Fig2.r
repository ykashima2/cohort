library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(purrr)
library(cowplot)
library(patchwork)
set.seed(1234)

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
#rmv_sample <- c("Todai_484",  "Todai_267")
Pearson_dir <-  "output/Figure2/"

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
RNAfile = paste0(RNAdir ,"/RNA_Bulk_mat.txt" )
Exp_mat <- read.table(file = RNAfile, header=T,stringsAsFactors=F,sep ="\t" , row.names=1)
# selcet col of  exprresion values
colnames(Exp_mat) <- gsub( paste0("_",celltype,"_Mean") ,"" , colnames(Exp_mat) )
colnames(Exp_mat) <- gsub( paste0("_Mean") ,"" , colnames(Exp_mat) )

plot_Expl <- Exp_mat %>% 
	mutate( Gene = rownames(Exp_mat)) %>% 
	pivot_longer(-Gene , names_to = "sample" ,values_to="exp" ) %>% 
	left_join( . ,metaInf , by=c("sample"= "Sample") ) %>%
	filter( RNA == 1 ) 

## corrlation 
cor_mat.rna.all <-plot_Expl %>% 
	group_by(Gene) %>% 
	summarise(correlation = cor( exp, age)) %>% 
	ungroup()

cor_mat.rna.old1 <-plot_Expl %>% 
	filter( groupe == "old1" )  %>%
	group_by(Gene) %>% 
	summarise(correlation = cor( exp, age)) %>% 
	ungroup()

cor_mat.rna.old2 <-plot_Expl %>% 
	filter( groupe == "old2" )  %>%
	group_by(Gene) %>% 
	summarise(correlation = cor( exp, age)) %>% 
	ungroup()

cor_mat.rna.young <-plot_Expl %>% 
	filter( groupe == "young" )  %>%
	group_by(Gene) %>% 
	summarise(correlation = cor( exp, age)) %>% 
	ungroup()

cor_mat.rna.old_1_2 <- plot_Expl %>%
        filter( groupe == "old2" |  groupe == "old1" )  %>%
        group_by(Gene) %>%
        summarise(correlation = cor( exp, age)) %>%
        ungroup()

cor_mata.rna_old_F <- plot_Expl %>%
        filter( (groupe == "old2" |  groupe == "old1" )  & sex == "F" )  %>%
        group_by(Gene) %>%
        summarise(correlation = cor( exp, age)) %>%
        ungroup()

cor_mata.rna_old_M <- plot_Expl %>%
        filter( (groupe == "old2" |  groupe == "old1") & sex == "M" )  %>%
        group_by(Gene) %>%
        summarise(correlation = cor( exp, age)) %>%
        ungroup()

cor_mata.rna_F <- plot_Expl %>%
        filter(  sex == "F" )  %>%
        group_by(Gene) %>%
        summarise(correlation = cor( exp, age)) %>%
        ungroup()

cor_mata.rna_M <- plot_Expl %>%
        filter(  sex == "M" )  %>%
        group_by(Gene) %>%
        summarise(correlation = cor( exp, age)) %>%
        ungroup()

dfs <- list(cor_mat.rna.all,cor_mat.rna.old1,cor_mat.rna.old2,cor_mat.rna.young , cor_mat.rna.old_1_2 ,cor_mata.rna_old_M,cor_mata.rna_old_F , cor_mata.rna_M ,cor_mata.rna_F)
# reduce()で繰り返し結合
cor_mat.rna.join <- reduce(dfs, left_join, by = "Gene")
colnames(cor_mat.rna.join) <- c( "Gene" , "all","old1","old2","young" ,"old" ,"Male" ,"Female","Male.all" ,"Female.all")
cor.file.rna = paste0( Pearson_dir,"rna_cor_score.txt")
write.table( cor_mat.rna.join , file=cor.file.rna  , quote=F, col.names=NA , sep="\t")
cor_gene_list <- cor_mat.rna.join %>% dplyr::filter( abs(old1) >0.5 ) %>% select(Gene) %>% unlist


GerRwithP <- function( test_mat, type ,age  ){
	cor_res_list <- apply( test_mat ,2, function(x){
                                   res<-cor.test(x,age);
                                   res_spearman <- cor.test(x,age ,method ="spearman");
                                   res_kendall <- cor.test(x,age ,method ="kendall");
                                   return( c( res$estimate , res$p.value , res_spearman$estimate ,res_spearman$p.value ,res_kendall$estimate ,res_kendall$p.value ))
	} ) %>% t %>% as.data.frame
	colnames( cor_res_list ) <- paste0( rep( paste0(type,".") ,n=6) , c("pearson","pearson.p", "spearman","spearman.p" , "kendall","kendall.p" )  )
	return( cor_res_list)
}



test_mat <-plot_Expl %>% 
	filter(   groupe == "old1" ) %>% 
	pivot_wider( names_from= Gene ,values_from= exp  ) %>% as.data.frame
rownames(test_mat) <- test_mat$sample
age<- test_mat$age
test_mat <- test_mat[, 9:ncol(test_mat) ]
cor_res_list.old1 <- GerRwithP( test_mat, "old1" ,age  )
       
test_mat <-plot_Expl %>%
        filter(   groupe == "old2" ) %>%
        pivot_wider( names_from= Gene ,values_from= exp  ) %>% as.data.frame
rownames(test_mat) <- test_mat$sample
age<- test_mat$age
test_mat <- test_mat[, 9:ncol(test_mat) ]
cor_res_list.old2 <- GerRwithP( test_mat, "old2" ,age  )


test_mat <-plot_Expl %>%
        filter(  groupe == "old1" | groupe == "old2" ) %>%
        pivot_wider( names_from= Gene ,values_from= exp  ) %>% as.data.frame
rownames(test_mat) <- test_mat$sample
age<- test_mat$age
test_mat <- test_mat[, 9:ncol(test_mat) ]
cor_res_list.old <-  GerRwithP( test_mat, "old" ,age  )

cor_res_list.join = cbind(  cor_res_list.old1,cor_res_list.old2,cor_res_list.old )
cor_res_list.join.file = paste0( Pearson_dir,"rna_cor_score_withP.txt")
write.table( cor_res_list.join  %>% as.data.frame %>% mutate(Gene = rownames(cor_res_list.join)) %>% filter( abs(old1.pearson) > 0.5), file= cor_res_list.join.file  , quote=F, col.names=NA , sep="\t")

### M / F
test_mat <-plot_Expl %>%
        filter(  (groupe == "old1" | groupe == "old2" ) & sex =="F" ) %>%
        pivot_wider( names_from= Gene ,values_from= exp  ) %>% as.data.frame
rownames(test_mat) <- test_mat$sample
age<- test_mat$age
test_mat <- test_mat[, 9:ncol(test_mat) ]
cor_res_list.F <-  GerRwithP( test_mat, "Female" ,age  )

test_mat <-plot_Expl %>%
        filter(  (groupe == "old1" | groupe == "old2" ) & sex =="M" ) %>%
        pivot_wider( names_from= Gene ,values_from= exp  ) %>% as.data.frame
rownames(test_mat) <- test_mat$sample
age<- test_mat$age
test_mat <- test_mat[, 9:ncol(test_mat) ]
cor_res_list.M <-  GerRwithP( test_mat, "Male" ,age  )

cor_res_list.joinFM = cbind(  cor_res_list.old1,cor_res_list.F,cor_res_list.M )
cor_res_list.join.file = paste0( Pearson_dir,"rna_cor_score_sex_withP.txt")
cor_res_list.joinFM <-cor_res_list.joinFM %>% 
	as.data.frame %>% 
	mutate(Gene = rownames(cor_res_list.joinFM)) %>% 
	filter( abs(old1.pearson) > 0.5 | abs(Female.pearson) > 0.5 |abs(Male.pearson) > 0.5 )
write.table(  cor_res_list.joinFM, file= cor_res_list.join.file  , quote=F, col.names=NA , sep="\t")


##################################
## ATAC
##################################

# get ron num
ATACfile = paste0(ATACdir , "TSS_", celltype ,"_mat.txt" )
TSS_mat <- read.table(file = ATACfile,  header=T,stringsAsFactors=F,sep ="\t",row.names=1 )
colnames(TSS_mat) <- gsub("_[TBNM].*" ,"" ,colnames(TSS_mat) )

plot_Tss <- TSS_mat %>% 
	    filter( !TSSid %in%  c(TSS_mat$TSSid %>% table %>% as.data.frame %>% filter(Freq > 1) %>% pull(".")) ) %>%
	    pivot_longer( colnames(TSS_mat)[6:ncol(TSS_mat)], names_to = "sample" ,values_to="TSSopen" ) %>% 
	    left_join(.,metaInf , by=c("sample"= "Sample") ) %>%
	    filter( ATAC == 1 )

## corrlation 
cor_mat.atac.all <-plot_Tss %>% 
	group_by(TSSid) %>% 
	summarise(correlation = cor( TSSopen, age)) %>% 
	ungroup()

cor_mat.atac.old1 <-plot_Tss %>% 
	filter( groupe == "old1" )  %>%
	group_by(TSSid) %>% 
	summarise(correlation = cor( TSSopen, age)) %>% 
	ungroup()

cor_mat.atac.old2 <-plot_Tss %>% 
	filter( groupe == "old2" )  %>%
	group_by(TSSid) %>% 
	summarise(correlation = cor( TSSopen, age)) %>% 
	ungroup()

cor_mat.atac.young <-plot_Tss %>% 
	filter( groupe == "young" )  %>%
	group_by(TSSid) %>% 
	summarise(correlation = cor( TSSopen, age)) %>% 
	ungroup()

cor_mat.atac.old_1_2 <-plot_Tss %>% 
	filter( groupe == "old2" | groupe == "old1" )  %>%
	group_by(TSSid) %>% 
	summarise(correlation = cor( TSSopen, age)) %>% 
	ungroup()

cor_mat.atac.F <-plot_Tss %>%
        filter( (groupe == "old2" | groupe == "old1") & sex =="F" )  %>%
        group_by(TSSid) %>%
        summarise(correlation = cor( TSSopen, age)) %>%
        ungroup()

cor_mat.atac.M <-plot_Tss %>%
        filter( (groupe == "old2" | groupe == "old1")  & sex =="M" )  %>%
        group_by(TSSid) %>%
        summarise(correlation = cor( TSSopen, age)) %>%
        ungroup()

cor_mat.atac.Fall <-plot_Tss %>%
        filter(  sex =="F" )  %>%
        group_by(TSSid) %>%
        summarise(correlation = cor( TSSopen, age)) %>%
        ungroup()

cor_mat.atac.Mall <-plot_Tss %>%
        filter(    sex =="M" )  %>%
        group_by(TSSid) %>%
        summarise(correlation = cor( TSSopen, age)) %>%
        ungroup()



dfs <- list(cor_mat.atac.all,cor_mat.atac.old1,cor_mat.atac.old2,cor_mat.atac.young , cor_mat.atac.old_1_2 , cor_mat.atac.M, cor_mat.atac.F,cor_mat.atac.Mall,cor_mat.atac.Fall)
# reduce()で繰り返し結合
cor_mat.atac.join <- reduce(dfs, left_join, by = "TSSid")
colnames(cor_mat.atac.join) <- c( "TSSid" , "all","old1","old2","young" , "old", "Male" ,"Female" ,"Male.all","Female.all" )
cor_mat.atac.join <- left_join ( cor_mat.atac.join , TSS_mat[,1:5] , by = "TSSid")

colnames(cor_mat.atac.join)[1] <- c( "Gene" )
colnames(plot_Tss)[5] <- c( "Gene" )

cor.file.atac = paste0( Pearson_dir,"atac_cor_score.txt")
write.table( cor_mat.atac.join , file=cor.file.atac  , quote=F, col.names=NA , sep="\t")
cor_gene_list.atac <- cor_mat.atac.join %>% dplyr::filter( abs(old1) >0.5 ) %>% select(Gene) %>% unlist

test_mat <-plot_Expl %>%
        filter(   groupe == "old1" ) %>%
        pivot_wider( names_from= Gene ,values_from= exp  ) %>% as.data.frame

test_mat <-plot_Tss %>%
        filter(   groupe == "old1" ) %>%
	select( Gene ,sample ,TSSopen ,age, sex) %>%
        pivot_wider( names_from= Gene ,values_from= TSSopen  ) %>% as.data.frame
rownames(test_mat) <- test_mat$sample
age<- test_mat$age
test_mat <- test_mat[, 4:ncol(test_mat) ]
cor_res_list.atac.old1 <-  GerRwithP( test_mat, "old1" ,age  )

test_mat <-plot_Tss %>%
        filter(   groupe == "old2" ) %>%
	select( Gene ,sample ,TSSopen ,age, sex) %>%
        pivot_wider( names_from= Gene ,values_from= TSSopen  ) %>% as.data.frame
rownames(test_mat) <- test_mat$sample
age<- test_mat$age
test_mat <- test_mat[, 4:ncol(test_mat) ]
cor_res_list.atac.old2 <-  GerRwithP( test_mat, "old2" ,age  )

test_mat <-plot_Tss %>%
        filter( groupe == "old1" | groupe == "old2" ) %>%
	select( Gene ,sample ,TSSopen ,age, sex) %>%
        pivot_wider( names_from= Gene ,values_from= TSSopen  ) %>% as.data.frame
rownames(test_mat) <- test_mat$sample
age<- test_mat$age
test_mat <- test_mat[,4:ncol(test_mat) ]
cor_res_list.atac.old <-  GerRwithP( test_mat, "old" ,age  )



cor_res_list.atac.join = cbind(  cor_res_list.atac.old1,cor_res_list.atac.old2,cor_res_list.atac.old )
cor_res_list.atac.join.file = paste0( Pearson_dir,"atac_cor_score_withP.txt")
write.table( cor_res_list.atac.join  %>% as.data.frame %>% mutate(Gene = rownames(cor_res_list.atac.join)) %>% filter( abs(old1.pearson) > 0.5), file= cor_res_list.atac.join.file  , quote=F, col.names=NA , sep="\t")

### M / F

## Female
test_mat <-plot_Tss %>%
	filter( (groupe == "old2" | groupe == "old1") & sex =="F" )  %>%
        select( Gene ,sample ,TSSopen ,age, sex) %>%
        pivot_wider( names_from= Gene ,values_from= TSSopen  ) %>% as.data.frame
rownames(test_mat) <- test_mat$sample
age<- test_mat$age
test_mat <- test_mat[,4:ncol(test_mat) ]
cor_res_list.atac.F <-  GerRwithP( test_mat, "Female" ,age  )

#Make
test_mat <-plot_Tss %>%
	filter( (groupe == "old2" | groupe == "old1") & sex =="M" )  %>%
        select( Gene ,sample ,TSSopen ,age, sex) %>%
        pivot_wider( names_from= Gene ,values_from= TSSopen  ) %>% as.data.frame
rownames(test_mat) <- test_mat$sample
age<- test_mat$age
test_mat <- test_mat[,4:ncol(test_mat) ]
cor_res_list.atac.M <-  GerRwithP( test_mat, "Male" ,age  )

cor_res_list.atac.joinFM = cbind(  cor_res_list.atac.old1,cor_res_list.atac.F,cor_res_list.atac.M )
cor_res_list.atac.join.file = paste0( Pearson_dir,"atac_cor_score_sex_withP.txt")
cor_res_list.atac.joinFM <-cor_res_list.atac.joinFM %>% 
	as.data.frame %>% 
	mutate(Gene = rownames(cor_res_list.atac.joinFM)) %>% 
	filter( abs(old1.pearson) > 0.5 | abs(Female.pearson) > 0.5 |abs(Male.pearson) > 0.5 )
write.table(  cor_res_list.atac.joinFM, file= cor_res_list.atac.join.file  , quote=F, col.names=NA , sep="\t")



###################################
# correlation Table
###################################

overlap_table <- c()
for( th in c( 0.2 , 0.3 , 0.4 ,0.5 )) {
	tmp_overlap_table <- cbind( 
		       old1_old2 = c(RNA = cor_mat.rna.join  %>% dplyr::filter( Gene %in% cor_gene_list)       %>% dplyr::filter( abs(old2) > th ) %>% nrow()   ,
				     ATAC = cor_mat.atac.join %>% dplyr::filter( Gene %in% cor_gene_list.atac ) %>% dplyr::filter( abs(old2) > th ) %>% nrow() ) ,
			old1_old  = c(RNA = cor_mat.rna.join  %>% dplyr::filter( Gene %in% cor_gene_list)       %>% dplyr::filter( abs(old)  > th  ) %>% nrow()  , 
				     ATAC = cor_mat.atac.join %>% dplyr::filter( Gene %in% cor_gene_list.atac ) %>% dplyr::filter( abs(old)  > th  ) %>% nrow() ),
			merges    = c(RNA = cor_mat.rna.join  %>% dplyr::filter( Gene %in% cor_gene_list)       %>% dplyr::filter( abs(all)  > th  ) %>% nrow()  ,
				     ATAC = cor_mat.atac.join %>% dplyr::filter( Gene %in% cor_gene_list.atac ) %>% dplyr::filter( abs(all)  > th  ) %>% nrow() ),
			old1_F    = c(RNA = cor_mat.rna.join  %>% dplyr::filter( Gene %in% cor_gene_list)       %>% dplyr::filter( abs(Female)  > th  ) %>% nrow()  ,
				     ATAC = cor_mat.atac.join %>% dplyr::filter( Gene %in% cor_gene_list.atac ) %>% dplyr::filter( abs(Female)  > th  ) %>% nrow() ),
			old1_M    = c(RNA = cor_mat.rna.join  %>% dplyr::filter( Gene %in% cor_gene_list)       %>% dplyr::filter( abs(Male)  > th  ) %>% nrow()  ,
				     ATAC = cor_mat.atac.join %>% dplyr::filter( Gene %in% cor_gene_list.atac ) %>% dplyr::filter( abs(Male)  > th  ) %>% nrow() ) ,
			parm      = c(RNA= th , ATAC=th)
			)
	overlap_table <-rbind( overlap_table ,tmp_overlap_table)

}

write.table(  overlap_table , file=  paste0( Pearson_dir,"/overlap_matrix.txt")   , quote=F, col.names=NA , sep="\t")


################ correlation ここまで

plot_gene <- c( "CD248","RALGPS2", "KLHL29"  , "CDKN2A", "CDKN1A" )
for ( gene in plot_gene ) {
#	gene = "KLHL29"
	r_val = cor_mat.rna.join %>%  filter(Gene == gene) %>% pull(all) %>% round(.,2)
	g_p1.r <-plot_Expl %>% 
		filter(Gene == gene )  %>% 
	#	filter( groupe == "old2" |  groupe == "old1" )  %>% 
		ggplot(aes(x=age  , y= exp , col=groupe )) + 
			geom_smooth(method = lm, se = T,  col = "black", aes(group = NULL))+
			#geom_smooth(method = lm, se = T)+
			labs(title = paste0( "Gene:", gene ,", R= " , r_val )) +
			geom_point() +
			theme_classic()+
			theme(text = element_text(size = 20) , axis.title.x = element_blank())

	r_val = cor_mat.atac.join %>%  filter(Gene == gene) %>% pull(all) %>% round(.,2)

# atac
	g_p1.a =plot_Tss%>%
	        filter(Gene == gene )  %>%
	        ggplot(aes(x=age  , y= TSSopen , col=groupe )) +
		        geom_smooth(method = lm, se = T,  col = "black", aes(group = NULL))+
		        #geom_smooth(method = lm, se = T)+
			labs(title = paste0( "Gene:", gene ,", R= " , r_val )) +
	                geom_point() +
                	theme_classic()+
        	        theme(text = element_text(size = 20) , axis.title.x = element_blank())
	p <- g_p1.r + g_p1.a
	fname <- paste0( Pearson_dir,"Gene_scatPlot_", gene,".pdf")
	ggsave(p , file =fname  , height= 4,width= 10  )

}




