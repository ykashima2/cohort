library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(purrr)
library(cowplot)
library(patchwork)
set.seed(1234)
options(scipen=100) 

###########################################
#Function
###########################################

# Create Randam Sample Combination 
# example
# age_ranges <- list(c(60, 69), c(70, 79), c(80, 89))
# sample_size <- 5
# num_sets <- 100#
# randam_sample_set <- random_sample(metaInf, age_ranges, sample_size, num_sets,seed = 123)

random_sample <- function(df, age_ranges, sample_size, num_sets, seed = 123) {
	set.seed(seed)

	# 年代ごとのレコード数を確認
	age_counts <- sapply(age_ranges, function(range) {
				     sum(df$age >= range[1] & df$age <= range[2])})

	# 抽出する人数がレコード数を超えている場合の警告
	if (any(age_counts < sample_size)) {
		warning("抽出する人数がレコード数を超えています。")
		return(NULL)
	}

	# 組み合わせの上限を確認
	max_combinations <- prod(sapply(age_ranges, function(range) {
						choose(sum(df$age >= range[1] & df$age <= range[2]), sample_size)}))
	if (max_combinations < num_sets) {
		warning("作成するセット数が組み合わせの上限を超えています。")
		return(NULL)
	}

	# 無作為抽出セットを作成
	sets <- list()
	for (i in 1:num_sets) {
		set <- list()
		for (range in age_ranges) {
			age_subset <- df %>% filter(age >= range[1] & age <= range[2])
			set[[paste0(range[1], "-", range[2])]] <- age_subset[sample(nrow(age_subset), sample_size), ]
		}
		sets[[i]] <- set
	}

	# ユニークセットを確認
	unique_sets <- unique(sets)
	if (length(unique_sets) < num_sets) {
		warning("ユニークなセットが作成できませんでした。")
	}
	return(unique_sets)
}


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
Pearson_dir <- paste0( DATAdir , "Figure2_CV/")

for(dir_path in c( Pearson_dir  )){
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

metaInfF <- metaInf %>% dplyr::filter( sex == "F"  )
metaInfM <- metaInf %>% dplyr::filter( sex == "M"  )

Type_list <- list( "Mono" ,"Tcell","Bcell" ,"NKcell" , "Bulk" )

# cor_gene_list <- old1_corlate %>% filter( abs(age ) > 0.5) %>% select(X) %>% unlist

#---------------------------------------------------------------------------
# Read input file
#---------------------------------------------------------------------------
age_data.rna  <- metaInf  %>% filter( RNA == 1 ) %>% dplyr::select(age) %>% unlist
age_data.atac <- metaInf  %>% filter( ATAC == 1 ) %>% dplyr::select(age) %>% unlist
names(age_data.rna)  <- metaInf  %>% filter( RNA == 1 ) %>% dplyr::select(Sample) %>% unlist
names(age_data.atac) <- metaInf  %>% filter( ATAC == 1 ) %>% dplyr::select(Sample) %>% unlist


#---------------------------------------------------------------------------
## AGE Randam Sampling
#---------------------------------------------------------------------------
age_ranges <- list(c(60, 89))
sample_size <- 15
num_sets <- 50

randam_sample_set <- random_sample(metaInf, list( c(60 ,69) , c(70,79) ,c(80 ,89) ) , 10 , num_sets,seed = 123)
randam_sample_set.F <- random_sample(metaInfF, age_ranges, sample_size, num_sets,seed = 123)
randam_sample_set.M <- random_sample(metaInfM, age_ranges, sample_size, num_sets,seed = 123)


##################
## RNA correlation
################
celltype <- "Bulk"
RNAfile = paste0(RNAdir ,"RNA_", celltype ,"_mat.txt" )
Exp_mat <- read.table(file = RNAfile,  header=T,stringsAsFactors=F,sep ="\t" , row.names=1 )

colnames(Exp_mat) <- gsub( paste0("_",celltype,"_Mean") ,"" , colnames(Exp_mat) )
colnames(Exp_mat) <- gsub( paste0("_Mean") ,"" , colnames(Exp_mat) )

plot_Expl <- Exp_mat %>% 
	mutate( Gene = rownames(Exp_mat)) %>% 
	pivot_longer( -Gene , names_to = "sample" , values_to="exp" ) %>% 
	left_join( . ,metaInf , by=c("sample"= "Sample") ) %>%
	filter( RNA == 1 ) 

## corrlation 
cor_mat.rna.old1 <-plot_Expl %>%
        filter( groupe == "old1" )  %>%
        group_by(Gene) %>%
        summarise(correlation = cor( exp, age)) %>%
        ungroup()

cor_gene_list <- cor_mat.rna.old1 %>% dplyr::filter( abs(correlation) >0.5 ) %>% select(Gene) %>% unlist
plot_Expl <- plot_Expl  %>% filter( Gene %in% cor_gene_list )

make_cv_p <-  function( randam_sample_set , plot_Expl ,out ){
	cor_mat.rna.va <- lapply( randam_sample_set , function(x, plot_Expl){
		sample_list <- do.call(rbind , x) %>% select(Sample) %>% unlist
		cor_mat.rna <-plot_Expl %>% 
			filter( sample %in%  sample_list) %>% 
			group_by(Gene) %>%
			summarise(correlation = cor( exp, age )) %>%
			ungroup()%>% 
			return()
		}, plot_Expl=plot_Expl ) %>% do.call(cbind ,. )

	cor_mat.rna.va <- cor_mat.rna.va[ , c(1 , seq(2,ncol(cor_mat.rna.va),2)) ]

	val_pval <- apply( cor_mat.rna.va ,1,function(x){ 
		xval <- x[2:length(x)] %>% as.numeric()
		mean_r = mean(xval)
		if (mean_r > 0){
		        res_t0.5 <- t.test( xval , mu = 0.5, alternative = "greater")
		        res_t0.4 <- t.test( xval , mu = 0.4, alternative = "greater")
		        res_t0.3 <- t.test( xval , mu = 0.3, alternative = "greater")
		        res_t0.2 <- t.test( xval , mu = 0.2, alternative = "less")
		} else {
	        	res_t0.5 <- t.test( xval , mu = -0.5, alternative = "less")
		        res_t0.4 <- t.test( xval , mu = -0.4, alternative = "less")
		        res_t0.3 <- t.test( xval , mu = -0.3, alternative = "less")
		        res_t0.2 <- t.test( xval , mu = -0.2, alternative = "less")
		}
		return(c (Gene = x[1], mean_r = mean_r ,
			  pval0.5= res_t0.5$p.value , num_over0.5 = sum(abs(xval)>0.5) , 
			  pval0.4= res_t0.4$p.value , num_over0.4 = sum(abs(xval)>0.4),
			  pval0.3= res_t0.3$p.value , num_over0.3 = sum(abs(xval)>0.3),
			  pval0.2= res_t0.2$p.value , num_over0.2 = sum(abs(xval)>0.2)
			  )
		) 
	} ) %>% t %>% as.data.frame
	write.table( val_pval , file= out  , quote=F, col.names=NA , sep="\t")
	return(val_pval )
}

cros.val.rna.file = paste0( Pearson_dir,"rna_cor_cross_validataon_100_Fe.txt")
val_pval.F <- make_cv_p( randam_sample_set.F , plot_Expl , cros.val.rna.file )
cros.val.rna.file = paste0( Pearson_dir,"rna_cor_cross_validataon_100_Me.txt")
val_pval.M <- make_cv_p( randam_sample_set.M , plot_Expl , cros.val.rna.file )
cros.val.rna.file = paste0( Pearson_dir,"rna_cor_cross_validataon_100_old.txt")
val_pval <- make_cv_p( randam_sample_set , plot_Expl , cros.val.rna.file )

##################################
## ATAC
##################################
# get ron num
ATACfile = paste0(ATACdir , "TSS_", celltype ,"_mat.txt" )
TSS_mat <- read.table(file = ATACfile,  header=T,stringsAsFactors=F,sep ="\t",row.names=1 )
colnames(TSS_mat) <- gsub("_[TBNM].*" ,"" ,colnames(TSS_mat) )

# get ron num
plot_Tss <- TSS_mat %>% 
	    filter( !TSSid %in%  c(TSS_mat$TSSid %>% table %>% as.data.frame %>% filter(Freq > 1) %>% pull(".")) ) %>%
	    pivot_longer( colnames(TSS_mat)[6:ncol(TSS_mat)], names_to = "sample" ,values_to="TSSopen" ) %>% 
	    left_join(.,metaInf , by=c("sample"= "Sample") ) %>%
	    filter( ATAC == 1 )

## corrlation 
cor_mat.atac.old1 <-plot_Tss %>% 
	filter( groupe == "old1" )  %>%
	group_by(TSSid) %>% 
	summarise(correlation = cor( TSSopen, age)) %>% 
	ungroup()
cor_gene_list.atac <- cor_mat.atac.old1 %>% dplyr::filter( abs(correlation) >0.5 ) %>% select(TSSid) %>% unlist

plot_Tss <- plot_Tss %>% filter( TSSid %in% cor_gene_list.atac)



make_cv_p.atac <-  function( randam_sample_set , plot_Tss ,out ){

	cor_mat.atac.va<- lapply( randam_sample_set, function(x, plot_Tss){
               sample_list <- do.call(rbind , x) %>% select(Sample) %>% unlist
               cor_mat.rna <-plot_Tss %>%
                       filter( sample %in%  sample_list) %>%
                       group_by(TSSid) %>%
                       summarise(correlation = cor( TSSopen, age )) %>%
                       ungroup()%>%
                       return()
        }, plot_Tss=plot_Tss ) %>% do.call(cbind ,. )

	cor_mat.atac.va <- cor_mat.atac.va[ , c(1 , seq(2,ncol(cor_mat.atac.va),2)) ]

	val_pval.atac <- apply( cor_mat.atac.va,1,function(x){
              xval <- x[2:length(x)] %>% as.numeric()
	      xval[is.na(xval)] <- 0
	      mean_r = mean(xval)
              if (mean_r > 0){
                      res_t0.5 <- t.test( xval , mu = 0.5, alternative = "greater")
                      res_t0.4 <- t.test( xval , mu = 0.4, alternative = "greater")
                      res_t0.3 <- t.test( xval , mu = 0.3, alternative = "greater")
                      res_t0.2 <- t.test( xval , mu = 0.2, alternative = "greater")
              } else {
                      res_t0.5 <- t.test( xval , mu = -0.5, alternative = "less")
                      res_t0.4 <- t.test( xval , mu = -0.4, alternative = "less")
                      res_t0.3 <- t.test( xval , mu = -0.3, alternative = "less")
                      res_t0.2 <- t.test( xval , mu = -0.2, alternative = "less")
              }
	      return(c (Gene = x[1], mean_r = mean_r ,
			pval0.5= res_t0.5$p.value , num_over0.5 = sum(abs(xval)>0.5) , 
			pval0.4= res_t0.4$p.value , num_over0.4 = sum(abs(xval)>0.4) ,
			pval0.3= res_t0.3$p.value , num_over0.3 = sum(abs(xval)>0.3) ,
			pval0.2= res_t0.2$p.value , num_over0.2 = sum(abs(xval)>0.2) 
			)
	      )
	} ) %>% t %>% as.data.frame

	write.table( val_pval.atac , file= out  , quote=F, col.names=NA , sep="\t")
	return(val_pval.atac)
}

cros.val.atac.file = paste0( Pearson_dir,"atac_cor_cross_validataon_100_Fe.txt")
val_pval.atac.F <- make_cv_p.atac( randam_sample_set.F , plot_Tss , cros.val.atac.file )
cros.val.atac.file = paste0( Pearson_dir,"atac_cor_cross_validataon_100_Me.txt")
val_pval.atac.M <- make_cv_p.atac( randam_sample_set.M , plot_Tss , cros.val.atac.file )
cros.val.atac.file = paste0( Pearson_dir,"atac_cor_cross_validataon_100_old.txt")
val_pval.atac <- make_cv_p.atac( randam_sample_set , plot_Tss , cros.val.atac.file )

## Table
pva_th <- 0.5
rna_size <- length(cor_gene_list)
atac_size <-length(cor_gene_list.atac)
th_list <- rev(c("0.5"=3,"0.4"= 5,"0.3"=7,"0.2" = 9))
i =3
table_data <- c()

for( i in c(1:length(th_list))  ) {

	th = names(th_list)[i]
	pos = th_list[i]
	tmp_table_data <- cbind(
		cate ="RNA" ,
		th =  th ,
		Female = sum(val_pval.F[,pos] < pva_th) ,
		Female_r = sum(val_pval.F[,pos] < pva_th)/rna_size,
		Male   = sum(val_pval.M[,pos] < pva_th),
		Male_r   = sum(val_pval.M[,pos] < pva_th) /rna_size,
		"All(old1+old2)" = sum(val_pval[,pos] < pva_th) ,
		"All(old1+old2)_r" = sum(val_pval[,pos] < pva_th) /rna_size
	)
	table_data = rbind( table_data , tmp_table_data)

	tmp_table_data <- cbind(
		cate = "ATAC" ,
		th =  th ,
		Female = sum(val_pval.atac.F[,pos] < pva_th) ,
		Female_r = sum(val_pval.atac.F[,pos] < pva_th)/atac_size,
		Male   = sum(val_pval.atac.M[,pos] < pva_th),
		Male_r   = sum(val_pval.atac.M[,pos] < pva_th)/atac_size,
		"All(old1+old2)" =sum(val_pval.atac[,pos] < pva_th) ,
		"All(old1+old2)_r" =sum(val_pval.atac[,pos] < pva_th)/atac_size
	)
	table_data = rbind( table_data , tmp_table_data)
}
write.table(  table_data , file=  paste0( Pearson_dir,"/cross_validation_pval.txt")   , quote=F, col.names=NA , sep="\t")

