#library(dplyr)
#library(ggplot2)
#library(tidyverse)
#library(ggridges)
#library(tibble)
#library(patchwork)
#library(ggpmisc)
#library(ggrepel)
set.seed(1234)


# この後に控えているのは相関解析
# だから、ある細胞種の RNA とATACのファイルを返して、抽出する必要があるなら、抽出したmatrixを返す、このmatrixは次の関数の引数になる。

Data.util <-  function(DATAdir , sample_set, CellType){
	
	file.rna  <- paste0( DATAdir,"/RNA_", CellType , "_mat.txt" )
        file.atac <- paste0( DATAdir, "/TSS_", CellType , "_mat.txt" )
        RNA_mat <- read.table(file = file.rna  , header=T,stringsAsFactors=F,sep ="\t" ,row.names=1)
        ATAC_mat <- read.table(file = file.atac , header=T,stringsAsFactors=F,sep ="\t" , row.names=1)

	cname_forcheck.rna <- gsub( paste0("_",CellType,".*")  , "" , colnames(RNA_mat)[2:ncol(RNA_mat)])
	RNA_mat <-RNA_mat[ ,c( rep(TRUE,1) ,cname_forcheck.rna %in% sample_set)]

	cname_forcheck.atac <- gsub( paste0("_",CellType,".*")  , "" , colnames(ATAC_mat)[8:ncol(ATAC_mat)])
	ATAC_mat <-ATAC_mat[,c( rep(TRUE,7) ,cname_forcheck.atac %in% sample_set)]
	
	ATAC_mat <- ATAC_mat %>% 
			  pivot_longer(cols = starts_with("Todai"), names_to = "Sample", values_to = "Score") %>%
			  group_by(Gene, TSSid) %>% 
			  summarise(mean_score = mean(Score), .groups = 'drop') %>% 
			  group_by(Gene) %>%
			  slice(which.max(mean_score)) %>% 
			  ungroup() %>% 
			  select(Gene, TSSid) %>% 
			  left_join(ATAC_mat, by = c("Gene", "TSSid")) %>% 
			  as.data.frame

	## Gene preparations
	rownames(ATAC_mat) <-ATAC_mat$Gene
	rownames(RNA_mat)  <- RNA_mat$Gene
	MultiomeGene = intersect( rownames(ATAC_mat) , rownames(RNA_mat) )

	ATAC_mat <- ATAC_mat[ MultiomeGene, ] 
	RNA_mat  <- RNA_mat[ MultiomeGene, ] 

	return( list(RNA_mat , ATAC_mat ) )
}

MakeGeneStats <- function( mat_data, metaInf , CellType , param   ) {

	RNA_filt_cellrate  = param$RNA_filt_cellrate
	RNA_Filter1  = param$RNA_Filter1
	ATAC_Open_th = param$ATAC_Open_th

	# Gene Selection
	RNA_mat.use <- mat_data[[1]]
	ATAC_mat.use <- mat_data[[2]]
	ATAC_mat.use <- ATAC_mat.use[,c(1 , 8:ncol(ATAC_mat.use)) ]

	RNAstats  <- t(apply(RNA_mat.use[ , 2:ncol(RNA_mat.use) ] , 1, function(x ) {
			return( c( RNAave = median(x[seq(1 , length(x),2 )])  ,
					   RNAnumCount = sum( x[seq(2 , length(x),2 )] >= RNA_filt_cellrate  ) # 10% 以上のcellが発現している検体の数
			)) 
	}))

	ATACstats <-  t(apply(ATAC_mat.use[,2:ncol(ATAC_mat.use)] ,1, function(x){
				return( c( ATACave= median(x), ATACnumCount = sum(x >= ATAC_Open_th)  ) )
			}  ))

	# RNAstats と ATACstats で発現している状態を確認 
	GeneStats <- cbind( ATACstats,RNAstats) %>% as.data.frame()
	print("dim...")
	print(dim(ATAC_mat.use))
	print(dim(RNA_mat.use))
	###
	tmpATAC <- ATAC_mat.use %>% 
		dplyr::filter( Gene %in% rownames(GeneStats %>% filter( ATACnumCount >= 4 & RNAnumCount >=4 )) ) %>%
		pivot_longer( cols = -Gene ,values_to = "ATAC" ) %>% mutate(Sample = gsub( paste0("_",CellType),"",name) )
	tmpRNA  <- RNA_mat.use  %>% 
		 dplyr::select( Gene , matches("_Mean$")) %>% 
		 dplyr::filter( Gene %in% rownames(GeneStats %>% filter( ATACnumCount >= 4 & RNAnumCount >=4 )) ) %>%
		pivot_longer( cols = -Gene ,values_to = "RNA"  ) %>% mutate(Sample = gsub( paste0("_",CellType,"_Mean")  ,"",name) )
	#print( cbind( tmpATAC , RNA=tmpRNA$RNA ) %>% select(Sample) %>%  unique() %>% unlist)

	print("dim...")
	print(dim(tmpATAC))
	print(dim(tmpRNA))

	print("test")
	tmpGEX_ATAC <- full_join(tmpATAC , tmpRNA[,c("Gene","Sample","RNA" )]  ,by= c("Gene","Sample")) %>% 
		rowwise() %>% mutate( Age =  metaInf[ metaInf$Sample == Sample ,"age" ] ) 
	print("mod")
	#tmpGEX_ATAC <- cbind( tmpATAC , RNA=tmpRNA$RNA ) %>% rowwise() %>% mutate( Age =  metaInf[ metaInf$Sample == Sample ,"age" ] )
	print("mod")
	return( list( GeneStats , tmpGEX_ATAC ))
}


LoesssAss <- function( mat_data , GeneStats_res , metaInf , CellType){
	RNA_mat.use <- mat_data[[1]]
	ATAC_mat <- mat_data[[2]]
	ATAC_mat.use <- ATAC_mat[,c(1 , 8:ncol(ATAC_mat)) ]

	GeneStats <- GeneStats_res[[1]]

	## Loesso
	test.ATAC_mat.use <- ATAC_mat.use[ rownames(GeneStats %>% filter( ATACnumCount >= 4 & RNAnumCount >=4 ))  ,c(-1) ] 
	metaInf.tmp <- metaInf
	rownames(metaInf.tmp) <- metaInf.tmp$Sample
	Age <- metaInf.tmp[gsub(  paste0("_",CellType,"$")  ,"", colnames(test.ATAC_mat.use)) , "age" ]

	# ATC
	ATACfitted_range <- t(apply( test.ATAC_mat.use  , 1 ,
		function(x){
			dat <- cbind( ATAC = x , Age =Age  ) %>% as.data.frame()
			loessRes <- loess(formula =  ATAC ~ Age, dat  , model = TRUE,span = 0.75 ,method = "loess" , degree = 1)
			return(c( ATAC=max(loessRes$fitted) - min(loessRes$fitted) , ATACmedian = median(x) , ATACrse=loessRes$s ,ATACmin = min(loessRes$fitted) ))
		}
	) )
	ATACfitted_range <- cbind( TSSid = ATAC_mat[ ATAC_mat[,1] %in%  rownames(ATACfitted_range) ,2 ] , ATACfitted_range ) %>% as.data.frame
 
	## RNA
	test.RNA_mat.use <- RNA_mat.use[ rownames(GeneStats %>% filter( ATACnumCount >= 4 & RNAnumCount >=4 ))  ,-1]
	num_sample <- ncol(test.RNA_mat.use)
	Age <- metaInf.tmp[gsub( paste0("_",CellType,".*") ,"", colnames(test.RNA_mat.use)[seq(1,num_sample,2)] ) , "age" ]
	RNAfitted_range <- t(apply( test.RNA_mat.use  , 1 ,
		function(x){
			dat <- cbind( RNA = x[ seq(1,num_sample,2)] , Age =Age  ) %>% as.data.frame()
			loessRes <- loess(formula =  RNA ~ Age, dat  , model = TRUE,span = 0.75 ,method = "loess" , degree = 1)
			return( c(range= max(loessRes$fitted) - min(loessRes$fitted) , median = median(loessRes$fitted) ,
				fc = max(loessRes$fitted)/min(loessRes$fitted), max = max(loessRes$fitted) ,  rario_median = median( x[ seq(2,num_sample,2) ]  )  ))
		}
	))
	AgeAsso.summa <- cbind(  ATACfitted_range ,  RNAfitted_range ) %>% as.data.frame
	return( AgeAsso.summa )
}

GenePlot <- function(pGene , tmpGEX_ATAC){
	p1<-tmpGEX_ATAC %>% 
		filter( Gene == pGene ) %>% 
		ggplot(aes( x= Age,y=ATAC ,label =Sample )) + 
			geom_point( colour = "red" ) +
			geom_smooth(method = "loess",method.args= list(degree =1), span = 0.75, colour ="red" ) + 
			theme_classic() +
			theme(text = element_text(size = 20) , axis.title.x = element_blank())


	p2<-tmpGEX_ATAC %>% 
		filter( Gene == pGene ) %>% 
		ggplot(aes( x= Age,y=RNA  ,label =Sample )) +
		       	geom_point( colour = "blue" )   +
			geom_smooth(method = "loess",method.args= list(degree = 1), span = 0.75 , colour = "blue" ) + 
			theme_classic() +
			theme(text = element_text(size = 20) , axis.title.x = element_blank())
	
	p3<-tmpGEX_ATAC %>% 
		filter( Gene == pGene ) %>% 
		ggplot(aes( x= ATAC,y=RNA,label =Sample )) + 
			geom_point( ) + 
			theme_classic() +
			stat_smooth(method = "lm" , colour = "#990099")+ 
			stat_poly_eq(formula = x ~ y , aes(label = paste(stat(eq.label))),parse = TRUE)

	print(cor( (tmpGEX_ATAC %>% filter( Gene == pGene ))$ATAC , (tmpGEX_ATAC %>% filter( Gene == pGene ))$RNA, method ="spearman"))
	print(cor( (tmpGEX_ATAC %>% filter( Gene == pGene ))$ATAC , (tmpGEX_ATAC %>% filter( Gene == pGene ))$RNA, method ="pearson"))
	return( p1+p2 + p3)
}


#doall <-  function( mat_data , metaInf  ,Type_list , OUTdir  , RNA_filt_cellrate , RNA_Filter1 , ATAC_Open_th , ATAC_Filter1 ,Ass_ATAC_delt,Ass_RNAmed,Ass_RNAfc ){
doall <-  function(  metaInf , sample_set ,DATAdir ,Type_list , OUTdir ,param ){
	RNA_filt_cellrate = param$RNA_filt_cellrate
	RNA_Filter1 = param$RNA_Filter1
	ATAC_Open_th = param$ATAC_Open_th
	ATAC_Filter1 = param$ATAC_Filter1
	Ass_ATAC_delt = param$Ass_ATAC_delt
	Ass_RNAmed = param$Ass_RNAmed
	Ass_RNAfc = param$Ass_RNAfc

	Stats1 <-c()
	Stats2 <-c()
	Asso_Gene <- c()
	All_Gene_cand<-c()
	forPlot_GeneStats_res <- vector("list", length(Type_list))

	#for ( CellType in Type_list ) {
	for ( i in c(1:length(Type_list))) {
		CellType = Type_list[i]
		print( paste0("doing ...",  CellType))
		mat_data <-  Data.util(DATAdir, sample_set, CellType)
		print( "finish reading input file")
		GeneStats_res <- MakeGeneStats ( mat_data, metaInf = metaInf , CellType , param )
		print( "finish GeneStats")
		AgeAsso.summa <-  LoesssAss( mat_data , GeneStats_res , metaInf, CellType )
		print( "finish Loess")
		GeneStats     <-  GeneStats_res[[1]]
		tmpGEX_ATAC   <-  GeneStats_res[[2]] 
		forPlot_GeneStats_res[[i]] <- tmpGEX_ATAC

		gene_stats_sum <-c(
		    ATAC_Open   =  GeneStats %>% filter( ATACnumCount >= ATAC_Filter1 ) %>% nrow ,
		    RNA_Express =  GeneStats %>% filter( RNAnumCount  >=RNA_Filter1 )  %>% nrow ,
			ATAC_Open_RNA_Expression = GeneStats %>% filter( ATACnumCount >= ATAC_Filter1 & RNAnumCount >= RNA_Filter1 ) %>% nrow ,
			ATAC_Open_RNA_not_Expression = GeneStats %>% filter( ATACnumCount >= ATAC_Filter1 & RNAnumCount <  RNA_Filter1 ) %>% nrow ,
			ATAC_Close_RNA_not_express    = GeneStats %>% filter( ATACnumCount <  ATAC_Filter1 & RNAnumCount <  RNA_Filter1 ) %>% nrow ,
			ATAC_Close_RNA_Expression     = GeneStats %>% filter( ATACnumCount <  ATAC_Filter1 & RNAnumCount >= RNA_Filter1 ) %>% nrow
		)

		gene_stats_sum2 <- c(
			ATAC = AgeAsso.summa %>% filter(  ATAC > Ass_ATAC_delt  )  %>% nrow,
			RNA  = AgeAsso.summa %>% filter(  median > Ass_RNAmed   & fc > Ass_RNAfc  ) %>% nrow,
			Both = AgeAsso.summa %>% filter(  median > Ass_RNAmed & ATAC > Ass_ATAC_delt  & fc>Ass_RNAfc )  %>% nrow 
		)
		Asso_Gene <- rbind( Asso_Gene , 
				    AgeAsso.summa %>% 
					   filter(  median > Ass_RNAmed & ATAC > Ass_ATAC_delt  & fc>Ass_RNAfc ) %>% 
					   mutate( Type=CellType ) %>% 
					   rownames_to_column(var = "Gene") 
				   )
		All_Gene_cand <- rbind( All_Gene_cand , AgeAsso.summa %>% mutate( Type=CellType) %>% rownames_to_column(var = "Gene") )

		Stats1 <- rbind( Stats1 ,`CellType` = gene_stats_sum)
		Stats2 <- rbind( Stats2 ,`CellType` = gene_stats_sum2)
	}

	rownames(Stats1) <-Type_list
	rownames(Stats2) <-Type_list

	write.table( cbind( Stats1, Stats2) , file=paste0(OUTdir,"AssotiationSummary.txt") , quote=F, col.names=NA , sep="\t")
	write.table( All_Gene_cand , file= paste0(OUTdir,"AgeAsso_AllGeneMAt.txt") , quote=F, col.names=NA , sep="\t")
	write.table( Asso_Gene , file= paste0(OUTdir,"AgeAsso_FilteredGeneMAt.txt") , quote=F, col.names=NA , sep="\t")

	write.table( All_Gene_cand %>% filter(  ATAC > Ass_ATAC_delt  ) , file=paste0(OUTdir,"AgeAsso_FilteredGeneMAt_ATAC.txt") , quote=F, col.names=NA , sep="\t" )
	write.table( All_Gene_cand %>% filter(median > Ass_RNAmed   & fc > Ass_RNAfc ) , file=paste0(OUTdir,"AgeAsso_FilteredGeneMAt_RNA.txt") , quote=F, col.names=NA , sep="\t")
	write.table( All_Gene_cand %>% filter(median > Ass_RNAmed & ATAC > Ass_ATAC_delt  & fc>Ass_RNAfc ) , file=paste0(OUTdir,"AgeAsso_FilteredGeneMAt_Both.txt" ), quote=F, col.names=NA , sep="\t")
	return(forPlot_GeneStats_res)
}


chcekAssoOverlap <- function( Out1 , Out2 , param , param2){
	RNA_filt_cellrate = param$RNA_filt_cellrate
        RNA_Filter1 = param$RNA_Filter1
        ATAC_Open_th = param$ATAC_Open_th
        ATAC_Filter1 = param$ATAC_Filter1
        Ass_ATAC_delt = param$Ass_ATAC_delt
        Ass_RNAmed = param$Ass_RNAmed
        Ass_RNAfc = param$Ass_RNAfc

	Out1.file = paste0(Out1,"AgeAsso_AllGeneMAt.txt")
	Out2.file = paste0(Out2,"AgeAsso_AllGeneMAt.txt")

	Out1.data <- read.table(file = Out1.file  , header=T,stringsAsFactors=F,sep ="\t" ,row.names=1)
        Out2.data <- read.table(file = Out2.file, header=T,stringsAsFactors=F,sep ="\t" , row.names=1)

	gene_stats_sum <- c()
	for ( CellType in Type_list ) {
		tmp_gene_stats_sum <- c( 
			ATAC = intersect( 
				 Out1.data %>% filter( Type == CellType &  ATAC > param$Ass_ATAC_delt ) %>% pull(Gene) , 
				 Out2.data %>% filter( Type == CellType &  ATAC > param2$Ass_ATAC_delt ) %>% pull(Gene)
				) %>% length ,
			RNA = intersect(
				Out1.data %>% filter( Type == CellType & median > param$Ass_RNAmed  & fc > param$Ass_RNAfc ) %>% pull(Gene) ,
				Out2.data %>% filter( Type == CellType & median > param2$Ass_RNAmed & fc > param2$Ass_RNAfc ) %>% pull(Gene) 
				) %>% length ,
			Both = intersect(
				 Out1.data %>% filter( Type == CellType & median > param$Ass_RNAmed  & ATAC > param$Ass_ATAC_delt  & fc > param$Ass_RNAfc ) %>% pull(Gene) ,
				 Out2.data %>% filter( Type == CellType & median > param2$Ass_RNAmed & ATAC > param2$Ass_ATAC_delt & fc > param2$Ass_RNAfc) %>% pull(Gene) 
				) %>% length
		)
		gene_stats_sum <-rbind(gene_stats_sum , tmp_gene_stats_sum)
	}
	rownames(gene_stats_sum) <- Type_list
	return(gene_stats_sum)

}





