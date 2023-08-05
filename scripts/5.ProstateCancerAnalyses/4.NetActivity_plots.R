#'#################################################################################
#'#################################################################################
#'                     NetActivity on PROMOTE vis2 Figures
#'                        Figure 5 and Supp. Figure 14   
#'                       
#' 1. Load data from differential gene set scores analysis (ProstateCancerAnalyses/3.NetActivity_PROMOTEvis2)
#' 2. Load data from differential gene analysis (ProstateCancerAnalyses/2.DEA_PROMOTEvis2R)
#' 3. Load input data (Prepare/Data9.create_SE)
#' 4. Plot Function
#' 5. Create NetActivity plots
#' 
#'#################################################################################
#'#################################################################################

library(SummarizedExperiment)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(reshape2)
library(forcats)
library(stringr)
library(ggtext)
library(broom)
library(org.Hs.eg.db)



# Load data ---------------------------------------------------------------
message("Loading data files")
load("results/mCRPC_analyses/NetActivity_GTEX_PROMOTEv2.Rdata") # scores, input_SE

topTab <- read.csv('results/mCRPC_analyses/SupTable6_promoteV2_diff_paths_TTCajdustlog_GTEx.csv',row.names = 1)

ddsSE <- readRDS('data/promoteV2_SE.rds')
ddsSE <- subset(ddsSE, select =  Biopsy_site_Visit_2  == "bone")
rownames(ddsSE) <- rowData(ddsSE)$new_symbol

vst_out <- varianceStabilizingTransformation(ddsSE)

DEgenes_DESeq <- read.csv("results/mCRPC_analyses/SupTable7_promoteV2bone_DEgenes_TTC&2LFC_DESeq.csv",header = T)[,1]
# Also update symbol genes of DE genes (does not match the order)
DEgenes_DESeq2 <- rowData(ddsSE)$new_symbol[which(names(rowData(ddsSE)$new_symbol)%in%DEgenes_DESeq)]


# PLOT FUNCTION ####

plot_weights_scatter <- function(nr_path=1,
                                 topTab=NULL,
                                 scores=NULL,
                                 DEgenes_DESeq=NULL,
                                 vst_out=NULL,
                                 nr_genes=3,
                                 LM=T){
  # Input checks
  if (is.null(topTab)){stop('Need to set a topTable result')}
  if (is.null(scores)){stop('Need to set a scores result')}
  if (is.null(DEgenes_DESeq)){stop('Need to set a DEgenes_DESeq vector')}
  if (is.null(vst_out)){stop('Need to set a vst_out object')}
  path <- rownames(topTab)[nr_path]
  path_name <- str_to_title(topTab$GeneSetName[nr_path])
  weights <- rowData(scores)[path, ]$Weights_SYMBOL[[1]]
  patt <- character(length = length(weights))
  
  
  df <- data.frame(weight=weights,
                   gene=names(weights)) %>% 
    mutate(Direction = factor(ifelse(weight > 0, "Positive", "Negative"))) %>%
    mutate(weight=abs(weight)) %>% 
    arrange(desc(weight)) %>% 
    mutate(gene=factor(gene,levels = gene)) %>% 
    mutate(Direction=factor(Direction,levels = c("Negative","Positive")))
  
  # BOXPLOT gene expression
  keep_genes <- as.character(df$gene)[1:nr_genes]
  keep_genes <- keep_genes[keep_genes%in%rownames(vst_out)]
  df_exp <- assay(vst_out[keep_genes,])
  df_exp_r <- melt(df_exp,value.name = "Expression")
  df_exp_r$Var1 <- factor(df_exp_r$Var1,levels=keep_genes)
  condition <- data.frame(sample=rownames(colData(scores)),
                          TTC.adjust=scores$TTC.adjust,
                          TTC.adjust_log=scores$TTC.adjust_log,
                          condition=scores$condition)
  
  df_boxplot1 <- df_exp_r %>% left_join(condition,by=c("Var2"="sample"))
  colnames(df_boxplot1) <- c("gene","sample","Expression","TTC.adjust","TTC.adjust_log","condition")
  
  ## FOR ALL GENES GRID and LM
  if (nr_genes>3){
    
    # ARRANGEMENTS FOR PLOTTING:
    # factor to show low-middle-high from left to right
    df_boxplot1$condition <- relevel( df_boxplot1$condition, ref="Middle")
    df_boxplot1$condition <- relevel( df_boxplot1$condition, ref="Low")
    
    # Get lm fit for all genes separatelly
    fit <- list()
    for (gene in keep_genes){
      fit[[gene]] <- summary(lm(data=df_boxplot1[df_boxplot1$gene==gene,] , formula= TTC.adjust_log~Expression))
      
    }
    # Get adj.r.sq
    rsq <- do.call(rbind,as.data.frame(do.call(rbind,fit))$adj.r.squared)
    # Get intercept pvalue and significant genes
    coef <- as.data.frame(do.call(rbind,fit))$coefficients
    significant_models <- sapply(coef, function(x) {x[2,4]<0.05 })
    significant_models <- names(significant_models)[significant_models]
    # Create new gene labels with significance for plot (facet_wrap)
    df_boxplot1$gene_sig <-  as.character(df_boxplot1$gene)
    df_boxplot1$gene_sig <- ifelse(df_boxplot1$gene_sig%in%significant_models,
                                   paste(keep_genes,"** Adj.R squared",signif(rsq,2)),
                                   paste(df_boxplot1$gene_sig,"Adj.R squared",signif(rsq,2)))
    keep_genes_sig <- ifelse(keep_genes%in%significant_models,
                             paste(keep_genes,"** Adj.R squared",signif(rsq,2)),
                             paste(keep_genes,"Adj.R squared",signif(rsq,2)))
    
    df_boxplot1$gene_sig <- factor(df_boxplot1$gene_sig, levels=keep_genes_sig)
    
    # ALL GENES MODEL FIT
    fit_all <- lm(data=df_boxplot1 , formula= TTC.adjust_log~Expression)
    # f <- summary(fit_all)$fstatistic
    # model_pval <- pf(f[1],f[2],f[3],lower.tail=F)
    
    caption=paste0("Adj R2 = ", signif(summary(fit_all)$adj.r.squared, 3),";",
                   " Intercept = ", signif(fit_all$coef[[1]],3 ),";\n",
                   " Slope = ", signif(fit_all$coef[[2]], 3), ";",
                   " Pvalue = ", signif(summary(fit_all)$coef[2,4], 3),";")
    # BOXPLOT 
    p1 <- df_boxplot1 %>%
      ggplot(aes(x = condition, y = Expression, fill = condition)) +
      geom_boxplot(show.legend = F) +
      theme_bw() +
      theme(axis.title.x =  element_text(size=7.5,hjust=-0.05,vjust=6),
            strip.text = element_text(size = 6))+
      # scale_fill_brewer(palette = "Oranges")+
      scale_fill_manual(values=c("#C18541", "#838581", "#4484C1"))+
      labs(x="",y="Normalized Gene Expression",
           title="Distribution of Gene Expression by response group",
           subtitle =paste(path_name,"-",path),
           caption=caption)+
      facet_wrap(vars(gene_sig),scales = "free",ncol = 3)
    
    # SCATTERPLOT
    p2 <- df_boxplot1 %>%
      ggplot(aes(x = TTC.adjust, y = Expression, colour = condition) )+
      geom_point(size=1, show.legend = T)+
      theme_bw() +
      theme(legend.position = c(0.82, 0.05), legend.key.size = unit(0.7,"cm"),
            strip.text = element_text(size = 10),
            axis.title.x =  element_text(size=10,hjust=0.1,vjust=0),
            legend.background = element_rect(fill="white",
                                             size=0.5))+
      scale_colour_manual(values=c("#C18541", "#838581", "#4484C1"))+
      labs(x="Time To Treatment Change [days]", y="Normalized Gene Expression",
           title= "Pathway Genes' expression",
           subtitle = paste(path_name,"-",path),
           colour="TTTC group")+
      stat_smooth(data= df_boxplot1, method = "lm", col = "turquoise3", se = F)+
      facet_wrap(vars(gene_sig),scales = "free_y", ncol= 3)
    
    # Terminate  
    return(list(boxlpota_all=p1, scatterplot_all=p2))
  }
  
  
  # FOR 3item-GRID
  # Top 3 genes scatterplot
  p <- df_boxplot1 %>%
    ggplot(aes(x = TTC.adjust, y = Expression, colour = condition) )+
    geom_point(show.legend = T)+
    theme_bw() +
    scale_x_continuous(trans='log10') +
    scale_colour_manual(values=c("#4484C1","#838581" ,"#C18541"))+
    labs(x="Time To Treatment Change [days]",y="Normalized Gene Expression",
         title="Top genes' expression",
         colour="TTTC group")+
    stat_smooth(data= df_boxplot1, method = "lm", col = "turquoise3", se = F)+
    # stat_smooth(data=subset(df_boxplot1,gene==keep_genes[1]),method = "lm", col = "turquoise3",se = F)+
    facet_wrap(vars(gene), scales = "free")
  
  # LM to get Adjusted R2 for first gene in pathway
  if (var(df_boxplot1[df_boxplot1$gene==keep_genes[1],'Expression'])==0){
    print("var 0")
    adj_r_df <- data.frame(topgene=NA)
    adj_r_df$gene_pval<- NA
    adj_r_df$signf_gene <- NA
    
  }else{
    fit <- lm(data=df_boxplot1[df_boxplot1$gene==keep_genes[1],] , formula= TTC.adjust_log~Expression)
    adj_r_df <- data.frame(topgene=as.numeric(signif(summary(fit)$adj.r.squared, 3)))
    adj_r_df$gene_pval<- as.numeric(summary(fit)$coef[2,4])
    adj_r_df$signf_gene <- as.numeric(summary(fit)$coef[2,4]<0.05)
  }
  # ann_text <- data.frame(gene = factor(keep_genes[1],levels = keep_genes),
  #                        sample=condition$sample[1],
  #                        Expression =  max(df_boxplot1[df_boxplot1$gene==keep_genes[1],]$Expression),
  #                        TTC.adjust_log =  max(df_boxplot1$TTC.adjust_log)-1,
  #                        condition=condition$condition[1])
  # if(summary(fit)$coef[2,4]<0.05){
  #   label= paste("Adj R2 =",signif(summary(fit)$adj.r.squared, 3),";\n",
  #                     "pvalue =",signif(summary(fit)$coef[2,4], 3),"*")
  #   }else{label= paste("Adj R2 =",signif(summary(fit)$adj.r.squared, 3),";\n Non significant")}
  # plabel <- p + geom_text(data = ann_text,label=label, color="black")
  
  
  # Adj R2 for genes plot
  # if(LM){
  #   adj_r <- c()
  #   for (gene in keep_genes){
  #     fit <- lm(data=df_boxplot1[df_boxplot1$gene==gene,] , formula= TTC.adjust_log~Expression)
  #     adj_r <- c(adj_r,as.numeric(signif(summary(fit)$adj.r.squared, 3)))
  #   }
  #   names(adj_r) <- keep_genes
  # }
  
  # PATHWAY SCORES Scatter plot
  df_boxplot2 <- data.frame(Expression=assay(scores[path, ]),
                            condition = scores$condition,
                            TTC.adjust = scores$TTC.adjust,
                            TTC.adjust_log = scores$TTC.adjust_log)
  
  
  p2 <- df_boxplot2 %>%
    ggplot(aes(x = TTC.adjust, y = Expression, colour = condition)) +
    geom_point(show.legend = T) +
    theme_bw() +
    scale_x_continuous(trans='log10') +
    scale_colour_manual(values = c("#4484C1", "#838581","#C18541"))+
    labs(x="Time To Treatment Change [days]",y="Activity scores",
         title=paste("Linear model on", path),
         colour="TTTC group")
  
  p3 <- p2
  # If want LM plot in scores
  if(LM){
    fit_score <- lm(data=df_boxplot2 , formula= TTC.adjust_log~Expression)
    
    p3 <- p2+stat_smooth(method = "lm", col = "turquoise3",se = F)+
      theme(plot.caption.position = "plot",
            plot.caption = element_text(vjust=10))+
      xlab("Time To Treatment Change [days]")+
      scale_x_continuous(trans='log10') +
      ggtitle(paste("Linear model on", path)) +
      ylab(paste("Activity scores"))+
      guides(color=guide_legend(title="TTTC group",
                                override.aes=list(3)))
    adj_r_df$score <- as.numeric(signif(summary(fit_score)$adj.r.squared, 3))
    adj_r_df$score_pval <- as.numeric(summary(fit_score)$coef[2,4])
    adj_r_df$signf_path <- as.numeric(summary(fit_score)$coef[2,4]<0.05)
    adj_r_df$path <- path
    adj_r_df$pathname <- path_name
    adj_r_df$shape <- paste(as.character(adj_r_df$signf_gene),as.character(adj_r_df$signf_path), collapse = "/")
    
  }
  
  # DATA Arrangements for BARPLOT
  if (sum(df$gene%in%DEgenes_DESeq)){
    idx <- which(df$gene%in%DEgenes_DESeq)
    patt[idx] <- "DE"
    df <- df %>% mutate(DEgenes=patt)
    print('hey')
  }else{
    df <- df %>% mutate(DEgenes=patt)
  }
  bold <- ifelse(df$DEgenes == "DE", "bold", "plain")
  size <- ifelse(df$DEgenes == "DE", 10, 7)
  
  if(length(unique(df$Direction))==2){
    values_color <- c("#fe4a49","#2ab7ca")
  }else{
    values_color <- ifelse(df$Direction[1]=="Positive","#fe4a49","#2ab7ca")
  }
  
  # BARPLOT scores
  q <-ggplot(df,aes(x = gene, y = weight, fill= Direction)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    scale_fill_manual(values=values_color)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = c(.85, .8),
          legend.box.background = element_rect(color="black", size=1),)+
    labs(x="Genes",y="Weight",title="Gene weights in Gene Set") +
    guides(x =  guide_axis(angle = 45))+
    # Add label to DE genes
    geom_point( shape=25, size=4,
                data=df %>% filter(DEgenes=="DE"), # Filter data first
                aes(x=gene, y = weight+0.03), show.legend = F)+
    # Add black margin to DE genes
    geom_bar(data=df %>% filter(DEgenes=="DE"),colour= "black",
             aes(x = gene, y = weight),size=1,
             stat = "identity",show.legend = F)+
    # Make labels of DE genes bold
    theme(axis.text.x = element_markdown(face = bold,size=size))
  
  # MAIN TITLE
  title <- ggdraw() +
    draw_label(
      paste0("PATHWAY: ", path_name, " (",path,")"),
      fontface = 'bold',
      x = 0,
      hjust = -0.05
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 0)
    )
  
  # END
  return(list(genes=p, scores=p3, weights=q, title=title, adj_r_df = adj_r_df))
  
  
}

# Plotting ####


#### GO:0051988
# 3D GRID layout
nr_path <- 1
path <- gsub(":","_" ,rownames(topTab)[nr_path])
plots_out <- plot_weights_scatter(nr_path = nr_path,
                                  topTab = topTab,
                                  scores = scores,
                                  DEgenes_DESeq = DEgenes_DESeq2,
                                  vst_out=vst_out,
                                  nr_genes=3,
                                  LM=T)
plots <- plot_grid(plots_out$scores,plots_out$genes,nrow=2,align = "h")
plots2 <- plot_grid(plots,plots_out$weights,rel_widths = c(1,1),align = "v", ncol=2)


pdf(paste0('figures/promoteV2_3Dgrid__',path,'.pdf'),width = 13)
print(plot_grid(plots_out$title,plots2,ncol=1,# rel_heights values control vertical title margins
                rel_heights = c(0.1, 1)))
dev.off()

# **Modified layout/legend in Illustrator (legend) then saved into .png from there

#### All genes scater plot 

nr_path <- 1
path <- gsub(":","_" ,rownames(topTab)[nr_path])
nr_genes <-   length(rowData(scores[rownames(topTab)[nr_path]])$Weights_SYMBOL[[1]])
plots_out <- plot_weights_scatter(nr_path = nr_path,
                                  topTab = topTab,
                                  scores = scores,
                                  DEgenes_DESeq = DEgenes_DESeq2,
                                  vst_out=vst_out,
                                  nr_genes=nr_genes,
                                  LM=T)

plots_out$scatterplot_all

pdf(paste0('figures/all_genes_pathway_',path,'_scatter.pdf'),width = 10, height = 10)
plots_out$scatterplot_all
dev.off()

# **Modified layout/legend in Illustrator (moved last plot to center) then saved into .png from there


