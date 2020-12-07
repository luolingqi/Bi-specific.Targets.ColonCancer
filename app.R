
options(repos = BiocManager::repositories())

library(shiny)
# Data preprocessing
library(GEOquery)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(reshape2)

################################
# Data extraction from GEO and GTEx
#################################
types = c("primary", "metastasis.liver", "metastasis.lung")

gsms_metastasis.lung <- "XXXXXXXXXXXXXXXXXXXXX1XX1X0X0XXXXXXXX0XX0XX1XXXXX0XX0XX0X0X0XX0X01X0X0X1X0XXX0X0X0XXXXXXXXX0XX1XXXX1XX0XX1XXXXXX1XXXXXXXXXXXXX1XXXXXXXXXXXXXXXXXXXX1XXXX11X1XX1X1XXXXX1XXX111XXXXXXXX0XXXXX00XXXXXXXXXXXXXXXX0XXX0XXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX0XXX0XX0XX0XXXXXXX0X0XXX0XXXXXXXX0XXXX0XXXX0X0XXXXX0XXX0XXXX0X0X0XXXXXXXX0XX0XXXXX0XXXXXXXX0XXX0XX0X0XXXXX0XXX0XXXXX0X0XXX00XXXXXXXXXXXXXX"

gsms_metastasis.liver <- "X1111111X1111111XX1XXXXXXX0X0XXXXXXXX0XX0XXXXXXXX0XX0XX0X0X0XX0X0XX0X0XXX0XXX0X0X0XXXXXXXXX0XXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX1XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX11XXX1X01XXXX00XXX1X1XX11XXX1XX011X0XXXXXXXXXX0XX1XXXX1XXXXXXXXXXX1X11111XXXXX0XXX0XX0XX0XXXXXXX0X0X1X0XX1X11XX0XX1X0XXXX0X0XXXXX0XXX0XXXX0X0X0XXXXXXXX0XX0XXX1X0X1XXXXXX0XXX0XX0X0XX1XX0X1X0X1XXX0X0X1X00X1XXXXXXXXXXXX"

gsms_primary <- "XXXXXXXX1XXXXXXXXXX1XXXXX101011X1X1X1011011X1X11101X01X0101011010X1010XX101X101010X1X11X111011X1111XX101XX1111X1X111X1111111XXX1X11111111X1XX1X1X11X1111XX1X1XX1X1111XX111XXX1XXX11XX0X111X00111X1XX1XXX11XXX0XX10111X1X11XX011X1XXXX1X111X1XXXXXXXXXXX11111011101101101111X11010XX10XXX1XX11011X101111010111110X110111101010111XX1X10110XXXX10XX11111101X10X10101XX110XX10XX1X1010XX100XXXXXXXXXXXXXX"

if (! file.exists("data.RData")) {
  # load series and platform data from GEO
  gset <- getGEO("GSE41258", GSEMatrix =TRUE, AnnotGPL=FALSE)
  if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  
  # make proper column names to match toptable 
  fvarLabels(gset) <- make.names(fvarLabels(gset))
  
  
  for (type in types){
    sml <- c()
    for (i in 1:nchar(get(paste0("gsms_",type)))) { sml[i] <- substr(get(paste0("gsms_",type)),i,i) }
    
    # eliminate samples marked as "X"
    sel <- which(sml != "X")
    sml <- sml[sel]
    #assign(paste0("gset_",type), gset[ ,sel])
    temp <- gset[ ,sel]
    
    # log2 transform
    ex <- exprs(temp)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
      (qx[6]-qx[1] > 50 && qx[2] > 0) ||
      (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    
    if (LogC) { ex[which(ex <= 0)] <- NaN
    exprs(temp) <- log2(ex) }
    
    sml <- paste("G", sml, sep="")    # set group names
    fl <- as.factor(sml)
    
    tmp_grp <- gsub("G1",paste(type,"Tumor", sep = "."),sml)
    assign(paste("groupLabel",type,sep = "_"), gsub("G0","Normal",tmp_grp))
    
    temp$description <- fl
    #temp$group <- groupLabel
    
    assign(paste0("gset_",type), temp)
  }
  # extract the TPM median values from GTEx normals
  gtex_exp <- read.table(file = file.path(getwd(),"GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"),  header = T, sep = "\t")
  # Save all the processed data objects in RData to save time next time
  save(types, gtex_exp, 
       gset_metastasis.liver, gset_metastasis.lung, gset_primary, gset, 
       groupLabel_primary, groupLabel_metastasis.liver, groupLabel_metastasis.lung, 
       file = "data.RData")
}else {
  load("data.RData")
}


df_pairs <- read.table(file = file.path(getwd(),"gene_pairs_for_exp.txt"), header = F, sep = "\t", quote = "", stringsAsFactors = F)

########################################################
# Define the funtion of draw histogram and density plot
########################################################
drawExpHist <- function (x,n) { 
  # declare an empty dataframe
  merged.df <- data.frame()
  merged.groupLabel <- c()
  # for a gene pair, subset and combine expression DF for normal, primary, met.liver, met.lung
  for (type in types) {
    mat <- exprs(get(paste0("gset_",type)))[fData(get(paste0("gset_",type)))$Gene.Symbol %in% x,]
    # median out all the probe sets in one genes and pad with zeros for not expressed genes
    probe_IDs <- row.names(mat)
    mat <- as.data.frame(mat)
    pheno <- fData(get(paste0("gset_",type)))
    mat$row.names <- pheno[pheno$Gene.Symbol %in% x,]$Gene.Symbol
    df_genes_list <- as.data.frame(mat %>% group_by(row.names) %>% 
                                     summarise_all(list(~median(.)))
    ) %>% column_to_rownames(var = "row.names")
    
    groupLabel <- get(paste("groupLabel",type,sep = "_"))
    # Convert to TPM space from intensity space
    for (idx in 1:length(x)) { 
      # array expression is in log2 scale, so TPM from GTEX needs to be transformed to log2
      # minus used instead of divide, due to the log scale
      factor <- apply(df_genes_list[x[idx],groupLabel=="Normal"],1,median) -
        log2(apply(gtex_exp[gtex_exp$Description==x[idx],
                            grepl("Colon",names(gtex_exp), ignore.case = T)], 1, median))
      # the factor converting intensity to TPM for each gene
      df_genes_list[x[idx],] <- df_genes_list[x[idx],] - factor
    }
    
    
    # Merge both dataframe and sample grouplabel
    if (dim(merged.df)[1] == 0){
      merged.df <- df_genes_list
      merged.groupLabel <- groupLabel
    }else {
      merged.groupLabel <- c(merged.groupLabel,groupLabel[!names(df_genes_list) %in% names(merged.df)])
      
      merged.df <- cbind(merged.df, df_genes_list[,!names(df_genes_list) %in% names(merged.df)])
      
    }
    # convert value from log2 scale to normal scale in TPM
    #merged.df <- 2^merged.df
  }
  
  # melt the dataframe
  df.plot <- data.frame(t(merged.df),group=merged.groupLabel) %>%
    mutate(group=factor(group, levels=c("Normal","primary.Tumor",
                                        "metastasis.liver.Tumor",
                                        "metastasis.lung.Tumor"))) %>%
    mutate(sampleID=names(merged.df)) %>%
    melt(variable.name= "gene",
         value.name = "expression",
         id.vars = c("sampleID","group"))
  
  # Draw histogram by cancer types
  theme_set(
    theme_classic() + 
      theme(legend.position = "top") +
      theme(legend.text = element_text(size = 7))
  )
  
  gene1 <- x[1]
  gene2 <- x[2]
  
  p_hist_g1 <- ggplot(df.plot[df.plot$gene == gene1,], aes(x = expression)) + 
    geom_histogram(aes(color=group, fill=group),
                   position = "identity", bins = n, alpha = 0.4) +
    xlab("Expression, log2(TPM)") + 
    ggtitle(paste0("Histogram of expression for the gene ",gene1)) + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12)) +
    scale_color_manual(values = c("#0000FF", "#FF0000","#A0A0A0","#FF9933")) +
    scale_fill_manual(values = c("#0000FF", "#FF0000","#A0A0A0","#FF9933")) 
  
  p_density_g1 <- ggplot(df.plot[df.plot$gene == gene1,], aes(x = expression)) + 
    geom_density(aes(color=group, fill=group), alpha=.2) +  
    xlab("Expression, log2(TPM)") + 
    ggtitle(paste0("Kernel Density Estimation of expression for the gene ",gene1)) + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12)) +
    scale_color_manual(values = c("#0000FF", "#FF0000","#A0A0A0","#FF9933")) +
    scale_fill_manual(values = c("#0000FF", "#FF0000","#A0A0A0","#FF9933")) 
  
  p_hist_g2 <- ggplot(df.plot[df.plot$gene == gene2,], aes(x = expression)) + 
    geom_histogram(aes(color=group, fill=group),
                   position = "identity", bins = n, alpha = 0.4) +
    xlab("Expression, log2(TPM)") + 
    ggtitle(paste0("Histogram of expression for the gene ",gene2)) + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12)) +
    scale_color_manual(values = c("#0000FF", "#FF0000","#A0A0A0","#FF9933")) +
    scale_fill_manual(values = c("#0000FF", "#FF0000","#A0A0A0","#FF9933")) 
  
  p_density_g2 <- ggplot(df.plot[df.plot$gene == gene2,], aes(x = expression)) + 
    geom_density(aes(color=group, fill=group), alpha=.2) +  
    xlab("Expression, log2(TPM)") + 
    ggtitle(paste0("Kernel Density Estimation of expression for the gene ",gene2)) + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12)) +
    scale_color_manual(values = c("#0000FF", "#FF0000","#A0A0A0","#FF9933")) +
    scale_fill_manual(values = c("#0000FF", "#FF0000","#A0A0A0","#FF9933")) 
  
  ggarrange(p_hist_g1,p_hist_g2,p_density_g1,p_density_g2, ncol = 2, nrow = 2)
  
}

###########################################################
# Shiny App Development
###########################################################


ui <- fluidPage(
  
   # App title -----
   titlePanel("Expression Histogram for a SELECTED gene pair among sample types"),

   
		# Input: Slider for the gene pair selection ----
		inputPanel(
			selectInput("genePair", label = "Select the pair of genes:",
				    choices = paste(df_pairs$V1,df_pairs$V2,sep = "_vs_"),
				    selected = "CDH11_vs_MET"),
			sliderInput("n", "Bins", 30, 100, 50)
			),
		 
	  # Main pabel for displaying output
	  mainPanel(
		  # Oupput: Expresion Histogram & density plot -----
		  plotOutput(outputId = "ExpHist")
		)
)

server <- function(input, output){
       output$ExpHist <- renderPlot({
				drawExpHist(unlist(strsplit(input$genePair, "_vs_")), input$n)
       		     		}
		      	)

}

# Call the app
shinyApp(ui,server)

