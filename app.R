setwd("C:/shiny2")

library(shiny)
library(DT)
library(ggplot2)
library(shinydashboard)
library(shinyBS)
library(shinythemes)
library(d3heatmap)
library(survival)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(gtable)
library(reshape)
library(sm)
library(vioplot)
library(e1071)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(ClassDiscovery)


header<- dashboardHeader(
  
  title = "GSMP: Tissue-specific, gene-specific methylation profile",
  titleWidth = 700
)


sidebar<- dashboardSidebar(
  width = 200,
  sidebarMenu(
    id="tabs",
    menuItem("Search by pathway", tabName = "pwlevel", icon = icon("hand-right", lib = "glyphicon")),
    menuItem("Search by gene", tabName = "genelevel", icon = icon("hand-right", lib = "glyphicon"))
  )
)

body<- dashboardBody(
  
  tabItems(
    
    
    tabItem(tabName = "pwlevel",  
            fluidRow(box(width=4,
                         selectInput(inputId = "cancer", label = "Cancer Type", choices = c("Colorectal cancer","Breast cancer","Pancreatic cancer"))
            ),box(width=4,
                  selectInput(inputId="checkpw", label="Pathway Categories", choices=c("Biocarta","Hallmark","KEGG","Oncogenic","Reactome","Genomic_Locus"))
            )),
            fluidRow(column(8,offset=2,DT::dataTableOutput('pathwaytable'))),
            fluidRow(actionButton("pw_infor",'GO!')),
            bsModal("result", "Pathway information","pw_infor",uiOutput("pw_result"),size="large")
    ),
    tabItem(tabName = "genelevel",  
            fluidRow(box(width=4,
                         selectInput(inputId = "cancer_by_gene", label = "Cancer Type", choices = c("Colorectal cancer","Breast cancer","Pancreatic cancer"))
            )),
            fluidRow(column(11,DT::dataTableOutput('genetable'))),
            fluidRow(actionButton("gw_infor",'GO!')),
            bsModal("result3", "GSMP information","gw_infor",uiOutput("gw_result"),size="large")
    )
    
  ))


ui<-dashboardPage(skin="blue",header,sidebar,body)


# load data
load("pathway.RData")
load("GSMP_brca.RData")
load("GSMP_crc.RData")
load("GSMP_paad.RData")
load("genelocation.RData")


#========== server function ===============
server<-function(input,output,session){
  
  
  output$pathwaytable = DT::renderDataTable({
    target_pw_table=result_pw_table[result_pw_table$pathway_type==input$checkpw,]
    target_pw_table
  },
  options=list(lengthMenu = c(20, 30, 50),autoWidth = FALSE,scrollX=TRUE,scrollY = "500px"), selection="single")
  
  
  observeEvent(input$pw_infor, {
    target_cancer=reactive({
      as.character(input$cancer)
    })
    
    target_pw=reactive({
      target_index <<-input$pathwaytable_rows_selected[length(input$pathwaytable_rows_selected)]
      target_index <<- as.numeric(target_index)
      
      target_pw_table=result_pw_table[result_pw_table$pathway_type==input$checkpw,]
      as.character(target_pw_table[target_index,1])})
    
    output$pw_result=renderUI({navbarPage(id="targetpw",
                                          title = target_pw(),
                                          tabPanel("Genes",DT::dataTableOutput('GSMP'),
                                                   fluidRow(actionButton("gene_info","GSMP")),
                                                   bsModal("result2", "GSMP information","gene_info",uiOutput("gene_result"),size="large")
                                          ),
                                          tabPanel("GSMS_heatmap",plotOutput("score_heatmap"),
                                                   downloadButton("download_GSMS_heatmap", "Download Plot"))
    )})
    
    
    target_genelist=reactive({
      genelist=all_DB_gene_list[[as.character(target_pw())]]
      
      if(target_cancer()=="Colorectal cancer")
      {
        genelist=as.character(genelist[genelist %in% names(GSMP_crc)])
        corr=rep(0,length(genelist))
        num=rep(0,length(genelist))
        for(i in 1:length(genelist))
        {
          if(!is.null(GSMP_crc[[genelist[i]]]$CpGs))
          {
            corr[i]=GSMP_crc[[genelist[i]]]$cv_TCGA
            num[i]=nrow(GSMP_crc[[genelist[i]]]$CpGs)
          }
        }
        genelist=data.frame(gene_symbol=as.character(genelist),corr_GSMS_exp=round(corr,3),num_CpGs=num)
      }
      
      else if(target_cancer()=="Breast cancer")
      {
        genelist=as.character(genelist[genelist %in% names(GSMP_brca)])
        corr=rep(0,length(genelist))
        num=rep(0,length(genelist))
        for(i in 1:length(genelist))
        {
          if(!is.null(GSMP_brca[[genelist[i]]]$CpGs))
          {
            corr[i]=GSMP_brca[[genelist[i]]]$cv_TCGA
            num[i]=nrow(GSMP_brca[[genelist[i]]]$CpGs)
          }
        }
        genelist=data.frame(gene_symbol=as.character(genelist),corr_GSMS_exp=round(corr,3),num_CpGs=num)
      }
      
      else 
      {
        genelist=as.character(genelist[genelist %in% names(GSMP_paad)])
        corr=rep(0,length(genelist))
        num=rep(0,length(genelist))
        for(i in 1:length(genelist))
        {
          if(!is.null(GSMP_paad[[genelist[i]]]$CpGs))
          {
            corr[i]=GSMP_paad[[genelist[i]]]$cv_TCGA
            num[i]=nrow(GSMP_paad[[genelist[i]]]$CpGs)
          }
        }
        genelist=data.frame(gene_symbol=as.character(genelist),corr_GSMS_exp=round(corr,3),num_CpGs=num)
      }
      
      genelist
    })
    
    
    plot_GSMS <- function(){
      genelist=target_genelist()
      
      if(nrow(genelist)<20)
      {
        font <- 12
      }
      
      else if(nrow(genelist)<40)
      {
        font <- 6
      }
      
      else
      {
        font <- 3
      }
      
      if(target_cancer()=="Colorectal cancer")
      {
        GSMS_matrix <- matrix(0,nrow=nrow(genelist),ncol=369)
        
        for(i in 1:nrow(genelist))
        {
          if(!is.null(GSMP_crc[[as.character(genelist[i,1])]]$meth_score_TCGA))
          {
            tmp_score <- GSMP_crc[[as.character(genelist[i,1])]]$meth_score_TCGA[,"scaled2"]
            GSMS_matrix[i,] <- tmp_score
          }
        }
        
        if(sum(genelist$corr_GSMS_exp>=0.4)>=3)
        {
          GSMS_matrix0 <- GSMS_matrix[genelist$corr_GSMS_exp>=0.4,]
        }
        else
        {
          GSMS_matrix0 <- GSMS_matrix
        }
        
        tmp1 <- hclust(distanceMatrix(GSMS_matrix0[, CMS=="CMS1"], 'pearson'), method="ward.D")
        tmp.idx1 = which(CMS=="CMS1")[tmp1$order]
        
        tmp2 <- hclust(distanceMatrix(GSMS_matrix0[, CMS=="CMS2"], 'pearson'), method="ward.D")
        tmp.idx2 = which(CMS=="CMS2")[tmp2$order]
        
        tmp3 <- hclust(distanceMatrix(GSMS_matrix0[, CMS=="CMS3"], 'pearson'), method="ward.D")
        tmp.idx3 = which(CMS=="CMS3")[tmp3$order]
        
        tmp4 <- hclust(distanceMatrix(GSMS_matrix0[, CMS=="CMS4"], 'pearson'), method="ward.D")
        tmp.idx4 = which(CMS=="CMS4")[tmp4$order]
        
        tmp5 <- hclust(distanceMatrix(GSMS_matrix0[, CMS=="Indet"], 'pearson'), method="ward.D")
        tmp.idx5 = which(CMS=="Indet")[tmp5$order]
        
        GSMS_matrix <- GSMS_matrix[,c(tmp.idx1,tmp.idx2,tmp.idx3,tmp.idx4,tmp.idx5)]
        rownames(GSMS_matrix) <- as.character(genelist$gene_symbol)
        
        ha_column = HeatmapAnnotation(df = data.frame(CMS = c(rep("CMS1",length(tmp.idx1)), rep("CMS2",length(tmp.idx2)),
                                                              rep("CMS3",length(tmp.idx3)), rep("CMS4",length(tmp.idx4)),
                                                              rep("Indet",length(tmp.idx5)))), 
                                      col = list(CMS = c("CMS1"="orange","CMS2"="blue","CMS3"="palevioletred1","CMS4"="green","Indet"="gray")))
        
        GSMS_matrix <- data.frame(GSMS_matrix)
        GSMS_matrix$group <- ifelse(genelist$corr_GSMS_exp>=0.4,"greater than 0.4","less than 0.4")
        
        Heatmap(GSMS_matrix[,-ncol(GSMS_matrix)],name="GSMS",split=GSMS_matrix$group,show_column_names=FALSE, cluster_columns=FALSE,
                column_title=paste("TCGA",target_cancer(),"samples",sep=" "),column_title_side="bottom",
                row_title=target_pw(), row_title_side="left", row_names_gp = gpar(fontsize = font), top_annotation=ha_column,col = colorRamp2(c(-2, 0, 2), c("red", "black", "blue"))) 
      }
      
      else if(target_cancer()=="Breast cancer")
      {
        GSMS_matrix <- matrix(0,nrow=nrow(genelist),ncol=782)
        
        for(i in 1:nrow(genelist))
        {
          if(!is.null(GSMP_brca[[as.character(genelist[i,1])]]$meth_score_TCGA))
          {
            tmp_score <- GSMP_brca[[as.character(genelist[i,1])]]$meth_score_TCGA[,"scaled2"]
            GSMS_matrix[i,] <- tmp_score
          }
        }
        
        if(sum(genelist$corr_GSMS_exp>=0.4)>=3)
        {
          GSMS_matrix0 <- GSMS_matrix[genelist$corr_GSMS_exp>=0.4,]
        }
        else
        {
          GSMS_matrix0 <- GSMS_matrix
        }
        
        tmp <- hclust(distanceMatrix(GSMS_matrix0, 'pearson'), method="ward.D")
        
        GSMS_matrix <- GSMS_matrix[,tmp$order]
        rownames(GSMS_matrix) <- as.character(genelist$gene_symbol)
        
        GSMS_matrix <- data.frame(GSMS_matrix)
        GSMS_matrix$group <- ifelse(genelist$corr_GSMS_exp>=0.4,"greater than 0.4","less than 0.4")
        
        Heatmap(GSMS_matrix[,-ncol(GSMS_matrix)],name="GSMS",split=GSMS_matrix$group,show_column_names=FALSE, cluster_columns=FALSE,
                column_title=paste("TCGA",target_cancer(),"samples",sep=" "),column_title_side="bottom",
                row_title=target_pw(), row_title_side="left", row_names_gp = gpar(fontsize = font),col = colorRamp2(c(-2, 0, 2), c("red", "black", "blue")))
      }
      
      else if(target_cancer()=="Pancreatic cancer")
      {
        GSMS_matrix <- matrix(0,nrow=nrow(genelist),ncol=178)
        
        for(i in 1:nrow(genelist))
        {
          if(!is.null(GSMP_paad[[as.character(genelist[i,1])]]$meth_score_TCGA))
          {
            tmp_score <- GSMP_paad[[as.character(genelist[i,1])]]$meth_score_TCGA[,"scaled2"]
            GSMS_matrix[i,] <- tmp_score
          }
        }
        
        if(sum(genelist$corr_GSMS_exp>=0.4)>=3)
        {
          GSMS_matrix0 <- GSMS_matrix[genelist$corr_GSMS_exp>=0.4,]
        }
        else
        {
          GSMS_matrix0 <- GSMS_matrix
        }
        
        tmp <- hclust(distanceMatrix(GSMS_matrix0, 'pearson'), method="ward.D")
        
        GSMS_matrix <- GSMS_matrix[,tmp$order]
        rownames(GSMS_matrix) <- as.character(genelist$gene_symbol)
        
        GSMS_matrix <- data.frame(GSMS_matrix)
        GSMS_matrix$group <- ifelse(genelist$corr_GSMS_exp>=0.4,"greater than 0.4","less than 0.4")
        
        Heatmap(GSMS_matrix[,-ncol(GSMS_matrix)],name="GSMS",split=GSMS_matrix$group,show_column_names=FALSE, cluster_columns=FALSE,
                column_title=paste("TCGA",target_cancer(),"samples",sep=" "),column_title_side="bottom",
                row_title=target_pw(), row_title_side="left", row_names_gp = gpar(fontsize = font),col = colorRamp2(c(-2, 0, 2), c("red", "black", "blue")))
      }
      
    }
    
    output$score_heatmap=renderPlot({
      plot_GSMS()
    })
    
    output$download_GSMS_heatmap <- downloadHandler(
      filename = function() {
        paste0("GSMS_",target_cancer(),"_",target_pw(),".png")
      },
      content = function(file) {
        png(file)
        print(plot_GSMS())
        dev.off()
      },contentType = 'image/png'
    )
    
    output$GSMP = DT::renderDataTable({
      GSMP_table=target_genelist()
      GSMP_table
    },
    options=list(lengthMenu = c(10, 20, 30),autoWidth = FALSE,scrollX=TRUE, scrollY = "500px"), selection="single")
    
    
    observeEvent(input$gene_info, {
      
      target_gene=reactive({
        target_g_index <<-input$GSMP_rows_selected[length(input$GSMP_rows_selected)]
        target_g_index <<- as.numeric(target_g_index)
        
        GSMP_table=target_genelist()
        
        as.character(GSMP_table[target_g_index,1])})
      
      
      output$gene_result=renderUI({navbarPage(id="targetgene",
                                              title = target_gene(),
                                              tabPanel("GSMP_plot",plotOutput("GSMPplot"),
                                                       downloadButton("download_GSMP_plot", "Download Plot")),
                                              tabPanel("GSMP_table",DT::dataTableOutput('GSMPtable',width=500),
                                                       tags$style(type='text/css', "#download_factor {float:right}"),
                                                       downloadButton("download_GSMP_table", "Download Table"))
      )})
      
      
      plot_GSMP <- function()
      {
        gene <- as.character(target_gene())
        
        if(target_cancer()=="Colorectal cancer")
        {
          tmpdata <- GSMP_crc[[gene]]
        }
        
        else if(target_cancer()=="Breast cancer")
        {
          tmpdata <- GSMP_brca[[gene]]
        }
        
        else
        {
          tmpdata <- GSMP_paad[[gene]]
        }
        
        
        info <- tmpdata$CpGs  
        
        startbp <- tmpdata$startbp
        endbp <- tmpdata$endbp
        chr <- as.character(tmpdata$chr)
        
        x1 <- pmax(0,pmin(startbp,endbp)/1000000 - 0.5)
        x2 <- pmax(startbp,endbp)/1000000 + 0.5
        
        if(!is.null(info))
        {
          newinfo <- info[order(info$bpLocation),]
          
          y1 <- min(info$standard_coef)-0.05
          y2 <- max(info$standard_coef)+0.05
          
          plot(newinfo$bpLocation/1000000,newinfo$standard_coef,type="n",xlim=c(x1,x2),ylim=c(y1,y2),xlab="Mbp",ylab="standardized coef",main=gene,cex.lab=1)
        }
        
        else
        {
          plot(0,0,type="n",xlim=c(x1,x2),ylim=c(-1,1),xlab="Mbp",ylab="standardized coef",main=gene,cex.lab=1)
        }
        
        idx1 <- (data$start >= par("usr")[1]*1000000) & (data$start <= par("usr")[2]*1000000)
        idx2 <- (data$end >= par("usr")[1]*1000000) & (data$end <= par("usr")[2]*1000000)
        idx3 <- data$chr == chr
        
        tmp <- data[(idx1 | idx2) & idx3 ,]
        
        if(is.null(info) | nrow(tmp) <= 10)
        {
          for(k in 1:nrow(tmp))
          {
            polygon(c(tmp[k,"start"]/1000000,tmp[k,"end"]/1000000,tmp[k,"end"]/1000000,tmp[k,"start"]/1000000),c(-10,-10,10,10),border=NA,col="gray94")
            
            loc <- pmax(tmp[k,"start"],par("usr")[1]*1000000)
            loc <- pmin(loc,par("usr")[2]*1000000)
            
            if(loc!=startbp)
            {
              axis(side = 3, at = loc/1000000,labels = tmp[k,"gene"],las=2,cex.axis=0.7,col="blue")
            }
          }
        }
        
        else
        {
          best_loc <- info$bpLocation
          
          for(k in 1:nrow(tmp))
          {
            start_end <- range(c(tmp[k,"start"],tmp[k,"end"]))
            
            check <- unlist(lapply(best_loc,function(x){ (x >= start_end[1] & x <= start_end[2]) | abs(x - start_end[1]) < 2000 | abs(x - start_end[2]) < 2000}))
            
            if(sum(check)!=0)
            {
              polygon(c(tmp[k,"start"]/1000000,tmp[k,"end"]/1000000,tmp[k,"end"]/1000000,tmp[k,"start"]/1000000),c(-10,-10,10,10),border=NA,col="gray94")
              
              loc <- pmax(tmp[k,"start"],par("usr")[1]*1000000)
              loc <- pmin(loc,par("usr")[2]*1000000)
              
              if(loc!=startbp)
              {
                axis(side = 3, at = loc/1000000,labels = tmp[k,"gene"],las=2,cex.axis=0.7,col="blue")
              }
            }
          }
        }
        
        polygon(c(startbp/1000000,endbp/1000000,endbp/1000000,startbp/1000000),c(-10,-10,10,10),border=NA,col="gray80")
        
        abline(h=0,col="lightgray")
        
        if(!is.null(info))
        {
          points(info$bpLocation/1000000,info$standard_coef,pch=21,bg="red",cex=1)
          
          rug(info$bpLocation[info$cgi=="OpenSea"]/1000000,side=3,lwd=1,col="black")
          rug(info$bpLocation[info$cgi=="Shelf"]/1000000,side=3,lwd=1,col="green")
          rug(info$bpLocation[info$cgi=="Shore"]/1000000,side=3,lwd=1,col="hotpink")
          rug(info$bpLocation[info$cgi=="Island"]/1000000,side=3,lwd=1,col="red")
          
          arrows(startbp/1000000,y2+0.02,endbp/1000000,y2+0.02,length=0.1,xpd=TRUE)
          axis(side=3,at=c(startbp/1000000,endbp/1000000),labels=c("5\'","3\'"),tick=FALSE)
          
          legend("topleft",text.col=c("black","red","pink","green","black"),
                 legend=c("top marks","CpG Island","CpG Shore","CpG Shelf","Open Sea"),bty="n",cex=0.7)
          
          if(target_cancer()=="Breast cancer")
          {   
            rug(info$bpLocation[grep("TssA",info$chr_state_breast)]/1000000,lwd=1,col="red")
            rug(info$bpLocation[grep("TssAFlnk",info$chr_state_breast)]/1000000,lwd=1,col="orangered")
            rug(info$bpLocation[grep("7_Enh",info$chr_state_breast)]/1000000,lwd=1,col="orange")
            rug(info$bpLocation[grep("EnhG",info$chr_state_breast)]/1000000,lwd=1,col="greenyellow")
            
            legend("bottomleft",text.col=c("black","orange","greenyellow","orangered","red"),
                   legend=c("bottom marks","Enhancer","Genic Enhancer","Flanking Active TSS","Active TSS"),bty="n",cex=0.7)
          }
          
          else if(target_cancer()=="Pancreatic cancer")
          {   
            rug(info$bpLocation[grep("TssA",info$chr_state_pancreas)]/1000000,lwd=1,col="red")
            rug(info$bpLocation[grep("TssFlnk",info$chr_state_pancreas)]/1000000,lwd=1,col="orangered")
            rug(info$bpLocation[grep("EnhA",info$chr_state_pancreas)]/1000000,lwd=1,col="orange")
            rug(info$bpLocation[grep("EnhG",info$chr_state_pancreas)]/1000000,lwd=1,col="greenyellow")
            
            legend("bottomleft",text.col=c("black","orange","greenyellow","orangered","red"),
                   legend=c("bottom marks","Active Enhancer","Genic Enhancer","Flanking TSS","Active TSS"),bty="n",cex=0.7)
          }
          
          else if(target_cancer()=="Colorectal cancer")
          {
            rug(info$bpLocation[info$active_promoter_enhancer==1]/1000000,lwd=1,col="red")
            rug(info$bpLocation[info$active_enhancer==1]/1000000+0.002,lwd=1,col="orange")
            rug(info$bpLocation[info$transcribed_enhancer==1]/1000000-0.002,lwd=1,col="greenyellow")
            
            legend("bottomleft",text.col=c("black","orange","greenyellow","red"),
                   legend=c("bottom marks","Active Enhancer","Transcribed Enhancer","Active Promoter/Enhancer"),bty="n",cex=0.7)
          }
          
          legend("topright",legend=c(paste("cv_TCGA:",tmpdata$cv_TCGA,sep=" ")),bty="n",cex=0.8)
          
        }
        
        else
        {
          legend("topright",legend=c("No CpG is selected for this gene!"),bty="n",cex=1)
        }
        
      }
      
      
      output$GSMPplot=renderPlot({
        plot_GSMP()
      })
      
      output$download_GSMP_plot <- downloadHandler(
        filename = function(){
          paste0("GSMP_",target_cancer(),"_",target_gene(),".png")
        }, 
        content = function(file) {
          png(file)
          print(plot_GSMP())
          dev.off()
        },contentType = 'image/png'
      )
      
      
      output$GSMPtable=DT::renderDataTable({
        if(target_cancer()=="Colorectal cancer")
        {
          info <- GSMP_crc[[as.character(target_gene())]]$CpGs
        }
        else if(target_cancer()=="Breast cancer")
        {
          info <- GSMP_brca[[as.character(target_gene())]]$CpGs
        }
        else if(target_cancer()=="Pancreatic cancer")
        {
          info <- GSMP_paad[[as.character(target_gene())]]$CpGs
        }
        info$relative_dist <- round(info$relative_dist,3)
        info$cor_expression <- round(info$cor_expression,3)
        info$prior_prob <- round(info$prior_prob,3)
        info$standard_coef <- round(info$standard_coef,3)
        info$methsd <- round(info$methsd,3)
        info$unstandard_coef <- round(info$unstandard_coef,3)
        info
      },options = list(lengthMenu = c(5, 10, 20),autoWidth = FALSE,scrollX=TRUE, scrollY = "500px"),selection="single",rownames=FALSE)
      
      
      table_CpG <- reactive({
        if(target_cancer()=="Colorectal cancer")
        {
          info <- GSMP_crc[[as.character(target_gene())]]$CpGs
        }
        else if(target_cancer()=="Breast cancer")
        {
          info <- GSMP_brca[[as.character(target_gene())]]$CpGs
        }
        else if(target_cancer()=="Pancreatic cancer")
        {
          info <- GSMP_paad[[as.character(target_gene())]]$CpGs
        }
        info$relative_dist <- round(info$relative_dist,3)
        info$cor_expression <- round(info$cor_expression,3)
        info$prior_prob <- round(info$prior_prob,3)
        info$standard_coef <- round(info$standard_coef,3)
        info$methsd <- round(info$methsd,3)
        info$unstandard_coef <- round(info$unstandard_coef,3)
        info
      })
      
      # Downloadable csv of selected dataset ----
      output$download_GSMP_table <- downloadHandler(
        filename = function() {
          paste0("GSMP_", target_cancer(),"_", target_gene(),".csv")
        },
        content = function(file) {
          write.csv(table_CpG(), file, row.names = FALSE)
        }
      )
      
    })
    
  })
  

    
  target_cancer_by_gene=reactive({
    as.character(input$cancer_by_gene)
  })
  
  
  target_genelist_by_gene=reactive({
    if(target_cancer_by_gene()=="Colorectal cancer")
    {
      genelist=names(GSMP_crc)
      corr=rep(0,length(genelist))
      num=rep(0,length(genelist))
      for(i in 1:length(genelist))
      {
        if(!is.null(GSMP_crc[[i]]$CpGs))
        {
          corr[i]=GSMP_crc[[i]]$cv_TCGA
          num[i]=nrow(GSMP_crc[[i]]$CpGs)
        }
      }
      genelist=data.frame(gene_symbol=as.character(genelist),corr_GSMS_exp=round(corr,3),num_CpGs=num)
    }
    
    else if(target_cancer_by_gene()=="Breast cancer")
    {
      genelist=names(GSMP_brca)
      corr=rep(0,length(genelist))
      num=rep(0,length(genelist))
      for(i in 1:length(genelist))
      {
        if(!is.null(GSMP_brca[[i]]$CpGs))
        {
          corr[i]=GSMP_brca[[i]]$cv_TCGA
          num[i]=nrow(GSMP_brca[[i]]$CpGs)
        }
      }
      genelist=data.frame(gene_symbol=as.character(genelist),corr_GSMS_exp=round(corr,3),num_CpGs=num)
    }
    
    else 
    {
      genelist=names(GSMP_paad)
      corr=rep(0,length(genelist))
      num=rep(0,length(genelist))
      for(i in 1:length(genelist))
      {
        if(!is.null(GSMP_paad[[i]]$CpGs))
        {
          corr[i]=GSMP_paad[[i]]$cv_TCGA
          num[i]=nrow(GSMP_paad[[i]]$CpGs)
        }
      }
      genelist=data.frame(gene_symbol=as.character(genelist),corr_GSMS_exp=round(corr,3),num_CpGs=num)
    }
    
    genelist
    
  })
  
  
  output$genetable = DT::renderDataTable({
     target_gene_table=target_genelist_by_gene()
     target_gene_table
  },
  options=list(lengthMenu = c(30, 50, 80),autoWidth = FALSE,scrollX=TRUE,scrollY = "600px"), selection="single")
  
  
  observeEvent(input$gw_infor, {
    
    target_gw=reactive({
      target_gw_index <<-input$genetable_rows_selected[length(input$genetable_rows_selected)]
      target_gw_index <<- as.numeric(target_gw_index)
      
      target_gene_table=target_genelist_by_gene()
      as.character(target_gene_table[target_gw_index,1])})
    
    output$gw_result=renderUI({navbarPage(id="targetgene2",
                                            title = target_gw(),
                                            tabPanel("GSMP_plot",plotOutput("GSMPplot2"),
                                                     downloadButton("download_GSMP_plot2", "Download Plot")),
                                            tabPanel("GSMP_table",DT::dataTableOutput('GSMPtable2',width=500),
                                                     tags$style(type='text/css', "#download_factor {float:right}"),
                                                     downloadButton("download_GSMP_table2", "Download Table"))
    )})
    
    plot_GSMP_by_gene <- function()
    {
      gene <- as.character(target_gw())
      
      if(target_cancer_by_gene()=="Colorectal cancer")
      {
        tmpdata <- GSMP_crc[[gene]]
      }
      
      else if(target_cancer_by_gene()=="Breast cancer")
      {
        tmpdata <- GSMP_brca[[gene]]
      }
      
      else
      {
        tmpdata <- GSMP_paad[[gene]]
      }
      
      
      info <- tmpdata$CpGs  
      
      startbp <- tmpdata$startbp
      endbp <- tmpdata$endbp
      chr <- as.character(tmpdata$chr)
      
      x1 <- pmax(0,pmin(startbp,endbp)/1000000 - 0.5)
      x2 <- pmax(startbp,endbp)/1000000 + 0.5
      
      if(!is.null(info))
      {
        newinfo <- info[order(info$bpLocation),]
        
        y1 <- min(info$standard_coef)-0.05
        y2 <- max(info$standard_coef)+0.05
        
        plot(newinfo$bpLocation/1000000,newinfo$standard_coef,type="n",xlim=c(x1,x2),ylim=c(y1,y2),xlab="Mbp",ylab="standardized coef",main=gene,cex.lab=1)
      }
      
      else
      {
        plot(0,0,type="n",xlim=c(x1,x2),ylim=c(-1,1),xlab="Mbp",ylab="standardized coef",main=gene,cex.lab=1)
      }
      
      idx1 <- (data$start >= par("usr")[1]*1000000) & (data$start <= par("usr")[2]*1000000)
      idx2 <- (data$end >= par("usr")[1]*1000000) & (data$end <= par("usr")[2]*1000000)
      idx3 <- data$chr == chr
      
      tmp <- data[(idx1 | idx2) & idx3 ,]
      
      if(is.null(info) | nrow(tmp) <= 10)
      {
        for(k in 1:nrow(tmp))
        {
          polygon(c(tmp[k,"start"]/1000000,tmp[k,"end"]/1000000,tmp[k,"end"]/1000000,tmp[k,"start"]/1000000),c(-10,-10,10,10),border=NA,col="gray94")
          
          loc <- pmax(tmp[k,"start"],par("usr")[1]*1000000)
          loc <- pmin(loc,par("usr")[2]*1000000)
          
          if(loc!=startbp)
          {
            axis(side = 3, at = loc/1000000,labels = tmp[k,"gene"],las=2,cex.axis=0.7,col="blue")
          }
        }
      }
      
      else
      {
        best_loc <- info$bpLocation
        
        for(k in 1:nrow(tmp))
        {
          start_end <- range(c(tmp[k,"start"],tmp[k,"end"]))
          
          check <- unlist(lapply(best_loc,function(x){ (x >= start_end[1] & x <= start_end[2]) | abs(x - start_end[1]) < 2000 | abs(x - start_end[2]) < 2000}))
          
          if(sum(check)!=0)
          {
            polygon(c(tmp[k,"start"]/1000000,tmp[k,"end"]/1000000,tmp[k,"end"]/1000000,tmp[k,"start"]/1000000),c(-10,-10,10,10),border=NA,col="gray94")
            
            loc <- pmax(tmp[k,"start"],par("usr")[1]*1000000)
            loc <- pmin(loc,par("usr")[2]*1000000)
            
            if(loc!=startbp)
            {
              axis(side = 3, at = loc/1000000,labels = tmp[k,"gene"],las=2,cex.axis=0.7,col="blue")
            }
          }
        }
      }
      
      polygon(c(startbp/1000000,endbp/1000000,endbp/1000000,startbp/1000000),c(-10,-10,10,10),border=NA,col="gray80")
      
      abline(h=0,col="lightgray")
      
      if(!is.null(info))
      {
        points(info$bpLocation/1000000,info$standard_coef,pch=21,bg="red",cex=1)
        
        rug(info$bpLocation[info$cgi=="OpenSea"]/1000000,side=3,lwd=1,col="black")
        rug(info$bpLocation[info$cgi=="Shelf"]/1000000,side=3,lwd=1,col="green")
        rug(info$bpLocation[info$cgi=="Shore"]/1000000,side=3,lwd=1,col="hotpink")
        rug(info$bpLocation[info$cgi=="Island"]/1000000,side=3,lwd=1,col="red")
        
        arrows(startbp/1000000,y2+0.02,endbp/1000000,y2+0.02,length=0.1,xpd=TRUE)
        axis(side=3,at=c(startbp/1000000,endbp/1000000),labels=c("5\'","3\'"),tick=FALSE)
        
        legend("topleft",text.col=c("black","red","pink","green","black"),
               legend=c("top marks","CpG Island","CpG Shore","CpG Shelf","Open Sea"),bty="n",cex=0.7)
        
        if(target_cancer_by_gene()=="Breast cancer")
        {   
          rug(info$bpLocation[grep("TssA",info$chr_state_breast)]/1000000,lwd=1,col="red")
          rug(info$bpLocation[grep("TssAFlnk",info$chr_state_breast)]/1000000,lwd=1,col="orangered")
          rug(info$bpLocation[grep("7_Enh",info$chr_state_breast)]/1000000,lwd=1,col="orange")
          rug(info$bpLocation[grep("EnhG",info$chr_state_breast)]/1000000,lwd=1,col="greenyellow")
          
          legend("bottomleft",text.col=c("black","orange","greenyellow","orangered","red"),
                 legend=c("bottom marks","Enhancer","Genic Enhancer","Flanking Active TSS","Active TSS"),bty="n",cex=0.7)
        }
        
        else if(target_cancer_by_gene()=="Pancreatic cancer")
        {   
          rug(info$bpLocation[grep("TssA",info$chr_state_pancreas)]/1000000,lwd=1,col="red")
          rug(info$bpLocation[grep("TssFlnk",info$chr_state_pancreas)]/1000000,lwd=1,col="orangered")
          rug(info$bpLocation[grep("EnhA",info$chr_state_pancreas)]/1000000,lwd=1,col="orange")
          rug(info$bpLocation[grep("EnhG",info$chr_state_pancreas)]/1000000,lwd=1,col="greenyellow")
          
          legend("bottomleft",text.col=c("black","orange","greenyellow","orangered","red"),
                 legend=c("bottom marks","Active Enhancer","Genic Enhancer","Flanking TSS","Active TSS"),bty="n",cex=0.7)
        }
        
        else if(target_cancer_by_gene()=="Colorectal cancer")
        {
          rug(info$bpLocation[info$active_promoter_enhancer==1]/1000000,lwd=1,col="red")
          rug(info$bpLocation[info$active_enhancer==1]/1000000+0.002,lwd=1,col="orange")
          rug(info$bpLocation[info$transcribed_enhancer==1]/1000000-0.002,lwd=1,col="greenyellow")
          
          legend("bottomleft",text.col=c("black","orange","greenyellow","red"),
                 legend=c("bottom marks","Active Enhancer","Transcribed Enhancer","Active Promoter/Enhancer"),bty="n",cex=0.7)
        }
        
        legend("topright",legend=c(paste("cv_TCGA:",tmpdata$cv_TCGA,sep=" ")),bty="n",cex=0.8)
        
      }
      
      else
      {
        legend("topright",legend=c("No CpG is selected for this gene!"),bty="n",cex=1)
      }
      
    }
    
    
    output$GSMPplot2=renderPlot({
      plot_GSMP_by_gene()
    })
    
    output$download_GSMP_plot2 <- downloadHandler(
      filename = function(){
        paste0("GSMP_",target_cancer_by_gene(),"_",target_gw(),".png")
      }, 
      content = function(file) {
        png(file)
        print(plot_GSMP_by_gene())
        dev.off()
      },contentType = 'image/png'
    )
    
    output$GSMPtable2=DT::renderDataTable({
      if(target_cancer_by_gene()=="Colorectal cancer")
      {
        info <- GSMP_crc[[as.character(target_gw())]]$CpGs
      }
      else if(target_cancer_by_gene()=="Breast cancer")
      {
        info <- GSMP_brca[[as.character(target_gw())]]$CpGs
      }
      else if(target_cancer_by_gene()=="Pancreatic cancer")
      {
        info <- GSMP_paad[[as.character(target_gw())]]$CpGs
      }
      info$relative_dist <- round(info$relative_dist,3)
      info$cor_expression <- round(info$cor_expression,3)
      info$prior_prob <- round(info$prior_prob,3)
      info$standard_coef <- round(info$standard_coef,3)
      info$methsd <- round(info$methsd,3)
      info$unstandard_coef <- round(info$unstandard_coef,3)
      info
    },options = list(lengthMenu = c(5, 10, 20),autoWidth = FALSE,scrollX=TRUE, scrollY = "500px"),selection="single",rownames=FALSE)
    
    
    table_CpG_by_gene <- reactive({
      if(target_cancer_by_gene()=="Colorectal cancer")
      {
        info <- GSMP_crc[[as.character(target_gw())]]$CpGs
      }
      else if(target_cancer_by_gene()=="Breast cancer")
      {
        info <- GSMP_brca[[as.character(target_gw())]]$CpGs
      }
      else if(target_cancer_by_gene()=="Pancreatic cancer")
      {
        info <- GSMP_paad[[as.character(target_gw())]]$CpGs
      }
      info$relative_dist <- round(info$relative_dist,3)
      info$cor_expression <- round(info$cor_expression,3)
      info$prior_prob <- round(info$prior_prob,3)
      info$standard_coef <- round(info$standard_coef,3)
      info$methsd <- round(info$methsd,3)
      info$unstandard_coef <- round(info$unstandard_coef,3)
      info
    })
    
    # Downloadable csv of selected dataset ----
    output$download_GSMP_table2 <- downloadHandler(
      filename = function() {
        paste0("GSMP_", target_cancer_by_gene(),"_", target_gw(),".csv")
      },
      content = function(file) {
        write.csv(table_CpG_by_gene(), file, row.names = FALSE)
      }
    )
    
  })

  
}





shinyApp(ui=ui,server=server)





