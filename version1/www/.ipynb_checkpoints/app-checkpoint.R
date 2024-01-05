library(shiny)
library(shinythemes)
library(shinyWidgets)
library(ggrepel)
library(ggplot2)
library(dplyr)

LOCALPATH <- '/Users/Joana/Joana/app/'
print(list.dirs())
print(list.files())

if (Sys.info()['user']=='vibflysleep'){
  root_dir = LOCALPATH
} else {
  root_dir = "."
}

genes <- as.vector(unique(read.csv(gzfile(file.path(root_dir, "data", "ct_mms_long.csv.gz")))$gene))
ct <- as.vector(read.csv(gzfile(file.path(root_dir, "data", "clusters.csv.gz")))$cluster)
ct2 <- as.vector(read.csv(gzfile(file.path(root_dir, "data", "ct2.csv.gz")))$cluster)

ui <- navbarPage("Fly Sleep Single Cell", theme = shinytheme("flatly"),
                 
                 tabPanel("about",
                          p("Explore expression of your gene of interest in your cluster of interest across 
                            (1) different circadian Zeitgeber times, 
                            (2) different points of the sleep/wake cycle and
                            (3) different degrees of sleep pressure.")),
                 
                 tabPanel("circadian cyclers",
                          fluidRow(sidebarPanel(
                            pickerInput("cell_type1",
                                        "Choose cell type(s)",
                                        choices = ct,
                                        selected = 'EG_1',
                                        multiple = TRUE,
                                        options = list(`actions-box` = TRUE)),
                            
                            selectInput("transcript1",
                                        label = "Choose a transcript",
                                        choices = genes,
                                        selected = genes[1]))),
                          
                          plotOutput("lineplot", height = "700px", width = "1500px")
                 ),
                 
                 tabPanel("sleep vs. wake",
                          fluidRow(sidebarPanel(
                            selectInput("cell_type2",
                                        label = "Choose a cell type",
                                        choices = ct,
                                        selected = 'EG_1'),
                            
                            selectInput("transcript2",
                                        label = "Choose a transcript",
                                        choices = genes,
                                        selected = genes[1]))),
                          
                          plotOutput("volcano", width="1250px", height = "600px")
                 ),
                 
                 tabPanel("sleep drive",
                          tabsetPanel(
                            tabPanel("heatmap",
                              selectInput("cell_type3",
                                          label = "Choose a cell type",
                                          choices = ct2,
                                          selected = 'EG_1'),
                              
                              plotOutput("heatmap", height = "1000px")),
                            
                            tabPanel("lineplot",
                              pickerInput("cell_type4",
                                          "Choose cell type(s)",
                                          choices = ct,
                                          selected = 'EG_1',
                                          multiple = TRUE,
                                          options = list(`actions-box` = TRUE)),
                              
                              selectInput("transcript3",
                                          label = "Choose a transcript",
                                          choices = genes,
                                          selected = genes[1]),
                              
                              plotOutput("lineplot2", height = "700px", width = "1500px"))
                            ))
                 )

server <- function(input, output) {

  zt_mms <- read.csv(gzfile(file.path(root_dir, "data/zt_mms_long.csv.gz")))
  setDT(zt_mms)
  
  output$lineplot <- renderPlot({
    
    ggplot(zt_mms[index == input$transcript1 & cluster %in% input$cell_type1], aes(x=ZT, y=expression, color=cluster)) +
      geom_line()+
      theme_classic(base_size=25)+
      ggtitle(input$transcript1)+
      ylab("mean normalized counts")+
      scale_x_continuous(breaks=c(2, 8, 14, 20))
    
  }, res=96)
  
  sleep_wake <- read.csv(gzfile(file.path(root_dir, "data/combined_nobatch.csv.gz"))) 
  setDT(sleep_wake)
  cols <- c("adj p < 1e-10" = "#e377c2", 
            "adj p < 0.01" = "#9467bd", 
            "adj p < 0.05" = "#c5b0d5",
            "p < 0.05" = "#7f7f7f",
            "p > 0.05" = "black")
  
  output$volcano <- renderPlot({
    
    ggplot(sleep_wake[cluster %in% input$cell_type2], aes(x=logfoldchange, y=-log10(pval), col=pvalue)) +
      geom_point() + 
      geom_text_repel(data=sleep_wake[cluster %in% input$cell_type2] %>%
                        filter(gene %in% input$transcript2),
                      aes(label=gene), size=10, max.overlaps = 100, box.padding = 2) +
      scale_x_continuous(limits = c(-3, 3), name="log2 fold change") +
      scale_y_continuous(limits = c(0, 10), name="-log10 p-value") +
      theme_classic()+
      theme(axis.text.x=element_text(size=25, colour = "black"),
            axis.title.x=element_text(size=25),
            axis.text.y=element_text(size=25, colour = "black"),
            axis.title.y=element_text(size=25),
            legend.text = element_text(size=25),
            legend.title = element_blank(),
            plot.title = element_text(size = 25, hjust=0.5)) +
      scale_color_manual(values=cols)+
      guides(color = guide_legend(override.aes = aes(label = "")))
    
  }, res=96)

  output$heatmap <- renderPlot({
    
    cluster_mms <- read.csv(gzfile(file.path(root_dir, "data", paste0(input$cell_type3, ".csv.gz"))), row.names = 'index')
    clock <- read.csv(gzfile("/Users/Joana/Joana/app/data/clock.csv.gz"), row.names = 'index')
    colnames(cluster_mms)[1] <- "gene"
    cluster_mms <- as.matrix(cluster_mms)
    #z-score
    cluster_mms_z <- t(apply(cluster_mms, 1, scale))
    
    heatmap(cluster_mms_z,
            Colv = NA)
  }, res=96)
  
  ct_mms <- read.csv(gzfile(file.path(root_dir, "data/ct_mms_long.csv.gz")))
  
  setDT(ct_mms)
  
  output$lineplot2 <- renderPlot({
    
    ggplot(ct_mms[gene == input$transcript3 & cluster %in% input$cell_type4], aes(x=condition, y=expression, color=cluster)) +
      geom_line()+
      theme_classic(base_size=25)+
      ylab("mean normalized counts")+
      xlab("conditions ordered from low to high sleep drive")+
      scale_x_continuous(breaks=c(1, 2, 3, 4, 5, 6, 7), labels=c("GBX-S ZT12-8", "S ZT12-20", "S ZT12-14", "SD ZT12-14", "SD ZT12-20", "SD ZT12-2", "SD ZT12-8"))+ 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  }, res=96)
}

shinyApp(ui = ui, server = server)

