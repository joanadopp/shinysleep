library(plotly)
library(shinythemes)
library(dplyr)
library(shiny)
library(RColorBrewer)
library(heatmaply)
library(tidyr)
library(R.utils)
library(ggplot2)
library(data.table)

LOCALPATH <- '/Users/Joana/Joana/shiny/version2/'

if (Sys.info()['user']=='vibflysleep'){
  root_dir = LOCALPATH
} else {
  root_dir = "."
}

gene1 <- as.vector(unique(read.csv(file.path(root_dir, "data", "clock_all_sign.csv"))$gene))
ct1 <- as.vector(unique(read.csv(file.path(root_dir, "data", "clock_all_sign.csv"))$cluster))
gene2 <- as.vector(unique(read.csv(file.path(root_dir, "data", "combined_nobatch_sign.csv"))$gene))
ct2 <- as.vector(unique(read.csv(file.path(root_dir, "data", "combined_nobatch_sign.csv"))$cluster))
gene3 <- as.vector(unique(read.csv(file.path(root_dir, "data", "sleep_drive_exp_r2_p_wide_sign.csv"))$gene))
ct3 <- as.vector(unique(read.csv(file.path(root_dir, "data", "sleep_drive_exp_r2_p_wide_sign.csv"))$cluster))

ui <- navbarPage("Fly Sleep Single Cell", theme = shinytheme("flatly"),
                 
                 tabPanel("about",
                          mainPanel(
                            tags$h1("Explore expression of your gene of interest in your cell type of interest", style = "font-size:18px"),
                            tags$h1("(1) across different circadian Zeitgeber times,", style = "font-size:18px"),
                            tags$h1("(2) between sleep and wake states and", style = "font-size:18px"),
                            tags$h1("(3) across different degrees of sleep pressure.", style = "font-size:18px"),
                            
                            tags$br(),
                            tags$br(),
                            
                            tags$h1("There are two versions of this shiny app.", style = "font-size:18px"),
                            div(tags$a(href="https://joana-dopp.shinyapps.io/Fly_Sleep_Single_Cell_v1/", "version 1"),
                                ": all genes in all annotated clusters", style = "font-size:18px"),
                            
                            tags$p("version 2 (this one): significant genes in all clusters", style = "font-size:18px"),
                            
                            tags$br(),
                            tags$br(),
                            
                            div("Our single cell atlas can be accessed on", style = "font-size:18px", 
                                tags$a(href="https://scope.aertslab.org/#/Fly_Brain_Sleep/Fly_Brain_Sleep%2FFly_Sleep.loom/welcome", "SCope.", style = "font-size:18px")),
                            
                            img(src='tsne.png', width = 750, height = 600, style = "margin-bottom: 50px; margin-top: 50px"),
                            
                            div("Citation: the associated preprint is available", style = "font-size:18px",
                                tags$a(href="https://www.biorxiv.org/content/10.1101/2023.03.22.533150v1", "here.", style = "font-size:18px")),
                            
                            div("This app was built by Joana Dopp. If you have questions/comments, please", style = "font-size:18px",
                                tags$a(href="mailto:joana.dopp@kuleuven.be", "get in touch.", style = "font-size:18px")),
                            
                            tags$br(),
                            tags$br()
                          )
                 ),
                 
                 tabPanel("circadian cyclers", 
                          
                          tags$h1("Hover over the points of the lineplot and you'll see adjusted p-values of the three replicates.", style = "font-size:18px"),
                          img(src='graphic_clock.png', width = 864, height = 219.2, style = "margin-top: 30px; margin-bottom: 30px"),
                          
                          fluidRow(sidebarPanel(
                            
                            selectizeInput("cell_type1",
                                           "Choose cell type(s)",
                                           choices = ct1,
                                           multiple = TRUE),
                            
                            selectizeInput("gene1",
                                           label = "Choose a gene",
                                           choices = NULL))),
                          
                          plotlyOutput("lineplot1", width = "1000px")
                 ),
                 
                 tabPanel("sleep vs. wake",
                          fluidRow(sidebarPanel(
                            
                            selectizeInput("cell_type2",
                                           label = "Choose a cell type",
                                           choices = ct2),
                            
                            selectizeInput("gene2",
                                           label = "Choose gene(s)",
                                           multiple = TRUE,
                                           choices = NULL))),
                          
                          mainPanel(tags$h1("Hover over the genes in the volcano plot and you'll see log2 fold change and adjusted p-value.", style = "font-size:18px"),
                                    img(src='wake_vs_sleep2.png', width = 509, height = 86, style = "margin-left: 180px;")),
                          
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          
                          plotlyOutput("volcanoPlot", width="1000px", height = "600px"),
                          
                          mainPanel(img(src='wake_vs_sleep.png', width = 1153.33, height = 251.11, style = "margin-top: 30px; margin-bottom: 30px"))
                 ),
                 
                 tabPanel("sleep drive",
                          tabsetPanel(
                            
                            tabPanel("compare cell types",
                                     
                                     tags$h1("Hover over the points of the lineplot and you'll see Pearson r2 and adjusted p-value.", style = "font-size:18px"),
                                     img(src='graphic_sleepdrive.png', width = 773.6, height = 194, style = "margin-top: 30px; margin-bottom: 30px"),
                                     
                                     selectizeInput("cell_type4",
                                                    "Choose cell type(s)",
                                                    choices = ct3,
                                                    multiple = TRUE),
                                     
                                     selectizeInput("gene3",
                                                    label = "Choose a gene",
                                                    choices = NULL),
                                     
                                     plotlyOutput("lineplot2", width = "1000px"),
                                     tags$p("GBX-S ZT12-8: Gaboxadol-induced sleep from ZT12 - ZT8 (20h)",
                                            tags$p("S ZT12-20: Sleep from ZT12 - ZT20 (8h)"), 
                                            tags$p("S ZT12-14: Sleep from ZT12 - ZT14 (2h)"), 
                                            tags$p("SD ZT12-14: Sleep deprivation from ZT12 - ZT14 (2h)"), 
                                            tags$p("SD ZT12-20: Sleep deprivation from ZT12 - ZT20 (8h)"), 
                                            tags$p("SD ZT12-2: Sleep deprivation from ZT12 to ZT2 (14h)"), 
                                            tags$p("SD ZT12-8: Sleep deprivation from ZT12 - ZT8 (20h)"), style = "margin-top: 30px;")
                            ),
                            
                            tabPanel("compare genes",
                                     
                                     tags$h1("Hover over the tiles of the heatmap and you'll see Pearson r2 and adjusted p-value (also color-coded in the last column).", style = "font-size:18px"),
                                     img(src='graphic_sleepdrive.png', width = 773.6, height = 194, style = "margin-top: 30px; margin-bottom: 30px"),
                                     
                                     selectizeInput("cell_type3",
                                                    label = "Choose a cell type",
                                                    choices = ct3),
                                     
                                     div(plotlyOutput("heatmap"), style = "margin-bottom:50px;")
                            )
                          )
                 )
)

server <- function(input, output, session) {
  
  observe({
    genes_filtered <- gene1[gene1 %in% unique(clock_long$gene[clock_long$cluster %in% input$cell_type1])]
    updateSelectizeInput(session, "gene1", choices = genes_filtered, server = FALSE)
  })
  
  clock_long <- read.csv(file.path(root_dir, "data/clock_all_sign.csv"))
  setDT(clock_long)
  
  output$lineplot1 <- renderPlotly({
    color_scale <- brewer.pal(9, "Set1")
    
    plot_ly(clock_long[gene == input$gene1 & cluster %in% input$cell_type1], x = ~ZT, y = ~expression, type = 'scatter', mode = 'lines',
            text = ~paste("adj. p-value (rep1): ", ifelse(is.na(adj.p_1), "NA", round(adj.p_1, 4)),
                          "<br>adj. p-value (rep2): ", ifelse(is.na(adj.p_2), "NA", round(adj.p_2, 4)),
                          "<br>adj. p-value (rep3): ", ifelse(is.na(adj.p_3), "NA", round(adj.p_3, 4))),
            line = list(color = ~cluster),
            name = ~cluster,
            width = 1000) %>%
      layout(title = input$gene1, showlegend=TRUE,
             xaxis = list(title = 'Zeitgeber time (ZT)',
                          tickvals = list(2,8,14,20),
                          tickfont = list(size = 18),
                          titlefont = list(size = 20)),
             yaxis = list(title = 'mean normalized counts',
                          titlefont = list(size = 20)))
  })
  
  updateSelectizeInput(session, "gene2", choices = gene2, server = TRUE)
  
  sleep_wake <- read.csv(file.path(root_dir, "data/combined_nobatch_sign.csv"))
  setDT(sleep_wake)
  
  output$volcanoPlot <- renderPlotly({
    
    key <- highlight_key(sleep_wake, ~gene, group ="gene")
    
    p <- ggplot(sleep_wake[cluster %in% input$cell_type2], aes(label=gene, label2=pval_adj, label3=pval)) +
      geom_point(aes(x=logfoldchange, y=-log10(pval), 
                     color = ifelse(gene %in% input$gene2, "selected genes", "all genes")), 
                 size = 0.5) +
      scale_color_manual(values = c("all genes" = "grey", "selected genes" = "red")) + 
      geom_text(data = sleep_wake[cluster %in% input$cell_type2] %>%
                  filter(gene %in% input$gene2),
                aes(label = gene, color = ifelse(gene %in% input$gene2, "selected genes", "all genes"),
                    x = logfoldchange, y = -log10(pval)), size = 3) +
      scale_x_continuous(limits = c(-3, 3), name = "log2 fold change") +
      scale_y_continuous(limits = c(0, 10), name = "-log10 p-value") +
      theme_classic(base_size = 15) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
      labs(color = "selection")
    
    ggplotly(p,
             tooltip = c("pval_adj", "gene", "pval", "logfoldchange", "y"), highlight = "plotly_selected")
  })
  
  sleep_drive_wide <- read.csv(file.path(root_dir, "data/sleep_drive_exp_r2_p_wide_sign.csv"))
  
  output$heatmap <- renderPlotly({
    
    ct_select <- subset(sleep_drive_wide, cluster == input$cell_type3)
    ct_sign <- as.data.frame(subset(ct_select, p.adj < 0.05, select=2:11))
    rownames(ct_sign) <- ct_sign$gene
    ct_sign$gene <- NULL
    ct_matrix_z <- as.matrix(t(apply(ct_sign[,1:7], 1, scale)))
    ct_final <- cbind(ct_matrix_z, ct_sign[, c("Pearson_r2", "p.adj")])
    colnames(ct_final)[1:ncol(ct_final)] <- c("GBX-S ZT12-8", "S ZT12-20", "S ZT12-14", "SD ZT12-14", "SD ZT12-20", "SD ZT12-2", "SD ZT12-8", "Pearson_r2", "p.adj")
    
    anno_mat <- ct_final %>%
      mutate_all(~ paste("adjusted p-value:", round(ct_final[, "p.adj"],4), "\nPearson r2:", round(ct_final[, "Pearson_r2"],4)))
    
    heatmap_plot <- heatmaply(
      x = ct_final[, 1:7],
      Rowv = TRUE,  
      Colv = FALSE, 
      hclust_method = "single",
      col = brewer.pal(9, "PiYG"),
      plot_method = "ggplot",
      #plotly::layout(showLegend = FALSE, annotations = list(visibile=FALSE)),
      custom_hovertext = anno_mat[, 1:7],
      row_side_palette = colorRampPalette(c('yellow','blue')),
      row_side_colors = round(ct_final[, "p.adj"], 4),
      showticklabels = c(F,T),
      width = 1000,
      height = 1500
    )
    
    heatmap_plot <- heatmap_plot %>%
      layout(xaxis = list(title = "conditions ordered from low to high sleep drive", tickfont = list(size = 20)),
             font = list (size = 18))
    
    heatmap_plot
  })
  
  observe({
    genes_filtered <- gene3[gene3 %in% unique(sleep_drive_long$gene[sleep_drive_long$cluster %in% input$cell_type3])]
    updateSelectizeInput(session, "gene3", choices = genes_filtered, server = FALSE)
  })
  
  sleep_drive_long <- sleep_drive_wide %>%
    pivot_longer(cols = starts_with("expression"),
                 names_to = c(".value", "condition"),
                 names_sep = "_")
  
  setDT(sleep_drive_long)
  
  output$lineplot2 <- renderPlotly({
    color_scale <- brewer.pal(9, "Set1")
    
    plot_ly(sleep_drive_long[gene == input$gene3 & cluster %in% input$cell_type4], x = ~condition, y = ~expression, type = 'scatter', mode = 'lines',
            text = ~paste("adj. p-value: ", round(p.adj, 4), "<br>Pearson r2: ", round(Pearson_r2, 4)),
            line = list(color = ~cluster),
            name = ~cluster,
            width = 1000) %>%
      layout(title = input$gene3, showlegend = TRUE,
             xaxis = list(title = 'conditions ordered from low to high sleep drive',
                          tickvals = 1:length(c("GBX-S ZT12-8", "S ZT12-20", "S ZT12-14", "SD ZT12-14", "SD ZT12-20", "SD ZT12-2", "SD ZT12-8")),
                          ticktext = c("GBX-S ZT12-8", "S ZT12-20", "S ZT12-14", "SD ZT12-14", "SD ZT12-20", "SD ZT12-2", "SD ZT12-8"),
                          tickfont = list(size = 18),
                          titlefont = list(size = 20)),
             yaxis = list(title = 'mean normalized counts',
                          titlefont = list(size = 20))
      )
  })
}

shinyApp(ui = ui, server = server)

