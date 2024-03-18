#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)
library(DT)
library(ggplot2)
library(viridis)
library(tidyr)
library(scales)
library(dplyr)
library(stringr)

js <- '
$(document).on("plotly_hover", function(e, data) {
  Shiny.setInputValue("hover", data.points[0].pointNumber);
});
'

ui <- fluidPage(
  tags$head(tags$script(HTML(js))),
  h1("KIF1C-bioID Data"),
  a("Original publication", href='https://doi.org/10.1083/jcb.201812170', target = "_blank"),
  # titlePanel("KIF1C bioID"),
  fluidRow(
    plotlyOutput("scatterPlot", height = 500)
  ),
  fluidRow(
    textInput("search", "Search gene:", placeholder = "kif1c")
  ),
 
  # fluidRow(
  #   DTOutput("infoTable1")
  # ),
  fluidRow(
    DTOutput("infoTable2")
  )
)

server <- function(input, output) {
  
  # setwd("/Users/qgeng/Desktop/shiny/KIF1C_bioID")
  data <- read.csv("./tidy_bioID.csv", header = TRUE, fileEncoding="UTF-8-BOM")
  df <- data # df is for numeric plotting, data is for displaying
  
  data$Enrichment[data$Enrichment == 0] <- "Not present in control"
  df$Enrichment[df$Enrichment == 0] <- 1500
  names(data)[names(data) == 'pValue'] <- 'p.value'
  
  # Hook3_KIF1C <- df %>% filter(Genes %in% c("KIF1C", "HOOK3"))
  # RBP_subset <- df %>% filter(Genes %in% c("CASC3", "YTHD3", "TUT7", "TDRD3", "CTIF",
  #                                            "MEX3A", "TDRD7", "TOP3B", "IF4E2", "CCDC9",
  #                                            "GPBP1", "G3P", "E2AK2", "IF4A3", "SMG1", 
  #                                            "YTHD1", "RENT1", "PRC2A", "DDX3X", "ATX2", 
  #                                            "1433E", "STAU2", "FUBP3", "YTHD2"))
  
  
  output$scatterPlot <- renderPlotly({
    
    sub_df_init <-  df[str_detect(df$Description, regex("kif1c", ignore_case = T)), ]
    if(!is.null(input$search) && input$search !="") {
      sub_df1 <- df[str_detect(df$Description, regex(input$search, ignore_case = T)), ]
      sub_df2 <- df[str_detect(df$Genes, regex(input$search, ignore_case = T)), ]
      sub_df <- rbind(sub_df1, sub_df2)
      sub_df <- sub_df[!duplicated(sub_df),]
    } else {
      sub_df <- sub_df_init
    }
    
    # x axis threshold: 3-fold enrichment
    xThre = 3
    # y axis threshold: 0.05 p-value
    yThre = 1/0.05

    p <- ggplot() +
      geom_vline(xintercept=xThre, linetype="dashed", color = "red")+
      geom_hline(yintercept=yThre, linetype="dashed", color = "red")+
      annotate('rect', xmin=2^(-6), xmax=xThre, ymin=1, ymax=10^9, alpha=.1, fill='black')+
      annotate('rect', xmin=2^(-6), xmax=2^12, ymin=1, ymax=yThre, alpha=.1, fill='black')+
      geom_point(data=df, aes(x = Enrichment, y = pValue^(-1), text = Genes)) + 
      geom_point(data=sub_df, aes(x=Enrichment, y=pValue^(-1), text = Genes), colour="#00CED1")+
      # geom_point(data=h_df, aes(x=Enrichment, y=pValue^(-1), text = Genes), colour="orange")+
      scale_x_continuous(trans='log2', breaks=trans_breaks('log2', function(x) 2^x, 10),
                         limits = c(2^(-6), 2^12), expand = c(0, 0),
                         # labels=trans_format('log2', comma)) +
                         labels=c(-6, -4, -2, 0, 2, 4, 6, 8, 10, " ")) +
      scale_y_continuous(trans='log10', breaks=trans_breaks('log10', function(x) 10^x, 10),
                         limits = c(1, 10^9), expand = c(0, 0),
                         labels=trans_format('log10', comma)) +
      annotate(geom = "text", label= "Fold Enrichment > 3\np-value < 0.05", 
               x = 40, y = 50000000, size = 5) +
      theme_classic()+
      theme(
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(size = 18, hjust = 0.5)
      )+
      # xlab(bquote('log'[2]*' enrichment'))+
      # ylab(bquote('-log'[10]*' p-value'))+
      xlab('log2 enrichment')+
      ylab('-log10 p-value')+
      ggtitle("Protein enrichment")
    
    # process the plot to allow interaction
    ggplotly(p, tooltip = c("text")) %>% layout(hovermode = "closest")
  })
  
  output$infoTable2 <- renderDT({
    
    search_df_init <-  data[str_detect(data$Description, regex("kif1c", ignore_case = T)), ]
    
    if(!is.null(input$search) && input$search !=""){
      search_df1 <- data[str_detect(data$Description, regex(input$search, ignore_case = T)), ]
      search_df2 <- data[str_detect(data$Genes, regex(input$search, ignore_case = T)), ]
      search_df <- rbind(search_df1, search_df2)
      search_df <- search_df[!duplicated(search_df),]
    } else {
      search_df <- search_df_init
    }
    
    if (nrow(search_df) != 0) {
      geneName <- str_extract(search_df$Description, "(?<=GN=)[^ ]+")
      search_df$Genes <- paste0("<a href='",
                                "https://www.genecards.org/cgi-bin/carddisp.pl?gene=", geneName,
                                "' target='_blank'>",
                                search_df$Genes,
                                "</a>")
    }
    
    datatable(search_df[c("Genes", "Description", "Enrichment", "p.value")], 
              caption = "Search results", rownames = FALSE,
              options = list(columnDefs = list(list(width = '100px', targets = 0), 
                                               list(width = '400px', targets = 1),
                                               list(width = '100px', targets = 2),
                                               list(width = '100px', targets = 3)
              ), 
              dom = 't'), escape = FALSE)
  })
  
  # output$infoTable1 <- renderDT({
  #   
  #   if (is.null(input$hover)) {
  #     hover_df <-  data[data$Genes == "KIF1C", ]
  #   } else {
  #     hoverIndex <- req(input$hover)
  #     hover_df <- data.frame(data[hoverIndex+1, ])
  #   }
  #   
  #   datatable(hover_df, caption = "Hover results", rownames = FALSE,
  #             options = list(columnDefs = list(list(width = '100px', targets = 0), 
  #                                              list(width = '400px', targets = 1),
  #                                              list(width = '100px', targets = 2),
  #                                              list(width = '100px', targets = 3)), 
  #                            dom = 't'))
  # })
}

# Run the application 
shinyApp(ui = ui, server = server)


