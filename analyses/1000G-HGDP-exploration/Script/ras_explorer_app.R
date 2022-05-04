#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
source("Functions/load_data.R")
# ras_data <- read_ras_tables("../Data")
ras_data <- readr::read_rds("../Data/combined_data.rds")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("RAS Explorer"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          checkboxInput("TF", "1240K"),
          checkboxInput("TVonly", "transversions only", TRUE),
          checkboxInput("mapMasked", "mapping mask filter", TRUE),
          radioButtons("rightPop", "Right Population",
                       c("CEU2", "FIN2", "GBR2", "IBS2", "TSI2",
                       "Basque2", "BergamoItalian2", "French2", "Orcadian2", "Russian2", "Sardinian2", "Tuscan2")),
          radioButtons("rasAF", "Allele Frequency Cutoff", 
                       c("01", "02", "05", "10", "20", "All", "Common")),
          radioButtons("dataset", "Dataset", c("1000G", "HGDP"))
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$distPlot <- renderPlot({
      ras_data %>% dplyr::filter(Right == input$rightPop & dataset == input$dataset &
                                 !(Group %in% c("YRI", "PEL", "CHB", "Han", "Karitiana", "Mbuti")) &
                                 TF == input$TF & rasAF == input$rasAF &
                                 tvOnly == input$TVonly &
                                 mapMasked == input$mapMasked) %>%
        ggplot(aes(x = Left, y = RAS, col = Group)) + geom_point() +
          geom_errorbar(aes(ymin = RAS - StdErr, ymax = RAS + StdErr)) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
