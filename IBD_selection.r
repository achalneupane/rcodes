SELECT_IBD <- function(IBD_DF){
library(plotly)
library(ggplot2)
library(shiny)


IBD_DF
## Defining a key column in IBD which will be used for event handling in event_data()
rownames(IBD_DF) <- paste(IBD_DF$IID1,IBD_DF$IID2, sep = ":")
IBD_DF$key <- rownames(IBD_DF)


### Ui code begins below
ui <- fluidPage(
  h1("Demo - Interactive data point selection with ggplotly/plotly charts"),
  h4("Subset the dataset using the data points selected from the chart. Drag and select one or multi data points"),
  br(),
  ## Plotly plot display
  plotlyOutput("plot"),
  
  ## Data point information display post click
  verbatimTextOutput("click")
)


## Server side code begins below
server <- function(input, output, session) {
  
  ## ggplotly scatter plot
  output$plot <- renderPlotly({
    
    myplot <- ggplot(IBD_DF, aes(x = Z0, 
                                 y = Z1, 
                                 key = key)) + geom_point()
    ## in above code line, use the argument key inside ggplot which will be used for event handling
    ggplotly(myplot) %>% 
      layout(dragmode = "select")
    
  })
  
  ## returns the data related to data points selected by the user
  output$click <- renderPrint({
    
    
    ## Click_data will have the keys (row identifiers) corresponding to selected data points
    click_data <- event_data("plotly_selected")
    
    ## Event mode options. There are many more to experiment
    ## plotly_click - click on one data point
    ##  Plotly_selected - multi point select
    
    if(is.null(click_data)) 
      "No data points selected on scatter plot..." 
    else 
      # filter(IBD_DF, key %in% click_data$key) %>% select(-key)
      filter(IBD_DF, key %in% click_data$key) %>% select(key)
    ## Subsetting in above step based on selected data points and removing the key column
    
  })
}
shinyApp(ui, server)
}
