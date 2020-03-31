## This is the code to run script_RASo_class in a shiny enviornment

######## INSTALLATION ALL PACKAGES AND CONNECTIONS
library(shiny)
require("RPostgreSQL")
library(RMySQL)
require(stringr)



# Define UI for application that draws a histogram
ui <- fluidPage(
   
  
   # Application title
   titlePanel("RASopathy-related variant classification"),
   
   # Sidebar with a slider input for number of bins 
   # Numeric Input with variant identifier in Pandora 
  sidebarLayout(
    sidebarPanel(     
      helpText("Classify variants in RASopathy-related genes 
               automatically following ACMG guidelines"),
      
      numericInput(inputId = "variantId", 
                   label = "Specify variant iD in Pandora",
                   value = NA, min = 0)
   ), 
   
   sidebarLayout(
      sidebarPanel(
      
         sliderInput("bins",
                     "Number of bins:",
                     min = 1,
                     max = 50,
                     value = 30)
      ),
 
  # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot")
      )
   )
))

# Define server logic required to draw a histogram
server <- function(input, output) {
   
      output$tbl <- renderTable({
      conn <- dbConnect(
        drv = RMySQL::MySQL(),
        dbname = "shinydemo",
        host = "shiny-demo.csa7qlmguqrf.us-east-1.rds.amazonaws.com",
        username = "guest",
        password = "guest")
      on.exit(dbDisconnect(conn), add = TRUE)
      dbGetQuery(conn, paste0(
        "SELECT * FROM City LIMIT ", input$nrows, ";"))
    })
  
     output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      x    <- faithful[, 2] 
      bins <- seq(min(x), max(x), length.out = input$bins + 1)
      
      # draw the histogram with the specified number of bins
      hist(x, breaks = bins, col = 'darkgray', border = 'white')
 })
}

# Run the application 
shinyApp(ui = ui, server = server)

