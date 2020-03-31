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
      helpText("Please, indicate teh identifier at Pandora
               of the variant you want to classify [ctl + shift + j]"),
      
      numericInput(inputId = "variantId", 
                   label = "Specify variant iD in Pandora",
                   value = NA, min = 0)
   ), 
##example, to delete
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

# Define server logic required to connect to DB and run all functions to classifiy RASopathy variant
server <- function(input, output) {
   
  ## conection to UCSC
      my_connection <- dbConnect(
        MySQL(),
        user="genome",
        dbname="hg38",
        host="genome-euro-mysql.soe.ucsc.edu",
        port=3306
      )
    
  ## connection to Pandora DB
      get_db_parameters <- function(db_conf) {
        params <- read.table("~/GitHub/RASo_variantsClassification/params.csv", sep = ",", stringsAsFactors = FALSE)
        return(list(user = params$V2[1],
                    password = params$V2[2],
                    dbname = params$V2[3],
                    host = params$V2[4],
                    port = params$V2[5]))
      }
      
      ## connects to the NGS BD using a config file
      ## param db_conf: a full path to a config file (CSV)
      
      
      db_connect_postgres <- function(db_conf) {
        drv <- dbDriver("PostgreSQL")
        
        db_conf <- get_db_parameters(db_conf)
        
        con <- dbConnect(drv,
                         user = db_conf[['user']],
                         password = db_conf[['password']],
                         dbname = db_conf[['dbname']],
                         host = db_conf[['host']],
                         port = db_conf[['port']])
        
        return(con)
      }
      
      con <- db_connect_postgres('db_pandora.conf')
      
      ####### Stablished criteria are different depending on GOF genes and LOF genes. 
      
      #LOF Genes: NF1 & SPRED1. BP1 % PVS1 criteria are only applicable to LOF genes
      #GOF Genes: BRAF, CBL, HRAS, KRAS, LZTR1, NF1, NRAS, MAP2K1, MAP2K2, PTPN11, RAF1, RIT1, SHOC2, SOS1, SOS2, SPRED1 i RASA1,
      
      ##### Calling RASopathies gene information
      domain_groupRAF <- read.csv("~/GitHub/RASo_variantsClassification/domini_grupRAF.csv")
      domain_groupRAS <- read.csv("~/GitHub/RASo_variantsClassification/domini_grupRAS.csv")
      domain_groupSOS <- read.csv("~/GitHub/RASo_variantsClassification/domini_grupSOS.csv")
      domain_groupMAPK <- read.csv("~/GitHub/RASo_variantsClassification/domini_grupMAPK.csv")
      Transcripts_RASos <- read.csv("~/GitHub/RASo_variantsClassification/Transcripts_RASos.csv")
      
      
      ###outputs examples
      
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

