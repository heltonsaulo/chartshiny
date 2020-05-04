#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(chartslogsym)

# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Log-symmetric control charts"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      h3("Data Input"),
      h5("Upload a CSV file with in-control data in the first column and out-of-control data in the second column"),
      # Input: Select a file ----
      fileInput("file1", "Choose CSV File - Phases I and II",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      # Horizontal line ----
      #tags$hr(),
      
      # Input: Checkbox if file has header ----
      checkboxInput("header", "Header", TRUE),
      
      # Input: Select separator ----
      radioButtons("sep", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t"),
                   selected = ","),
      
      # Input: Select quotes ----
      # radioButtons("quote", "Quote",
      #              choices = c(None = "",
      #                         "Double Quote" = '"',
      #                         "Single Quote" = "'"),
      #             selected = '"'),
      
      # Horizontal line ----
      tags$hr(),
      
      # Input: Select number of rows to display ----
      radioButtons("disp", "Display",
                   choices = c(Head = "head",
                               All = "all"),
                   selected = "head"),
      

      
      # Horizontal line ----
      tags$hr(),
      
      h3("Control Parameters - Phase I"),
      
      numericInput('n', 'Number of Samples (n):', 10, min = 1, max = 100), 
      numericInput('m', 'Sample size of each sample (m):', 5, min = 1, max = 15),
      numericInput('gamma', 'False Alarm Rare (FAR) - gamma:', 0.0027, min = 0.001, max = 0.999), 
      numericInput('p', 'Percentile of Interest:', 0.01, min = 0.01, max = 0.99),
      numericInput('Boot', 'Bootstrap replications:', 10000, min = 100, max = 20000),
      
      # Horizontal line ----
      tags$hr(),
      
      h3("Control Parameters - Phase II"),
      numericInput('n2', 'Number of Samples (n):', 10, min = 1, max = 100), 
      numericInput('m2', 'Sample size of each sample (m):', 5, min = 1, max = 15),
      numericInput('linf', 'Lower grid limit for the extra parameter xi:', 1, min = 0.01, max = 100),
      numericInput('lsup', 'Upper grid limit for the extra parameter xi:', 20, min = 1, max = 100),
      numericInput('by', 'grid: increment of the sequence:', 0.1, min = 0.01, max = 1),
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Data file ----
      h3("Data sets"),
      tableOutput("contents"),
      h3("Control Chart - Phase I"),
      h5("Estimation results"),
      tableOutput("resultsA"),
      plotOutput("plot1"),
      h3("Control Chart - Phase II"),
      plotOutput("plot2")
      
    )
    
  )
  
  
  
  
)

# Define server logic to read selected file ----
server <- function(input, output) {
  
  
  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, head of that data file by default,
  # or all rows if selected, will be shown.
  
  data1 <- reactive({
  req(input$file1)
  
  
  
  df <- read.csv(input$file1$datapath,
                 header = input$header,
                 sep = input$sep)
  # quote = input$quote)
  
  if(input$disp == "head") {
    return(head(df))
  }
  else {
    return(df)
  }
  
  })
  
  data2 <- reactive({
    req(input$file1)
    
    
    
    df <- read.csv(input$file1$datapath,
                   header = input$header,
                   sep = input$sep)
    # quote = input$quote)
    
  
      return(df)
    
    
  })
  
  output$contents <- renderTable({
    
    data1()
    
  },digits=4)
  
  
  
  results1 <- reactive({
  
    n <- input$n
    m <- input$m
    gamma <- input$gamma
    p     <- input$p
    boot <- input$Boot
    x <- data2()[,1]
    y <- data2()[,2]
    #print("Control Chart - Phase I")
    resI <- chart.logsym1(data=x, p=p, n=n, m=m, gamma=gamma, boot=boot)  
    return(resI)
  })
    
  
  
  output$resultsA <- renderTable({ 
    x <- data2()[,1]
    
  resmle <- best.logsym(y = x)
    
  mylist <-   list(family = results1()$family,
                   eta    = results1()$eta,
                   phi    = results1()$phi,
                   xi     = results1()$xi,
                   LCL    = results1()$lcl,
                   CL     = results1()$cl,
                  UCL     = results1()$ucl,
                  AIC     = resmle$AIC,
                  BIC     = resmle$BIC)
    
  return(mylist)
  
  })
    
  
  output$plot1 <- renderPlot({
    
    n <- input$n
    m <- input$m
    gamma <- input$gamma
    p     <- input$p
    boot <- input$Boot
    x <- data2()[,1]
    y <- data2()[,2]
    
    familyhat = results1()$family
    xihat     = results1()$xi
    lcl    = results1()$lcl
    cl     = results1()$cl
    ucl    = results1()$ucl
    
    datam <- matrix (x,n,m,byrow = TRUE )
    Wt <- double()
    for(i in 1:nrow(datam)){
      parp <- mle.logsym2(datam[i,], family = familyhat, xi = xihat)
      Wt[i] <- quantile.logsym(p = p, parp, family = familyhat, xi = xihat)
    }
    
    bchart1(x=Wt, lcl=lcl, ucl=ucl, cl=cl, n=n)
    
  })
  
  

  output$plot2 <- renderPlot({
    
    n <- input$n2
    m <- input$m2
    gamma <- input$gamma
    p     <- input$p
    boot <- input$Boot
    
    linf <- input$linf
    lsup <- input$lsup
    by  <- input$by
        
    x <- data2()[,1]
    y <- data2()[,2]
    familyhat = results1()$family
    xihat     = results1()$xi
    lcl    = results1()$lcl
    cl     = results1()$cl
    ucl    = results1()$ucl
    
    chart.logsym2(data=y, p=p, n=n, m=m, family=familyhat,
                  lcl=lcl, ucl=ucl , cl=cl,
                  linf=1, lsup=20,by=.1)
  })
  
}
# Run the app ----
shinyApp(ui, server)
