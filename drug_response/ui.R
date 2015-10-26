library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    
    # Application title
    titlePanel("Drug response curve"),
    
    # Sidebar with a slider input for the number of bins
    sidebarLayout(
        sidebarPanel(
        #    sliderInput("bins",
        #                "Number of bins:",
        #                min = 1,
        #                max = 50,
        #                value = 30)
        
        fileInput("file", label = "File" ),
        textInput("xlab", label = "X axis label", value = "Log drug concentration"),
        textInput("ylab", label = "Y axis label", value = "viability"),
        checkboxInput("sf50", label = "Show SF50", value = TRUE)
        ),
        # Show a plot of the generated distribution
        mainPanel(
            fluidRow(
                plotOutput("drugresponsePlot", height="600px", width="800px")
            ),
            fluidRow(
                dataTableOutput(outputId="table")
            )   
        )
    )
))
