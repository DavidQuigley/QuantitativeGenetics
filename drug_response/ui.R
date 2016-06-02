library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    
    singleton(tags$head(HTML(
        '
        <script type="text/javascript">
        $(document).ready(function() {
        $("#download_PDF").attr("disabled", "true").attr("onclick", "return false;");
        
        Shiny.addCustomMessageHandler("download_PDF", function(message) {
        $("#download_PDF").removeAttr("disabled").removeAttr("onclick").html(
        "<i class=\\"fa fa-download\\"></i>Download PDF");
        });
        })
        </script>
        '
    ))),
    
    # Application title
    titlePanel("Drug response curve"),
    
    # Sidebar with a slider input for the number of bins
    sidebarLayout(
        sidebarPanel(
        fileInput("file", label = "File",width="300px" ),
        h4("Labels"),
        textInput("xlab", label = "X axis", value = "Log drug concentration", width="200px"),
        textInput("ylab", label = "Y axis", value = "viability", width="200px"),
        hr(),
        h4("X axis"),
        fluidRow(
            column( width=4, selectInput("min_x", label = "From", 
                    choices = list("0"=0, "10" = 10, "1e2" = 100,"1e3" = 1e3, 
                                   "1e4" = 1e4, "1e5" = 1e5, "1e6" = 1e6, "1e7" = 1e7), 
                    selected = 0, width="75px") ),
            column( width=4, selectInput("max_x", label = "To", 
                    choices = list("1e1" = 1e1, "1e2" = 1e2,"1e3" = 1e3, "1e4" = 1e4, 
                                   "1e5" = 1e5, "1e6" = 1e6, "1e7" = 1e7), 
                    selected = 1e5, width="75px") )
        ),
        
        fluidRow(
            column( width=4, selectInput("show_x_logtic", label = "Inner ticks", 
                        choices = list("Yes"=TRUE, "No"=FALSE), 
                        selected = TRUE, width="150px") ),
            column( width=4, selectInput("show_x_exponent", label = "Values as", 
                        choices = list("10^N" = TRUE, "N" = FALSE ),
                        selected = TRUE, width="150px") )   
        ),
        fluidRow( 
            column( width=4, 
                    selectInput("axis_pointsize", label = "Axis Point size", 
                                choices = list("1" = 1, "1.5" = 1.5,"2" = 2), 
                                selected = 1, width="100px") ),
            column( width=4, 
                    selectInput("axis_labelsize", label = "Label size", 
                                choices = list("1" = 1, "1.5" = 1.5, "2" = 2, "2.5"=2.5), 
                                selected = 1.5, width="100px") )
        ),
        hr(),
        h4("Y axis"),
        fluidRow( 
            column( width=4, strong("From:"), br(), strong("0") ),
            column( width=4, selectInput("max_y", label = "To", 
                    choices = list("1" = 1, "1.05"=1.05, "1.1" = 1.1, "1.2" = 1.2, "1.3"=1.3, "1.4"=1.4, "1.5"=1.5), 
                    selected = 1, width="150px"))
        ),
        selectInput("barmultiple", label = "Error Bar SEM", 
                    choices = list("1" = 1, "2" = 2,"3" = 3), 
                    selected = 1, width="75px"),
        checkboxInput("sf50", label = "Show EC50", value = TRUE),
        hr(),
        h4("Legend location"),
        fluidRow(
            column( width=4,
                textInput("legend_x", label = "X", value = "10", width="150px")),
            column( width=4,
                textInput("legend_y", label = "Y", value = "0.2", width="150px") )
        ),
        downloadButton("download_PDF")
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
