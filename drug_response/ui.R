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
        fileInput("file", label = "File",width="200px" ),
        textInput("xlab", label = "X axis label", value = "Log drug concentration", width="150px"),
        textInput("ylab", label = "Y axis label", value = "viability", width="150px"),
        selectInput("xmax", label = "Maximum X value", 
                    choices = list("1e1" = 1, "1e2" = 2,"1e3" = 3, "1e4" = 4, 
                                   "1e5" = 5, "1e6" = 6, "1e7" = 7), 
                    selected = 5, width="150px"),
        selectInput("ymax", label = "Maximum Y value", 
                    choices = list("1" = 1, "1.1" = 1.1,"1.2" = 1.2), 
                    selected = 1, width="150px"),
        checkboxInput("sf50", label = "Show SF50", value = TRUE),
        textInput("legend_x", label = "Legend X location", value = "10", width="150px"),
        textInput("legend_y", label = "Legend Y location", value = "0.2", width="150px"),
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
