#' A graphical user interface for the package GFDsurv
#'
#'
#' @aliases GFDsurvGUI
#'
#' @export

GFDsurvGUI <- function() {
  requireNamespace("shiny", quietly = TRUE)
  if (!("package:shiny" %in% search())) {
    attachNamespace("shiny")
  }



      ui <- fluidPage(theme = shinythemes::shinytheme("cerulean"),
                      shinyjs::useShinyjs(),
                      titlePanel("Tests for GFDsurv?"),
                      sidebarLayout(
                        sidebarPanel(
                          splitLayout(
                            fileInput("infile", "Choose CSV File",
                                      accept = c(
                                        "text/csv",
                                        "text/comma-separated-values,text/plain",
                                        ".csv")),
                            checkboxInput("header", "Header", TRUE),
                            textInput("sep","Seperater")

                          ),

                          selectInput("Method", "Select Testing Method:",
                                      c("CASANOVA: Cumulative Aalen survival analyis-of-variance" = "casanova",
                                        "test-2 " = "test2",
                                        "test-3" = "test3")),

                          splitLayout(

                            uiOutput(outputId = 'dynamicInput'),


                            textInput("formula", "Formula ", "timeFactor ~ FactorA*FactorB")

                          ),


                            mainPanel(strong("Select Weights")),


                          splitLayout(cellWidths = c("20%","20%","60%"),
                            checkboxInput("crossing", "Crossing", TRUE),
                            checkboxInput("proportional", "Proportional", TRUE),
                            selectInput("weights1","User specifeid directions of form w(x) = x^r(1-x)^g",
                                        paste0("r = ",expand.grid(1:10,1:10)[,1],
                                               ", g = ",expand.grid(1:10,1:10)[,2]),
                                        multiple=TRUE,selectize = TRUE)

                          ),


                          splitLayout(
                            numericInput("nperm", "nperm/nbsp", value = 1999),

                            numericInput("alpha", "Alpha", value = 0.05, min = 0, max = 1)

                          ),

                          actionButton("process", "Calculate", class = "btn-primary")
                        ),


                        mainPanel(

                          verbatimTextOutput("result")


                        )
                      )
      )


       server <- function(input, output,session) {

         datasetInput <- reactive({

           req(input$infile)

           if (is.null(input$infile))
             return(NULL)
           read.csv(input$infile$datapath, header = input$header, sep = as.character(input$sep))
         })

         observeEvent(input$Method, {

           if (input$Method != "casanova") {
             # data  <- as.data.frame(datasetInput())

             shinyjs::hide("proportional")
             shinyjs::hide("crossing")

           }


           if (input$Method == "casanova") {
             # data  <- as.data.frame(datasetInput())

             shinyjs::show("proportional")
             shinyjs::show("crossing")


           }

          })## observeevent


         values <- reactiveValues()

         output$dynamicInput <- renderUI({

           # This input exists if the `static`
           # one is equal to `A` only
           if (input$Method == "casanova") {
             selectInput(inputId = 'dynamic',
                         label = "event",
                         choices = colnames(datasetInput()))
           } else {
             return(NULL)
           }

         })

         ## this bit fixes the issue
         observe({
           if (input$Method == "casanova") {
             values$dyn <- input$dynamic
           } else {
             values$dyn <- NULL
           }
         })




         observeEvent(input$process, {


          if (input$formula == "~ + ") {

            output$result <- renderPrint({
              "'formula' missing or invalid"
            })

          } else {
            data <- as.data.frame(datasetInput())
            output$result <- renderPrint({
              casanova(formula= isolate(input$formula),
                       event = input$dynamic,
                       data = isolate(data),
                       nperm = isolate(input$nperm),
                       alpha = isolate(input$alpha),
                       cross = isolate(input$crossing),
                       nested.levels.unique = FALSE,
                       rg = list(c(0,0)))
            })

          }

         }
         ) #end of observeEvent(input$process

      }


    shinyApp(ui = ui, server = server)

}
