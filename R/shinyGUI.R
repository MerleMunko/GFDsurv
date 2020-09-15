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
                      titlePanel("Tests for GFDsurv"),
                      sidebarLayout(
                        sidebarPanel(
                          splitLayout(
                            fileInput("infile", "Choose CSV File",
                                      accept = c(
                                        "text/csv",
                                        "text/comma-separated-values,text/plain",
                                        ".csv")),
                            checkboxInput("header", "Header", TRUE),
                            selectInput("sep","Seperator in csv", c(",",
                                                                    ";",
                                                                    ".",
                                                                    "|"))

                          ),

                          tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;
                              }
                              "))), #for selectinput in splitlayout with full dropdown view


                          selectInput("Method", "Select Testing Method:",
                                      c("CASANOVA: Cumulative Aalen survival analyis-of-variance" = "casanova",
                                        "medSANOVA: Median survival analyis-of-variance"= "medSANOVA",
                                        "copSANOVA: concordance probability SANOVA"="copSANOVA")),

                          splitLayout(

                            uiOutput(outputId = 'dynamicInput'),


                            textInput("formula", "Formula ", "timeFactor ~ FactorA*FactorB")

                          ),

                          splitLayout(cellWidths = c("20%","60%","20%"),

                            checkboxGroupInput("Weights", "Choose weights:",selected = c("crossing","proportional"),
                                      choiceNames = list("Crossing", "Proportional"),
                                      choiceValues = list("crossing", "proportional")),


                            selectInput("weights1","User specifeid directions of form w(x) = x^r(1-x)^g",
                                        paste0("r = ",expand.grid(1:10,1:10)[,1],
                                               ", g = ",expand.grid(1:10,1:10)[,2]),
                                        multiple=TRUE,selectize = TRUE)


                          ),


                          splitLayout(cellWidths = c("30%","70%"),

                            radioButtons("variante", "Variance estimation based on",
                                       c("one.sided"= "onesided",
                                         "two.sided" = "twosided"), inline = TRUE),

                            numericInput("var_level", "confidence intervalle with level",
                                                   min = 0, max = 1,
                                                   value = 0.05, width = "20%")
                          ),



                          splitLayout(

                            numericInput("nperm", "number of permutations", value = 1999)

                          ),

                          splitLayout(
                            numericInput("nboot", "number of bootstraps iterations", value = 99),

                            selectInput("bootstrapweights", "Select wildbootstrap weights:",
                                        c("pois" = "pois",
                                          "norm"= "norm",
                                          "weird"="weird",
                                          "corrLibPois"= "corrLibPois",
                                          "corrLibNorm" = "corrLibNorm",
                                          "corrLibWeird" = "corrLibWeird"
                                        ))
                          ),

                          splitLayout(
                            numericInput("tau", "Choose Tau",NULL),
                            actionButton("tau_suggest", "Calculate automaticly tau", class = "btn-primary")
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

             shinyjs::hide("Weights")
             shinyjs::hide("weights1")

           }


           if (input$Method == "casanova") {
             # data  <- as.data.frame(datasetInput())

             shinyjs::show("Weights")
             shinyjs::show("weights1")



           }


           if (input$Method != "medSANOVA") {
             # data  <- as.data.frame(datasetInput())
             shinyjs::hide("variante")
             shinyjs::hide("var_level")

           }
           if (input$Method == "medSANOVA") {
             # data  <- as.data.frame(datasetInput())
             shinyjs::show("variante")
             shinyjs::show("var_level")


           }
           if (input$Method != "copSANOVA") {
             # data  <- as.data.frame(datasetInput())
             shinyjs::hide("bootstrapweights")
             shinyjs::hide("nboot")
             shinyjs::hide("tau")
             shinyjs::hide("tau_suggest")
             shinyjs::show("nperm")

           }
           if (input$Method == "copSANOVA") {
             # data  <- as.data.frame(datasetInput())
             shinyjs::show("bootstrapweights")
             shinyjs::show("nboot")
             shinyjs::show("tau")
             shinyjs::show("tau_suggest")
             shinyjs::hide("nperm")



           }




          })## observeevent


         values <- reactiveValues()

         output$dynamicInput <- renderUI({

           # This input exists if the `static`
           # one is equal to `A` only
           if (input$Method == "casanova" || input$Method == "medSANOVA"|| input$Method == "copSANOVA") {
             selectInput(inputId = 'dynamic',
                         label = "Name of event",
                         choices = colnames(datasetInput()))
           } else {
             return(NULL)
           }

         })

         ## this bit fixes the issue
         observe({
           if (input$Method == "casanova" || input$Method == "medSANOVA"|| input$Method == "copSANOVA") {
             values$dyn <- input$dynamic
           } else {
             values$dyn <- NULL
           }
         })


         observeEvent(input$tau_suggest, {
           if (input$formula == "~ + ") {

             output$result <- renderPrint({
               "'formula' missing or invalid"
             })
           }
           updateNumericInput(session, "tau", value = copsanova_tau(formula= isolate(input$formula),
                                                                    event = input$dynamic,
                                                                    data = isolate(data)))

         }
         )



         observeEvent(input$process, {


          if (input$formula == "~ + ") {

            output$result <- renderPrint({
              "'formula' missing or invalid"
            })

          } else {
            if (input$Method == "casanova" ){
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

            if (input$Method == "medSANOVA" ){
              data <- as.data.frame(datasetInput())
              output$result <- renderPrint({
                medSANOVA(formula= isolate(input$formula),
                         event = input$dynamic,
                         data = isolate(data),
                         variant = isolate(input$variante),
                         nperm = isolate(input$nperm),
                         alpha = isolate(input$alpha),
                         nested.levels.unique = FALSE
                         )
              })
            }
            if (input$Method == "copSANOVA" ){
              data <- as.data.frame(datasetInput())
              output$result <- renderPrint({
                copsanova(formula= isolate(input$formula),
                          event = input$dynamic,
                          data = isolate(data),
                          BSiter = isolate(input$nboot),
                          weights = isolate(input$bootstrapweights),
                          tau = isolate(input$tau),
                          nested.levels.unique = FALSE
                )
              })


            }

          }

         }
         ) #end of observeEvent(input$process

      }


    shinyApp(ui = ui, server = server)

}
