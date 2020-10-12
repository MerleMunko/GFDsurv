#' A graphical user interface for the package GFDsurv
#'
#'
#' @aliases GFDsurvGUI
#'
#'
#' @import shiny
#'
#' @export

GFDsurvGUI <- function() {
  requireNamespace("shiny", quietly = TRUE)

  if (!("package:shiny" %in% search())) {
    attachNamespace("shiny")
  }
  requireNamespace("shinyWidgets", quietly = TRUE)
  if (!("package:shinyWidgets" %in% search())) {
    attachNamespace("shinyWidgets")
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
                          tags$style(HTML("
                                 input[type=number] {
                                                              -moz-appearance:textfield;
                                                    }
                                  input[type=number]::{
                                                  -moz-appearance:textfield;
                                                    }
                        input[type=number]::-webkit-outer-spin-button,
                        input[type=number]::-webkit-inner-spin-button {
                        -webkit-appearance: none;
                        margin: 0;
                        }
                        ")),

                        h3(id="titleLoadData","Load dataset first!", style = "color:red"),

                        shinyjs::hidden(
                          selectInput("Method", "Select Testing Method:",
                                      c("CASANOVA: Cumulative Aalen survival analyis-of-variance" = "casanova",
                                        "MedSANOVA: Median survival analyis-of-variance"= "medSANOVA",
                                        "CopSANOVA: Concordance probability survival analyis-of-variance"="copSANOVA"))
                        ),


                          splitLayout(
                            uiOutput(outputId = 'dynamicInput'),
                            shinyjs::hidden(
                             textInput("formula", "Formula ", "timeFactor ~ FactorA*FactorB")
                            )
                          ),

                        shinyjs::hidden(
                          h5(id="titleWeights",strong("Which weight functions w should be combined?"), style = "color:grey")
                        ),


                          splitLayout(cellWidths = c("35%","60%","5%"),
                                      shinyjs::hidden(
                            checkboxGroupInput("Weights", "Pre-specified weights",selected = c("crossing","proportional"),
                                      choiceNames = list("Crossing", "Proportional"),
                                      choiceValues = list("crossing", "proportional"))
                                      ),

                            shinyjs::hidden(
                            selectInput("weights1","Specify the exponents (r,g) of weights  w(x) = x^r(1-x)^g ",
                                        paste0("(",expand.grid(0:10,0:10)[,1],",",expand.grid(0:10,0:10)[,2],")"),
                                        multiple=TRUE,selectize = TRUE)
                            )

                          ),


                          splitLayout(cellWidths = c("30%","70%"),
                                      shinyjs::hidden(
                            radioButtons("variante", "Variance estimation based on",
                                       c("one.sided"= "onesided",
                                         "two.sided" = "twosided"), inline = TRUE)
                                      ),
                                  shinyjs::hidden(
                            numericInput("var_level", "confidence intervalle with level",
                                                   min = 0, max = 1,
                                                   value = 0.05, width = "20%")
                                  )
                          ),



                          splitLayout(
                            shinyjs::hidden(
                            numericInput("nperm", "Number of permutations", value = 1999)
                            )
                          ),

                          splitLayout(cellWidths = c("20%","10%","20%","10%","30%"),
                                      shinyjs::hidden(
                          shinyWidgets::sliderTextInput(
                              inputId = "sliderBoot",
                              label = "Bootstrap type:",
                              choices = c("wild","weird"),
                              selected = "wild")
                                      ),
                            shinyjs::hidden(
                              selectInput("Platz1", "Distribution:",
                                          c("Poisson"="pois","Normal"="norm"))
                            ),
                            shinyjs::hidden(
                            selectInput("methodBoot", "Distribution:",
                                        c("Poisson"="pois","Normal"="norm"))
                            ),
                            shinyjs::hidden(
                              selectInput("Platz2", "Distribution:",
                                          c("Poisson"="pois","Normal"="norm"))
                            ),
                            shinyjs::hidden(
                            checkboxInput("correction", "correction for liberal test", TRUE)
                            )

                          ),

                          splitLayout(cellWidths = c("15%","85%"),
                                      shinyjs::hidden(
                          numericInput("nboot", "Number of bootstrap iterations", value = 99)
                                      )
                          ),

                        splitLayout(cellWidths = c("15%","85%"),
                            shinyjs::hidden(
                            numericInput("tau", "Endpoint tau of the relevant time window [0,tau]",NULL)
                            )),

                        shinyjs::hidden(
                          actionButton("tau_suggest", "Specify tau automatically", class = "btn-primary")
                        ),

                        shinyjs::hidden(
                          actionButton("process", "Calculate", class = "btn-primary")
                        )

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


         observeEvent(input$infile, {

           if(is.null(input$infile)){
             shinyjs::hide(id = "Method")
             shinyjs::hide(id = "weights1")
             shinyjs::hide(id = "Weights")
             shinyjs::hide(id = "formula")
             shinyjs::hide(id = "variante")
             shinyjs::hide(id = "var_level")
             shinyjs::hide(id = "nperm")
             shinyjs::hide(id = "sliderBoot")
             shinyjs::hide(id = "methodBoot")
             shinyjs::hide(id = "correction")
             shinyjs::hide(id = "nboot")
             shinyjs::hide(id = "tau")
             shinyjs::hide(id = "tau_suggest")
             shinyjs::hide(id = "process")
             shinyjs::hide(id = "titleWeights")
             shinyjs::hide(id = "Platz1")
             shinyjs::hide(id = "Platz2")

           }else{
             shinyjs::show(id = "Method")
             shinyjs::show(id = "formula")
             shinyjs::show(id = "process")
             shinyjs::hide(id = "titleLoadData")


             observeEvent(input$Method, {

               if (input$Method != "casanova") {

                 shinyjs::hide(id = "titleWeights")
                 shinyjs::hide("Weights")
                 shinyjs::hide("weights1")

               }


               if (input$Method == "casanova") {

                 shinyjs::show(id = "titleWeights")
                 shinyjs::show("Weights")
                 shinyjs::show("weights1")

               }


               if (input$Method != "medSANOVA") {

                 shinyjs::hide("variante")
                 shinyjs::hide("var_level")

               }

               if (input$Method == "medSANOVA") {

                 shinyjs::show("variante")
                 shinyjs::show("var_level")

               }

               if (input$Method != "copSANOVA") {
                 # data  <- as.data.frame(datasetInput())
                 shinyjs::hide("methodBoot")
                 shinyjs::hide("nboot")
                 shinyjs::hide("sliderBoot")
                 shinyjs::hide("correction")
                 shinyjs::hide("tau")
                 shinyjs::hide("tau_suggest")
                 shinyjs::show("nperm")

               }
               if (input$Method == "copSANOVA") {
                 # data  <- as.data.frame(datasetInput())
                 shinyjs::show("methodBoot")
                 shinyjs::show("nboot")
                 shinyjs::show("sliderBoot")
                 shinyjs::show("correction")
                 shinyjs::show("tau")
                 shinyjs::show("tau_suggest")
                 shinyjs::hide("nperm")

               }




             })## observeevent




             observeEvent(input$sliderBoot, {

               if (input$sliderBoot != "weird" && input$Method == "copSANOVA") {
                 # data  <- as.data.frame(datasetInput())
                 shinyjs::show("methodBoot")

               }
               if (input$sliderBoot == "weird") {
                 # data  <- as.data.frame(datasetInput())
                 shinyjs::hide("methodBoot")

               }


             })

           }
         })












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
           data <- as.data.frame(datasetInput())
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

            rg <- list()
            givenWeights <- isolate(input$Weights)
            if("proportional" %in% givenWeights){
              rg[[length(rg)+1]] <- c(0,0)
            }
            if("crossing" %in% givenWeights){
              crossing <- TRUE
            } else {
              crossing <- FALSE
            }

            inputWeights <- isolate(input$weights1)

            if(is.null(inputWeights)){}else{
              expand.grid(0:10,0:10)
              kombiAll <- paste0("(",expand.grid(0:10,0:10)[,1],",",expand.grid(0:10,0:10)[,2],")")
              kombi <- expand.grid(0:10,0:10)[which(kombiAll %in% inputWeights),]

              for (i in 1:(dim(kombi)[1])){
                rg[[length(rg)+1]] <- as.numeric(kombi[i,])
              }
            }

            output$result <- renderPrint({
              casanova(formula= isolate(input$formula),
                       event = input$dynamic,
                       data = isolate(data),
                       nperm = isolate(input$nperm),
                       cross = crossing,
                       nested.levels.unique = FALSE,
                       rg = isolate(rg))
            })
            }

            if (input$Method == "medSANOVA" ){
              data <- as.data.frame(datasetInput())
              output$result <- renderPrint({
                medsanova(formula= isolate(input$formula),
                         event = input$dynamic,
                         data = isolate(data),
                         var_method = isolate(input$variante),
                         nperm = isolate(input$nperm),
                         nested.levels.unique = FALSE
                         )
              })
            }
            if (input$Method == "copSANOVA" ){
              data <- as.data.frame(datasetInput())
              bootstrapMethod <- isolate(input$sliderBoot)
              bootstrapDis <- isolate(input$methodBoot)
              correction <- isolate(input$correction)
             if(bootstrapMethod == "wild"){
               if(bootstrapDis =="pois"){
                 if(correction == TRUE){
                   weights <- "corrLibPois"
                 } else {
                   weights <- "pois"
                 }
               } else {
                 if(correction == TRUE){
                  weights <- "corrLibNorm"
                }else{
                  weights <- "norm"
               }
               }
             } else {
               if(correction == TRUE){
                 weights <- "corrLibWeird"
               }else{
                 weights <- "weird"
               }

             }

              output$result <- renderPrint({
                copsanova(formula= isolate(input$formula),
                          event = input$dynamic,
                          data = isolate(data),
                          BSiter = isolate(input$nboot),
                          weights = weights,
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

