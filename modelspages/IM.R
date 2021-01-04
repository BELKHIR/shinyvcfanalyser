    IMParams <- function(failed = FALSE) {  
           
    # params = (s,nu1,nu2,T,m12,m21)

    # Isolation-with-migration model with exponential pop growth.

    # s: Size of pop 1 after split. (Pop 2 has size 1-s.)
    # nu1: Final size of pop 1.
    # nu2: Final size of pop 2.
    # T: Time in the past of split (in units of 2*Na generations) 
    # m12: Migration from pop 2 to pop 1 (2*Na*m12)
    # m21: Migration from pop 1 to pop 2

    
      modalDialog(      
        fluidRow( width =largeur*2,style = "height:100px;",
        column(width = largeur, sliderInput("popSize", label = "% Size of pop 1 after split bounds", min = 1, max = 99, value = c(1, 99))),
        column(width = largeur, numericInput("popSize0", "Start value", value = 50))
        ),
        fluidRow(width =largeur*2,style = "height:100px;",
        column(width = largeur,sliderInput("nu1", label = "Final size of pop 1 bounds", min = 1, 
                     max = 100, value = c(1, 100))),
        column(width = largeur, numericInput("nu10", "Start value", value = 1))
        ),
        
        fluidRow( width =largeur*2, style = "height:100px;",
        column(width = largeur, sliderInput("nu2", label = "Final size of pop 2 bounds", min = 1, max = 100, value = c(1, 100))),
        column(width = largeur, numericInput("nu20", "Start value", value = 1))
        ),

        fluidRow(width =largeur*2,style = "height:100px;",
          column(width = largeur,sliderInput("Tsplit", label = "Time in the past of split bounds", min = 1, max = 10, value = c(1, 10))),
        column(width = largeur, numericInput("Tsplit0", "Start value", value = 1))
        ),

        fluidRow( width =largeur*2,style = "height:100px;",
        column(width = largeur,sliderInput("m12", label = "Migration from pop 2 to pop 1 bounds", min = 0, max = 20, value = c(0, 20))),
        column(width = largeur, numericInput("m120", "Start value", value = 1, min=0, max=20))
        ),
        fluidRow( width =largeur*2,style = "height:100px;",
        column(width = largeur,sliderInput("m21", label = "Migration from pop 1 to pop 2 bounds", min = 0, max = 20, value = c(0, 20)))  ,
        column(width = largeur, numericInput("m210", "Start value", value = 1, min=0, max=20))
        ),                   
        fluidRow(column(width = largeur,selectInput("IMOptiMethod", "Optimization Method", choices =c("optimize","optimize_log","optimize_log_lbfgsb","dual_anneal"))),
        column(width = largeur,sliderInput("IMMaxIter", label = "How long the optimizer will run", min = 1, max = 100, value = 1))
        ),
        tags$head(tags$style(".row{height:150px;}")),
        span('(These values will be used to create lower and upper_bounds  ',
             ' used in the optimization procedure)'),
        if (failed)
          div(tags$b("Invalid params.", style = "color: red;")),

        footer = tagList(
          modalButton("Cancel"),
          actionButton("IMok", "OK")
        )
      )
    }            

    observeEvent(input$IMok, {
      # Check the params are ok.
      if (!is.null(input$IMMaxIter) ) {
        Modelparams$data <- list(
            lower_bounds = c(input$popSize[1] / 100, input$nu1[1], input$nu2[1], input$Tsplit[1], input$m12[1], input$m21[1]), 
            upper_bounds = c(input$popSize[2]/100, input$nu1[2], input$nu2[2], input$Tsplit[2], input$m12[2], input$m21[2]),
            p0 = c(input$popSize0/100, input$nu10, input$nu20, input$Tsplit0, input$m120, input$m210),
            OptiMethod  = input$IMOptiMethod,
            MaxIter     = input$IMMaxIter,
            model       = "IM" 
            )
        removeModal()
      } else {
        showModal(IMParams(failed = TRUE))
      }
    })