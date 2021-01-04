    SIParams <- function(failed = FALSE) {  
    # Strict Isolation.

    # nu1: Size of population 1 after split.
    # nu2: Size of population 2 after split.
    # Ts: Duration of divergence in isolation.
   

      modalDialog(      

        fluidRow(width =largeur*2,style = "height:100px;",
        column(width = largeur,sliderInput("nu1", label = "size of pop 1 after split bounds", min = 1, 
                     max = 100, value = c(1, 100))),
        column(width = largeur, numericInput("nu10", "Start value", value = 40))
        ),
        
        fluidRow( width =largeur*2, style = "height:100px;",
        column(width = largeur, sliderInput("nu2", label = "size of pop 2 after split bounds", min = 1, max = 100, value = c(1, 100))),
        column(width = largeur, numericInput("nu20", "Start value", value = 40))
        ),

        fluidRow(width =largeur*2,style = "height:100px;",
          column(width = largeur,sliderInput("Ts", label = "Duration of divergence in isolation bounds", min = 1, max = 100, value = c(1, 100))),
        column(width = largeur, numericInput("Ts0", "Start value", value = 40))
        ),

        fluidRow(
          column(width = largeur,selectInput("SIOptiMethod", "Optimization Method", choices =c("optimize","optimize_log","optimize_log_lbfgsb","dual_anneal"))),
          column(width = largeur,sliderInput("SIMaxIter", label = "How long the optimizer will run", min = 10, max = 100, value = 10))
        ),
        tags$head(tags$style(".row{height:150px;}")),
        span('(These values will be used to create lower and upper_bounds  ',
             ' used in the optimization procedure)'),
        if (failed)
          div(tags$b("Invalid params.", style = "color: red;")),

        footer = tagList(
          modalButton("Cancel"),
          actionButton("SIok", "OK")
        )
      )
    }            

    observeEvent(input$SIok, {
      # Check the params are ok.
      if (!is.null(input$SIMaxIter) ) {
        Modelparams$data <- list(
            lower_bounds = c( input$nu1[1], input$nu2[1],  input$Ts[1]), 
            upper_bounds = c( input$nu1[2], input$nu2[2],  input$Ts[2]),
            p0 = c( input$nu10, input$nu20,  input$Ts0),
            OptiMethod  = input$SIOptiMethod,
            MaxIter     = input$SIMaxIter,
            model       = "SI" 
            )
        removeModal()
      } else {
        showModal(SIParams(failed = TRUE))
      }
    })