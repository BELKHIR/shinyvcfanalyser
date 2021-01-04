    SCParams <- function(failed = FALSE) {  
    
    # nu1: Size of population 1 after split.
    # nu2: Size of population 2 after split.
    # m12: Migration from pop 2 to pop 1.
    # m21: Migration from pop 1 to pop 2.
    # Ts: Duration of divergence in isolation.
    # Ti: Duration of divergence with migration.      

      modalDialog(      

        fluidRow(width =largeur*2,style = "height:100px;",
        column(width = largeur,sliderInput("nu1", label = "Size of pop1 after split bounds", min = 1, 
                     max = 100, value = c(1, 100))),
        column(width = largeur, numericInput("nu10", "Start value", value = 1))
        ),
        
        fluidRow( width =largeur*2, style = "height:100px;",
        column(width = largeur, sliderInput("nu2", label = "Size of pop2 after split bounds", min = 1, max = 100, value = c(1, 100))),
        column(width = largeur, numericInput("nu20", "Start value", value = 1))
        ),                 

        fluidRow(width =largeur*2,style = "height:100px;",
          column(width = largeur,sliderInput("Ts", label = "Duration of divergence in isolation bounds", min = 1, max = 100, value = c(1, 100))),
        column(width = largeur, numericInput("Ts0", "Start value", value = 40))
        ),

        fluidRow(width =largeur*2,style = "height:100px;",
          column(width = largeur,sliderInput("Ti", label = "Duration of divergence with migration bounds", min = 1, max = 100, value = c(1, 100))),
        column(width = largeur, numericInput("Ti0", "Start value", value = 40))
        ),

        fluidRow( width =largeur*2,style = "height:100px;",
        column(width = largeur,sliderInput("m12", label = "Migration from pop 2 to pop 1 bounds", min = 0, max = 20, value = c(0, 20))),
        column(width = largeur, numericInput("m120", "Start value", value = 1, min=0, max=20))
        ),
        fluidRow( width =largeur*2,style = "height:100px;",
        column(width = largeur,sliderInput("m21", label = "Migration from pop 1 to pop 2 bounds", min = 0, max = 20, value = c(0, 20)))  ,
        column(width = largeur, numericInput("m210", "Start value", value = 1, min=0, max=20))
        ),         

        fluidRow(
          column(width = largeur,selectInput("SCOptiMethod", "Optimization Method", choices =c("optimize","optimize_log","optimize_log_lbfgsb","dual_anneal"))),
          column(width = largeur,sliderInput("SCMaxIter", label = "How long the optimizer will run", min = 10, max = 100, value = 10))
        ),
        tags$head(tags$style(".row{height:150px;}")),
        span('(These values will be used to create lower and upper_bounds  ',
             ' used in the optimization procedure)'),
        if (failed)
          div(tags$b("Invalid params.", style = "color: red;")),

        footer = tagList(
          modalButton("Cancel"),
          actionButton("SCok", "OK")
        )
      )
    }            

    observeEvent(input$SCok, {
      # Check the params are ok.
      if (!is.null(input$SCMaxIter) ) {
        Modelparams$data <- list(
            lower_bounds = c( input$nu1[1], input$nu2[1],  input$m12[1], input$m21[1], input$Ts[1], input$Ti[1]), 
            upper_bounds = c( input$nu1[2], input$nu2[2],  input$m12[2], input$m21[2], input$Ts[2], input$Ti[2]),
            p0 = c( input$nu10, input$nu20,  input$m120, input$m210, input$Ts0, input$Ti0),
            OptiMethod  = input$SCOptiMethod,
            MaxIter     = input$SCMaxIter,
            model       = "SC" 
            )
        removeModal()
      } else {
        showModal(SCParams(failed = TRUE))
      }
    })