 prior_onegrow_migParams <- function(failed = FALSE) {  
    # nu1F: The ancestral population size after growth. (Its initial size is defined to be 1.)
    # nu2B: The bottleneck size for pop2
    # nu2F: The final size for pop2
    # m: The scaled migration rate
    # Tp: The scaled time between ancestral population growth and the split.
    # T: The time between the split and present
      modalDialog(      
        fluidRow(width =largeur*2,style = "height:100px;",
        column(width = largeur,sliderInput("nu1F", label = "The ancestral population size after growth", min = 0, max = 100, value = c(0, 100))),
        column(width = largeur, numericInput("nu1F0", "Start value", value = 2))
        ),
        fluidRow( width =largeur*2, style = "height:100px;",
        column(width = largeur, sliderInput("nu2B", label = "The bottleneck size for pop2", min = 0, max = 100, value = c(0, 100))),
        column(width = largeur, numericInput("nu2B0", "Start value", value = 1))
        ),
        fluidRow( width =largeur*2, style = "height:100px;",
        column(width = largeur, sliderInput("nu2F", label = "The final size for pop2", min = 0, max = 100, value = c(0, 100))),
        column(width = largeur, numericInput("nu2F0", "Start value", value = 2))
        ),
        fluidRow(width =largeur*2,style = "height:100px;",
          column(width = largeur,sliderInput("m", label = "The scaled migration rate", min = 0, max = 10, value = c(0, 10))),
        column(width = largeur, numericInput("m0", "Start value", value = 1))
        ),
        fluidRow(width =largeur*2,style = "height:100px;",
          column(width = largeur,sliderInput("Tp", label = "scaled time btw anc pop growth & split", min = 0, max = 100, value = c(0, 100))),
        column(width = largeur, numericInput("Tp0", "Start value", value = 100))
        ),
        fluidRow(width =largeur*2,style = "height:100px;",
          column(width = largeur,sliderInput("T", label = "time between the split and present", min = 0, max = 100, value = c(0, 100))),
        column(width = largeur, numericInput("T0", "Start value", value = 2))
        ),
            
        fluidRow(column(width = largeur,selectInput("OptiMethod", "Optimization Method", choices =c("optimize","optimize_log","optimize_log_lbfgsb","dual_anneal"))),
        column(width = largeur,sliderInput("prior_onegrowMaxIter", label = "How long the optimizer will run", min = 10, max = 100, value = 10))
        ),
        tags$head(tags$style(".row{height:150px;}")),
        span('(These values will be used to create lower and upper_bounds  ',
             ' used in the optimization procedure)'),
        if (failed)
          div(tags$b("Invalid params.", style = "color: red;")),

        footer = tagList(
          modalButton("Cancel"),
          actionButton("prior_onegrow_migok", "OK")
        )
      )
    }    
            
    observeEvent(input$prior_onegrow_migok, {
      # Check the params are ok.
      if (!is.null(input$prior_onegrowMaxIter) ) {
        Modelparams$data <- list(
            lower_bounds = c(input$nu1F[1], input$nu2B[1], input$nu2F[1], input$m[1] , input$Tp[1],input$T[1]), 
            upper_bounds = c(input$nu1F[2], input$nu2B[2], input$nu2F[2], input$m[2] , input$Tp[2],input$T[2]),
            p0 = c(input$nu1F0, input$nu2B0, input$nu2F0, input$m0, input$Tp0, input$T0),
            OptiMethod  = input$OptiMethod,
            MaxIter     = input$prior_onegrowMaxIter,
            model       = "prior_onegrow_mig"
            )
        removeModal()
      } else {
        showModal(prior_onegrow_migParams(failed = TRUE))
      }
    })

