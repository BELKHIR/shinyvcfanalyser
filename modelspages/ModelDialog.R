    ModelParams <- function(model, descr, failed = FALSE) {  
      limits = list(a_max = 10.0 , a_min = 0.1 , a_start = 1.0 ,#beta-binomial alpha
        b_max = 10.0 , b_min = 0.1 , b_start = 1.0 ,#beta-binomial beta
        nu1_max = 20.0 , nu1_min = 1e-3 , nu1_start = 1.0 ,#current size of pop1
        nu2_max = 20.0 , nu2_min = 1e-3 , nu2_start = 1.0 ,#current size of pop2
        nu1e_max = 20.0 , nu1e_min = 1e-3 , nu1e_start = 1.0 ,#current size of pop1
        nu2e_max = 20.0 , nu2e_min = 1e-3 , nu2e_start = 1.0 ,#current size of pop2

        ne1_max = 20.0 , ne1_min = 1e-3 , ne1_start = 1.0 ,#size of pop1 after exponential change
        ne2_max = 20.0 , ne2_min = 1e-3 , ne2_start = 1.0 ,#size of pop2 after exponential change
        Ts_max = 12.0 , Ts_min = 1e-10 , Ts_start = 1.0 ,#duration of the period after split
        Ti_max = 12.0 , Ti_min = 1e-10 , Ti_start = 1.0 ,#duration of the second period
        m12_max = 20.0 , m12_min = 1e-10 , m12_start = 1.0 ,#migration from 2 to 1
        m21_max = 20.0 , m21_min = 1e-10 , m21_start = 1.0 ,#migration from 1 to 2
        nr_max = 1.00 , nr_min = 1e-10 , nr_start = 0.5 ,#fraction of the genome with reduced Ne
        bf_max = 1.00 , bf_min = 1e-3 , bf_start = 0.5 ,#factor of Ne reduction
        P_max = 1.00 , P_min = 1e-10 , P_start = 0.5 ,#1 - fraction of the genome with Me=0
        O_max = 1.0 , O_min = 1e-10 , O_start=0.8 #fraction of SNPs accurately oriented relative to the outroup
        )  
      models_params = read.table("modelspages/models_params.txt", header=T, sep="\t")
      curModel <<- model
      curmodelParams <<- models_params[models_params$Model == model,]
      nbparamsModel <<- nrow(curmodelParams)

      pvars <- nbparamsModel
      modalDialog(title=paste0(model, ": ", descr),      
       if (pvars > 0) {
       div(
        lapply(seq(pvars), function(i) {
            param = curmodelParams$lib[i] 
            mini = limits[[paste0(param,"_min")]]
            maxi = limits[[paste0(param,"_max")]]
            start = limits[[paste0(param,"_start")]]
            fluidRow( width =largeur*2,style = "height:100px;",
            column(width = largeur, tags$div(title=paste0(curmodelParams$help[i], ": lower & upper bounds"),sliderInput(curmodelParams$lib[i], label =curmodelParams$lib[i] , min = mini, max = maxi, value = c(mini, maxi + 1)))),
            column(width = largeur, tags$div(title=paste0(curmodelParams$help[i], ": initial value"),numericInput(paste0(curmodelParams$lib[i],"0"), "Start value", value = start, min = mini, max = maxi + 1)))
            )
        }),
        # fluidRow(column(width = largeur,selectInput("OptiMethod", "Optimization Method", choices =c("optimize","optimize_log","optimize_log_lbfgsb","dual_anneal"))),
        # column(width = largeur,sliderInput("MaxIter", label = "How long the optimizer will run", min = 1, max = 100, value = 1))),
  
        #tags$head(tags$style(".row{height:150px;}")),
        
        if (failed)
          div(tags$b("Invalid params.", style = "color: red;"))
       )
      },
      footer = tagList( modalButton("Cancel"), actionButton("IMok", "OK")  )
      )     
    }      

    observeEvent(input$IMok, {
      # Check the params are ok.
      #if (!is.null(input$MaxIter) ) { have to decide wich value to check for validation
      if(TRUE){    
        lower_bounds <- rep(NA, length(nbparamsModel))
        upper_bounds <- rep(NA, length(nbparamsModel))
        p0 <- rep(NA, length(nbparamsModel))

        for(k in 1:nbparamsModel) { 
            inputName <- curmodelParams$lib[k]
            # only get a value if the  input exists
            if (!is.null(inputName))
            {
            lower_bounds[k] <- input[[inputName]][1]
            upper_bounds[k] <- input[[inputName]][2]
            }  
            startValInput =  paste0(inputName,"0")
            if (!is.null(startValInput) ) p0[k] <- input[[startValInput]]
        }
        
        Modelparams$data <- list(
            lower_bounds = lower_bounds, 
            upper_bounds = upper_bounds,
            p0 = p0,
            # OptiMethod  = input$OptiMethod,
            # MaxIter     = input$MaxIter,
            model       = curModel 
            )
        removeModal()
      } else {
        showModal(ModelParams(curModel, failed = TRUE))
      }
    })