
library(shinydashboard)
library(shiny)
library(threejs)
library(SNPRelate)
library(SeqArray)
library(Rtsne)
library(RColorBrewer)
#library(gplots)
library(ComplexHeatmap)

library("ggtree")
library(dplyr)
library(calibrate)

require(pairsD3)
library(crosstalk)
library(d3scatter)
library(ggplot2)
library(ggalt)
library(gridExtra)
library(plyr)
library(shinycssloaders)
library(shinyFiles)
library(pophelper)
library(ape)
library(snapper)
options(shiny.maxRequestSize=250*1024^2) 
library(patchwork)
library(reticulate)
library(shinyWidgets)

  use_python("/usr/bin/python3")
  
  moments = import("moments", convert=FALSE)
  momentsMisc = import("moments.Misc",convert = FALSE)
  momentsSpectrum = import("moments.Spectrum_mod",convert = FALSE)
  #source_python("momentsmodels/momentsCustomDemographics2D.py")
  source_python("momentsmodels/moments_models_2pop.py")
  source_python("local_Misc.py",convert = FALSE)
  myInference = import("moments.moments_inference_dualanneal", convert=FALSE) 

  dadi = import("dadi",convert = FALSE) 
  # dadiSpectrum = import("dadi.Spectrum_mod",convert = FALSE)
  dadiMisc = import("dadi.Misc",convert = FALSE)
  # source_python("dadimodels/dadiCustomDemographicModels.py")

  plt = import("matplotlib.pyplot")


read_data <- function(fic, ficname, maxLD, minMAF, maxMissing){
    
    outdir <- tempfile()
    outimage = paste0(outdir,"/bcftools_mqc.png")

    outstats = tempfile()
    # Generate the PNG
    cmd = paste0("bcftools stats --threads 8 -s - ", fic, "  > ", outstats )
    ttt = system(cmd, intern=T)

    cmd = paste0("plot-vcfstats -P -p ", outdir, " -t '",basename(ficname)," maxLD=",maxLD," minMAF=",minMAF," maxMissing=", maxMissing,"' ",  outstats )
    ttt = system(cmd, intern=T)

    cmd = paste0("convert -size 1200x1000 \\( " ,outdir, "/depth.0.png  ",outdir,"/substitutions.0.png ", outdir,"/hwe.0.png +smush +400 \\) \\
          \\( ",outdir,"/snps_by_sample.0.png ",outdir,"/indels_by_sample.0.png +smush +150 \\) \\
         \\( ",outdir,"/dp_by_sample.0.png ", outdir,"/singletons_by_sample.0.png +smush +150 \\) \\
          \\( ",outdir,"/hets_by_sample.0.png  ",outdir,"/tstv_by_sample.0.png   +smush +150 \\) \\
          -append   ", outimage)

    ttt = system(cmd, intern=T)
    
    cmd = paste0("grep 'number of samples:' ", outstats, " | cut -f4 ")
    samples = system(cmd, intern=T) 

    cmd = paste0("grep 'number of records:' ", outstats, " | cut -f4 ")
    records = system(cmd, intern=T) 

    cmd = paste0("grep 'number of SNPs:' ", outstats, " | cut -f4 ")
    snps = system(cmd, intern=T)

    cmd = paste0("grep 'number of indels:' ", outstats, " | cut -f4 ")
    indels = system(cmd, intern=T)

    return(list(samples = samples, records = records, snps = snps, indels= indels, bcfstatsimg= outimage))
}

body <- dashboardBody(height = '1200px',
  shinyjs::useShinyjs(),

  fluidRow(id = 'page',
    load_snapper(),
    
    shinyjs::extendShinyjs(text = '
          shinyjs.selectInput_tooltips = function(params){

          var defaultParams = {
            id : null,
            tooltips : null
          };
          params = shinyjs.getParams(params, defaultParams);

          var selectInput = $("#"+params.id).closest("div").find(".dropdown-menu").get(1);
          var element_selectInput = selectInput.childNodes;

          if(element_selectInput.length >0 && element_selectInput[0].title == ""){ // to be trigger only once
            for(var i = 0; i< element_selectInput.length; i++){
              element_selectInput[i].title = params.tooltips[i];
            }
          }

        }; 
      ',functions=c("selectInput_tooltips") ),
  
    tabBox(title = "VCF analysis", width=12, height = "1200px", id = "tabset1", 
      tabPanel("Summary", "You can apply filter to redraw summary stats", height = "900px",
        fluidRow(
           tags$head(tags$style(HTML(".small-box {height: 100px}"))),
           valueBoxOutput("nbSamples", width = 3),
           valueBoxOutput("nbRecords", width = 3),
           valueBoxOutput("nbSnps", width = 3),
           valueBoxOutput("nbIndels", width = 3)
        ),
        fluidRow( 
          # Use imageOutput to place the image on the page
          box(width=12, #this snippet generates a progress indicator for long anlyses
          div(class = "busy",  id="busy-img",
              p("In progress..."), 
              img(src="8puiO.gif", height = 100, width = 100,align = "center")
          ), withSpinner(imageOutput("bcftoolstats",  height = '1000px') ) )
        )
      ),
      tabPanel("Scan-geno",
      fluidRow(   
       box( width = 2, numericInput("chr_geno_plot", "Chr or Ctg", value = 1, min = 1, max =50,step = 1) ,solidHeader=TRUE),
       box( width = 1, numericInput("window", "Window", value = 1, min = 1, max = 10000,step = 1) ,solidHeader=TRUE),
       box( width = 2, numericInput("maxSnp_geno_plot", "max nb snps per window", value = 600, min = 10, max = 1000,step = 10),solidHeader=TRUE ),
       box( width = 7, radioButtons("codage", label = ("How to code genotype"), choices = list("REF/REF=2 REF/ALT=1 ALT/ALT=0" = 1, "REF/REF=0 REF/ALT=1 ALT/ALT=2 (biallelic only!)" = 2), selected = 2, inline=T),solidHeader=TRUE),
      tabBox(width=12,
      tabPanel("Genotypes of selected region using REF base counts", box(width=12, withSpinner(plotOutput("genoPlot", width="100%", height = '1000px') ))  ),
      tabPanel("Missingness for pairwise samples ordered by IBS",  box(width=12, withSpinner(plotOutput("missingPlot", width="100%", height = '1000px') )) ),
      tabPanel("FstScan",  box(width=12, withSpinner(plotOutput("SnpFstPlot", width="100%", height = '1000px') )) )
      )
      )
      ),
      tabPanel("PCA ",  height = "1200px",
       tabBox(width=12,
        tabPanel("PCA 2D scatterplot",
        fluidRow(
            box( width = 2, actionButton("runpca", "Go"), downloadButton("PCAdownloadData", "Download")) ,
            box( width = 3, selectInput('xcol', 'PCA First Axis', 1:10)),
            box( width = 3, selectInput('ycol', 'PCA Second Axis', 1:10,  selected=2) )
        ),
        fluidRow(
            column( 6, withSpinner(d3scatterOutput("Scatterplot",height = "700px" ) ) ),
            column( 6, plotOutput("Eigen",height = "400px" ), plotOutput("by_pop") ) 
        )
       ),
       tabPanel("PCA and Hclust scatterplot",'The PCA coordinate are clustered with the R hclust function ',
        fluidRow(
            box( width = 3,  selectInput('xcolhc', 'PCA First Axis', 1:10)),
            box( width = 3,  selectInput('ycolhc', 'PCA Second Axis', 1:10,  selected=2) ),

            box( width = 2, checkboxInput("pcashowLabels", "Show samples labels", value=F) ),
            box( width = 2, numericInput("pcahclustK", "Desired number of clusters from hclust res", value = 3, min = 2, max = 10,step = 1) ),
            box( width = 2, checkboxInput("pcashowClusters", "Delineate clusters", value=T) )
        ),
        fluidRow(
            box( width = 12, plotOutput("PcaHclust",height = "1000px" ) )
        )
       ),
       tabPanel("PCA 3D scatterplot", "3D plots with 3 PCA axis",
        fluidRow(
            box( width = 4, selectInput('xcol3d', 'PCA First Axis', 1:10)),
            box( width = 4, selectInput('ycol3d', 'PCA Second Axis', 1:10,  selected=2) ),      
            box( width = 4, selectInput('zcol3d', 'PCA Third Axis', 1:10,  selected=3)  )
        ),
        fluidRow(
            box( width=12, scatterplotThreeOutput("Scatterplot3d", width="100%", height = "900px" ) )
        )
      ),
      tabPanel("PCAPairewise Scatterplots", "Pairewise Scatterplots with the top 10 principal components", 
        fluidRow(
              box( width = 4, numericInput('nbaxes', 'Number of PCA Axis to plot', value=5, min=3, max=10, step=1))
        ),
        fluidRow(
              box(width=12, withSpinner(plotOutput("Scatterplotpairs", width="100%", height = '900px') ) )
        )
      )
      )),
      tabPanel("t-SNE", "t-SNE: a method for constructing a low dimensional embedding of high-dimensional data using Rtsne package (default params). https://cran.r-project.org/web/packages/Rtsne/index.html.", 
          fluidRow(   
            box( width = 2, numericInput("perplexity", "Perplexity", value = 5, min = 1, max =50,step = 1) ),
            box( width = 1, actionButton("runTSNE", "Go")) ,
            box( width = 3, checkboxInput("showLabels", "Show samples labels", value=F) )
          ),
          fluidRow( box(width=12, withSpinner(plotOutput("tsnePlot", height = '1000px') )) )
      ),
      tabPanel("Pairewise Fst", "Pairewise Fst estimate from Weir and Cockerham's 1984 paper. *Above digonal: mean Fst estimate  *Under diagonal: weighted Fst estimate", 
          fluidRow(   
            box( width = 4, actionButton("runfst", "Go")) 
          ),
          fluidRow( box(width=12, shinycssloaders::withSpinner(plotOutput("treePlot") ))),
          fluidRow( box(width=12, DT::dataTableOutput("Fstpairs", height = '1200px') ))
      ),
      tabPanel("IBS pairwise identities", "genome-wide average Identity-By-State (IBS) pairwise identities", 
        tabBox(width=12,
          tabPanel("Heatmap",
          fluidRow(   
            box( width = 4, actionButton("runibs", "Go")) 
          ),
          fluidRow( box(width=12, withSpinner(plotOutput("ibs",  height = '1200px') ) ) )
          ),
          tabPanel("Hclustering and automatic grouping",
            fluidRow( box(width=12, withSpinner(plotOutput("ibsClusters",  height = '1200px') ) ) )
          )
        )
      ),
      tabPanel("fastStructure", "algorithm for inferring population structure from large SNP genotype data. It is based on a variational Bayesian framework for posterior inference", 
        fluidRow(
            box( width = 3, sliderInput("K", label = "K: number of populations range ", min = 2, max = 10, value = c(2, 4)) ),
            box( width = 3, numericInput("fastcv", "test sets for cross-validation", value = 0, min = 0, max = 5,step = 1) ),
            box( width = 3, selectInput('Prior', 'Prior', c("simple","logistic")) ),
            box( width = 3, actionButton("runfastStructure", "Go")) 
             ),
            box(width = 9,div(style = 'overflow-x: scroll', withSpinner(imageOutput("distruct", width="1200px", height="800px")) )),
            box(width = 3, selectInput("kfile", "Mean of admixture proportions for k= :", choices = c("K=2"="2", "K=3"="3", "K=4"="4")),
            downloadButton("downloadData", "Download"),
            DT::dataTableOutput("meanQtable",  height="700px" )
            ) 
      ),
      tabPanel("Pairwise frequency Spectrum", 
        tabBox(width=12,
          tabPanel("Plot pairwise FS", "Use moments to extract a frequency spectrum for each population pair from the VCF file ",
          fluidRow( 
            tags$head(tags$style(HTML('
            .selectize-input {white-space: nowrap}
            #firstPop+ div>.selectize-dropdown{width: 200px !important;}
            #firstPop+ div>.selectize-input{width: 200px !important; }
            #downFirst+ div>.selectize-dropdown{width: 200px !important;}
            #downFirst+ div>.selectize-input{width: 200px !important; }
            #secondPop+ div>.selectize-dropdown{width: 200px !important; }
            #secondPop+ div>.selectize-input{width: 200px !important; }
            #downSencond+ div>.selectize-dropdown{width: 200px !important;}
            #downSencond+ div>.selectize-input{width: 200px !important; }
            
            '))),  
            box(width = 2, tags$div(title="Click here to create a FS from your VCF",actionButton("runsfs", "Create spectrum")), checkboxInput("fold","Data is polarized", value=FALSE )) ,
            box(width = 2, selectInput("firstPop", "First pop", choices =character(0))),
            box(width = 2,numericInput("downFirst", "down sample first pop to", value = 20, min=2, max=100)),
            box(width = 2, selectInput("secondPop", "Second pop", choices = character(0))),
            box(width = 2,numericInput("downSencond", "down sample second pop to", value = 20, min=2, max=100)),
            box(width=2, tags$div(title="Stop reading VCF when this amount of memory is reached",numericInput("maxMemPercent", "Max % memory", value = 50, min=2, max=80)))
          ),
          fluidRow( box(width=12, withSpinner(imageOutput("sfsPlot",  height = '800px') ) ) )
          ),
          tabPanel("moments Demographic Inference",
          fluidRow( 
            tags$head(tags$style(HTML('
            .selectize-input {white-space: nowrap}
            #firstPop1+ div>.selectize-dropdown{width: 200px !important;}
            #firstPop1+ div>.selectize-input{width: 200px !important; }
            #downFirst1+ div>.selectize-dropdown{width: 200px !important;}
            #downFirst1+ div>.selectize-input{width: 200px !important; }
            #secondPop1+ div>.selectize-dropdown{width: 200px !important; }
            #secondPop1+ div>.selectize-input{width: 200px !important; }
            #downSencond1+ div>.selectize-dropdown{width: 200px !important; }
            #downSencond1+ div>.selectize-input{width: 200px !important; }
            #Model+ div>.selectize-dropdown{width: 200px !important; }
            #Model+ div>.selectize-input{width: 200px !important; }
            '))),  
            box(width = 2, selectInput("firstPop1", "First pop", choices =character(0)),
            numericInput("downFirst1", "down sample to", value = 20, min=2, max=100)),
            box(width = 2, selectInput("secondPop1", "Second pop", choices = character(0)),
            numericInput("downSencond1", "down sample to", value = 20, min=2, max=100)),
            #box(width = 2, selectInput("Model", "Demgraphic Model", choices =c("SI","IM","SC","prior_onegrow_mig")),
            box(width = 2, pickerInput("Model", "Demgraphic Model", choices =c("SI","IM","SC","prior_onegrow_mig")),
            actionButton("setParams", "Set Model params")),
            box(width=2,selectInput("OptiMethod", "Optimization Method", choices =c("optimize","optimize_log","optimize_log_lbfgsb","dual_anneal")), div(style="height: 35px;",sliderInput("MaxIter", label = "How long the optimizer will run", min = 1, max = 100, value = 1))), 
            
            box(width=2,numericInput("generationTime", "Generation Time (Y)", value = 20, min = 1, max = 100,step = 1), numericInput("mutationRate", "Mutation rate (%)", value = 2, min = 1, max = 10,step = 1)),
            box(width = 2, actionButton("runInference", "Run Inference with selected Model"), checkboxInput("foldInference","Data is polarized", value=FALSE ))

          ),
          tabBox(title = "Graphs", width=12,  id = "tabsetModel", 
            tabPanel("Comparison", "plot comparisons between 2D FS for models and data", 
          fluidRow( box(width=12, withSpinner(imageOutput("demoPlot",  height = '800px'))))
            ),
            tabPanel("ModelDiagram", "schematic representation of the selected model", 
            fluidRow( box(width=12, withSpinner(imageOutput("PlotModel",  height = '800px'))))),
            tabPanel("Optimization", "Log-likelihood of the data given the model during the optimization procedure",
            fluidRow( box(width=12, withSpinner(imageOutput("logLikPlot",  height = '800px')))))
          )
          
          )
         )
    )
  )
  
)
)

shinyApp(
  ui = dashboardPage(
    dashboardHeader(title = "Analyse VCF file"),
    dashboardSidebar(
    
    # add a download link for the main panel
    snapper::download_link(
      ui = '#page',
      label = '',
      filename = 'main_panel.png'
    ),
        shinyFilesButton("servervcffile" ,label="Select a VCF in the server", title="", multiple=FALSE),
       # h5("Or"),
       # fileInput("vcffile", "Upload a gzipped VCF File", accept = c( "text/vcf",".vcf.gz") ),
        br(),
        shinyFilesButton("serverpopMapfile" ,label="Select popMap file in the server", title="", multiple=FALSE),
        br(),
        # h5("Or"),
        # fileInput("popMap", "Upload population Map", accept = c( "txt/tsv", "text/tab-separated-values,text/plain", ".tsv")),
        numericInput("minMAF", "min MAF", value = 0.0, min = 0, max = 1,step = 0.01) ,
        numericInput("maxMissing", "max missing snps(1 allows sites completely missing 0 no missing data allowed)", value = 1, min = 0, max = 1,step = 0.01) ,
        numericInput("maxLD", "Recursively removes SNPs with LD greater than", value = 1.0, min = 0, max = 1,step = 0.01) ,
        actionButton("applyFilter", "Apply filter"), 
        br(),
        h5("This tool is aimed at helping to analyse a small to medium vcf files of genotyped individuals"),
        br(),
        shinyDirButton("saveFilteredVCF", "Save filtered VCF", "Please select a folder")

    ),
    body
    
  ),

server = function(input, output, session) {

ficname = ""
gds_file = ""
curr.sample.id = 1:10
ibs = NULL
ibsDendro = NULL

filteredFile  = ""
globalmaxMissing = 0
globalminMAF     = 0
maxLD = 1
palette = "Set3" #"Accent" #"Spectral"
curModel <- "SI"
nbparamsModel <- 4
curmodelParams=data.frame()

shinyFileChoose(input, "servervcffile", root=c(Data="/Data",Results="/Results"),filetypes=c('vcf', 'gz'), session = session)
shinyFileChoose(input, "serverpopMapfile", root=c(Data="/Data",Results="/Results"),filetypes=c('txt', 'tsv','csv'), session = session)
#shinyFileSave(input,"saveFilteredVCF", root = c(Data="/Data",Results="/Results"),filetypes=c('vcf', 'gz'), session = session)
shinyDirChoose(input, "saveFilteredVCF", roots = c(Data="/Data",Results="/Results"), session = session,  allowDirCreate = TRUE)
Modelparams <- reactiveValues(data = NULL)
doInference <- FALSE
shinyjs::hide("busy-img")
shinyjs::disable("saveFilteredVCF")
shinyjs::disable("applyFilter")
# Return the UI for a modal dialog with data selection input. If 'failed' is
# TRUE, then display a message that the previous value was invalid.
largeur = 6

source("modelspages/ModelDialog.R", local=T)

models_params = read.table("modelspages/models_params.txt", header=T, sep="\t")
availableModels = unique(models_params$Model)
updatePickerInput(session, "Model",
      choices  = availableModels,
      selected = curModel
    )
models_tooltips = read.table("modelspages/listeModels.txt", header=F, sep="\t")    
rownames(models_tooltips) = models_tooltips[,1]
models_tooltips <<- models_tooltips[match(availableModels,models_tooltips[,1]),]

shinyjs::onclick("Model" ,shinyjs::js$selectInput_tooltips("Model",models_tooltips[,2])) 

dummy = list(tableau = data.frame(sample.id = 1,
                    pop = factor(1),
                    Axe1 = 0,    # the first eigenvector
                    Axe2 = 0,    # the second eigenvector
                    Axe3 = 0,    # the 3 eigenvector
                    Axe4 = 0,    # the 4 eigenvector
                    Axe5 = 0,
                    Axe6 = 0,
                    Axe7 = 0,
                    Axe8 = 0,
                    Axe9 = 0,
                    Axe10 = 0,
                    couleur = 'black',
                    stringsAsFactors = FALSE) ,
              valpropre = 0,
              ibs = NULL ,
              popc = 'black'
              )

summ = list(depth=NULL,minmafF=0,maxmissingF=0, missingness=NULL, maf=NULL, sitemissingness = NULL  )
summary <- reactiveValues()
summary(summ)

selectedData <- reactive({
    if (is.null(input$vcffile))
    {
        fics = parseFilePaths(c(Data="/Data",Results="/Results"),input$servervcffile)
        if (nrow(fics)>0) {
          vcf.fn = fics$datapath[1]
          ficname = vcf.fn
        }
        else return(NULL)
    }
    else
    {
      if (input$vcffile$datapath == "") return(NULL)
      vcf.fn = input$vcffile$datapath
      ficname = input$vcffile$name
    }

    filteredFile <<- vcf.fn #sans apply filter filteredFile est = au fichier uploader (ou choisi)
    #create seqGDS
    
    original_gds_file <<- paste0("/tmp/",basename(vcf.fn),".gds")
    gds_file <<- original_gds_file 
    seqVCF2GDS(vcf.fn, gds_file, verbose=FALSE, storage.option="LZ4_RA", parallel = 8,raise.error=F )
   
    # Create bed file for fastStructure
    genofile <- seqOpen(gds_file)
    curr.sample.id <<- seqGetData(genofile, "sample.id")

    #calc IBS over all SNP (this can be used for correlation with missingness over widows)
    ibs <<- snpgdsIBS(genofile,  snp.id=NULL, num.thread=8, autosome.only=FALSE)
    # dendrogram on IBS values:
    ibsDendro <<- snpgdsHCluster(ibs, need.mat = FALSE)$hclust

    #il faut upgrader la version rocker/r-ver à 4 pour avoir des versions plus recentes de SeqArray SNPRelate ...
    seqGDS2BED(genofile, out.fn=filteredFile,verbose=F)
    seqClose(genofile)

    resu = read_data(vcf.fn,ficname,1,0,0) 
    summary$minmafF = 0
    summary$maxmissingF=0

    summary$samples = resu$samples
    summary$snps = resu$snps
    summary$records = resu$records
    summary$indels = resu$indels
    summary$bcfstatsimg = resu$bcfstatsimg
    shinyjs::enable("applyFilter")
    return(summary)

})

popDeco <-reactive({
 # print("Inside popDeco")
  pop_color1= rep("grey", length(curr.sample.id))
  pop_code1 = rep("pop", length(curr.sample.id))
  names(pop_code1) = curr.sample.id
  popc="grey"
  names(popc)="pop"

  fics = parseFilePaths(c(Data="/Data",Results="/Results"),input$serverpopMapfile)
  if (nrow(fics)>0) {        
      ficpopmap = fics$datapath[1]
  
      popmap = read.table(ficpopmap, header=F, sep='\t', stringsAsFactors=F)
      if  (length(curr.sample.id) == length(popmap[,1]) & length(setdiff(curr.sample.id, popmap[,1])) == 0)
      {
        rownames(popmap) = popmap[,1]
        pop_code1 = popmap[popmap[,1] %in% curr.sample.id, 2]
        names(pop_code1)  = curr.sample.id
        ncolors = length(unique(pop_code1))
        popc <- colorRampPalette(brewer.pal(12, "Set3"))(ncolors) #brewer.pal(n = ncolors, name = palette)
        names(popc)       = unique(pop_code1)
        pop_color1        = popc[pop_code1]
        names(pop_color1) = curr.sample.id  
        
      }  else {
      showModal(modalDialog(title = "Warning", "The list of samples in popMap is different from samples in the VCF. Population info is ignored. !"))
      }
  }

  popDeco=list(pop_color1=pop_color1, popc = popc, pop_code1= pop_code1)
})

observeEvent(input$applyFilter,{
    globalmaxMissing <<- isolate(input$maxMissing)
    globalminMAF     <<- isolate(input$minMAF)
    maxLD            <<- isolate(input$maxLD)
    showModal(modalDialog(
        title = "Warning",
        "The filters will be applied to summary stats ! \n You have to re-run the other analysis to apply these filters"
      ))

    shinyjs::show("busy-img")

    if (is.null(input$vcffile))
    {
        fics = parseFilePaths(c(Data="/Data",Results="/Results"),input$servervcffile)
        if (nrow(fics)>0) {
          vcf.fn = fics$datapath[1]
          ficname = vcf.fn
        }        else return(dummy)
    }
    else
    {
      vcf.fn = input$vcffile$datapath
      ficname = input$vcffile$name
    }

    
    #filter
    filteredFile <<- paste0("/tmp/",basename(vcf.fn),"filtered-maxLD",isolate(input$maxLD),"-missing",isolate(input$maxMissing),"-maf",isolate(input$minMAF),".recode.vcf.gz") 
    if (! file.exists(filteredFile) ) 
    { 
      if (input$maxLD == 1 && input$maxMissing == 1 && input$minMAF==0)
       filteredFile <<- vcf.fn
      else 
      {
      genofile <- seqOpen(original_gds_file)
      snpset <- snpgdsLDpruning(genofile,autosome.only=FALSE, maf = input$minMAF, missing.rate=input$maxMissing, ld.threshold=input$maxLD, verbose=F)
          
      # get SNP ids
      snpset.id <- unlist(snpset, use.names=FALSE)

      #filter the snpset
      seqSetFilter(genofile, variant.id=snpset.id)

      gds_file <<- paste0("/tmp/",basename(vcf.fn),"-maxLD", maxLD,"-maxmiss",globalmaxMissing,"-maf",globalminMAF,".gds")
      seqExport(genofile, gds_file,verbose=F)

      # convert back to vcf
      seqGDS2VCF(genofile, filteredFile, verbose=F)

      seqClose(genofile)

      # Create bed file for fastStructure
      genofile <- seqOpen(gds_file)
      #il faut upgrader la version rocker/r-ver à 4 pour avoir des versions plus recentes de SeqArray SNPRelate ...
      seqGDS2BED(genofile, out.fn=filteredFile,verbose=F)
    
      #calc IBS over all SNP (this can be used for correlation with missingness over widows)
      ibs <<- snpgdsIBS(genofile,  snp.id=NULL, num.thread=8, autosome.only=FALSE)
      # dendrogram on IBS values:
      ibsDendro <<- snpgdsHCluster(ibs, need.mat = FALSE)$hclust
      seqClose(genofile)
      
      }
    }
    shinyjs::hide("busy-img")
    shinyjs::enable("saveFilteredVCF")
    resu = read_data(filteredFile, ficname, isolate(input$maxLD), isolate(input$minMAF), isolate(input$maxMissing)) 
    
    summary$minmafF = isolate(input$minMAF)
    summary$maxmissingF=isolate(input$maxMissing)

    summary$samples = resu$samples
    summary$snps = resu$snps
    summary$records = resu$records
    summary$indels = resu$indels
    summary$bcfstatsimg = resu$bcfstatsimg
    
})    

observe({
   if (! is.integer(input$saveFilteredVCF)) {
      sel_path=parseDirPath(c(Data="/Data",Results="/Results"), input$saveFilteredVCF)
      fileDest = paste0(sel_path,"/",basename(filteredFile))
      file.copy(filteredFile,fileDest)
      showModal(modalDialog(title = "Message", paste0("File saved to : ", fileDest)))
   }
})
 

output$bcftoolstats <- renderImage({
  if (is.null(selectedData()) ) {
      outimage = tempfile()
      png(outimage, width = 1200, height = 1000)
      dev.off()
    } else
    {
      outimage = selectedData()$bcfstatsimg
    }
    # Return a list containing the filename
    list(src = outimage,
         contentType = 'image/png',
         width = "100%",
         height = "100%",
         alt = "Bcftools stats")
  }, deleteFile = TRUE)

output$nbSamples <- renderValueBox({
    if (is.null(selectedData()) ) valueBox("?", "Samples", icon = icon("barcode", lib = "font-awesome"), color = "yellow"    )
    valueBox( selectedData()$samples, "Samples", icon("barcode", lib = "font-awesome"), color = "yellow"    )
  })

output$nbRecords <- renderValueBox({
    if (is.null(selectedData()) ) valueBox("?",  "Records", icon = icon("bars", lib = "font-awesome"), color = "yellow" )
    valueBox( selectedData()$records, "Records", icon = icon("bars", lib = "font-awesome"), color = "yellow"  )
  })

output$nbSnps <- renderValueBox({
    if (is.null(selectedData()) ) valueBox("?",  "SNP", icon = icon("bars", lib = "font-awesome"), color = "yellow" )
    valueBox( selectedData()$snps, "SNP", icon = icon("bars", lib = "font-awesome"), color = "yellow"  )
  })

output$nbIndels <- renderValueBox({
    if (is.null(selectedData()) ) valueBox("?",  "IDEL", icon = icon("bars", lib = "font-awesome"), color = "yellow" )
    valueBox( selectedData()$indels, "INDEL", icon = icon("bars", lib = "font-awesome"), color = "yellow"  )
  })  


tsneRun <- eventReactive(input$runTSNE,{
    if (is.null(input$vcffile))
    {
        fics = parseFilePaths(c(Data="/Data",Results="/Results"),input$servervcffile)
        if (nrow(fics)>0) vcf.fn = fics$datapath[1]
        else return(dummy)
    }
    else
    { 
      if (! is.null(input$vcffile$datapath) )      vcf.fn = input$vcffile$datapath
      else return(dummy)
    }
    
    if (! file.exists(gds_file))
    {
      seqVCF2GDS(filteredFile,gds_file, verbose=TRUE,  storage.option="LZ4_RA", parallel = 8)
    }
    genofile <- seqOpen(gds_file)
   
    # for now problem with missing data
    infos <- seqGetData(genofile, c("genotype","sample.id") )
    seqClose(genofile) 
  #this is for bi-allelic snp codage in 0 (ref/ref) , 1 (ref/alt) and 2 (alt/alt)
  if( isolate(input$codage) == "2")
  {
    geno_matrix <- infos[[1]][1,,] + infos[[1]][2,,]
  }else{
  #this is for multi-allelic snp codage in 2 (ref/ref) , 1 (ref/alt) and 0(alt/alt)
    geno_matrix <- (infos[[1]][1,,] == 0) + (infos[[1]][2,,] == 0)
  }
  samples.id = infos[[2]]
  rm(infos)
  gc(TRUE)
    rownames(geno_matrix) = samples.id
    tsne_out  <- Rtsne(geno_matrix, pca=TRUE, perplexity=input$perplexity, theta=0.5, dims=2, num_threads = 4)
    
    rm(geno_matrix) 
    gc(TRUE) 
    popDeco = popDeco() 
    tsneresu <-list(tsneRes = tsne_out ,              
              couleurs = popDeco()$pop_color1,
              pops = popDeco()$pop_code1,
              popc = popDeco()$popc,
              samples = samples.id,
              minMAF = input$TSENminMAF,
              maxMissing = input$TSNEmaxMissing
            )          
})

plot_cluster=function(data, var_cluster,  couleurs, shape)
{ #https://luisdva.github.io/rstats/Grouping-points/
  # find_hull <- function(df) df[chull(df$axe1, df$axe2), ]
  # hulls <- ddply(data, which_hull, find_hull)
  
  ggplot(data, aes_string(x="axe1", y="axe2", color=var_cluster, shape=shape)) + ggtitle(paste0("maxLD=",maxLD," maxMissing=", globalmaxMissing," minMaf=", globalminMAF)) +
  geom_point(size=3) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("") +
  #theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.direction = "horizontal", 
        legend.position  = "bottom",
        legend.box       = "horizontal") + 
  scale_colour_manual(values=couleurs)
  #scale_colour_brewer(palette = palette) 
}

genoWindow <- reactive({

  window= (input$window)
  wsize = (input$maxSnp_geno_plot)

  #if (! file.exists(gds_file)) {print("No GDS file!"); return(NULL) }

  showfile.gds(closeall=TRUE)
  
  if (! file.exists(gds_file)) {print("No GDS file!"); }
 
  genofile <- seqOpen(gds_file)

  l=seqGetData(genofile, c("chromosome",  "variant.id"))
  df <- data.frame(matrix(unlist(l), ncol=length(l), byrow=F), stringsAsFactors=F)
  colnames(df) = c("ctg","ident")
  
  #params$chr_geno_plot can be a contig name or a contig index in the list of ctgs (dans quel ordre ?!!!)
  if ((input$chr_geno_plot) %in% df$ctg) {
    chr = input$chr_geno_plot
  }else {
    chr = unique(df$ctg)[input$chr_geno_plot]
  }
  
  snps = (df %>% filter(ctg == chr ))$ident
  totsnps=length(snps)
  first= ((window-1) * wsize) + 1
  if(first >= totsnps){
    first=totsnps - wsize
    showModal(modalDialog(title = "Warning", "You reached the end of the contig !"))
    return(NULL)
  }
  last=first+wsize -1
  if (last > totsnps ){last=totsnps}
  
  seqSetFilter(genofile, variant.id=snps[first:last])
  infos = seqGetData(genofile, c("sample.id", "position", "variant.id","genotype"))
  seqClose(genofile)
  #geno <- infos[[4]]
  firstSnpPos=infos[[2]][1]
  lastSnpPos=infos[[2]][length(infos[[2]])]

  #this is for bi-allelic snp codage in 0 (ref/ref) , 1 (ref/alt) and 2 (alt/alt)
  if( input$codage == "2")
  {
    geno_matrix <- infos[[4]][1,,]+ infos[[4]][2,,]
    keylabels=c("REF/REF","REF/ALT","ALT/ALT","Missing")
    my_palette =c("#d4b9da","#e7298a","#980043")#,"#FFFFFF")
  }else{
  #this is for multi-allelic snp codage in 2 (ref/ref) , 1 (ref/alt) and 0(alt/alt)
    geno_matrix <- (infos[[4]][1,,] == 0) + (infos[[4]][2,,] == 0)
    keylabels=c("ALT/ALT","REF/ALT","REF/REF", "Missing")
    my_palette =c("#980043","#e7298a","#d4b9da")#,"#FFFFFF")
  }
  #geno_matrix[genos == "NA"] <- NA
  rownames(geno_matrix)=infos[[1]]
  colnames(geno_matrix)=infos[[2]]
  rm(infos)
  gc()
  # snp that are missing for all slected samples (they make heatmap.2 bugs while doing hsculst)
  snpallmissing = apply(geno_matrix, 2, function(snp) all(is.na(snp)))
  
  # if a sample is NA for all the snp 
  sampleallmissing = apply(geno_matrix, 1, function(snp) all(is.na(snp)))

  titre = paste0("Contig: ", chr, " window: ", (input$window), " snps: ", ncol(geno_matrix), " 1st: ", firstSnpPos, " last: ", lastSnpPos, " samples with 0 snps: ", length(which(sampleallmissing)), "  snp in 0 samples:", length(which(snpallmissing)) )

  popDeco = popDeco()

  genoWindow = list(geno_matrix=geno_matrix, pop_color1=popDeco()$pop_color1, my_palette=c(my_palette,"#FFFFFF"), keylabels=keylabels, popc=popDeco()$popc, titre=titre, pop_code1=popDeco()$pop_code1)

})

output$genoPlot <- renderPlot({
if (is.null(genoWindow()) ) return(NULL)

ddc = apply(genoWindow()$geno_matrix, 2, function(snp) sum(! is.na(snp)))
names(ddc)=""

ddr = apply(genoWindow()$geno_matrix, 1, function(snp) sum(! is.na(snp)))
names(ddr)=""

column_ha = HeatmapAnnotation(nonmissingsamples = anno_barplot(ddc) )
row_ha    = rowAnnotation( nonmissingsnps = anno_barplot(ddr))
pops_ha   = rowAnnotation(pops= genoWindow()$pop_code1, col=list(pops=genoWindow()$popc))

Heatmap(genoWindow()$geno_matrix, name = "Genotypes", row_names_side = "left", top_annotation = column_ha, right_annotation = row_ha, left_annotation =pops_ha ,show_column_dend = FALSE, show_row_dend = FALSE, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names=FALSE, column_title=genoWindow()$titre, na_col = "white", col=structure(  genoWindow()$my_palette,  names=c(2,1,0,NA)), row_names_gp = gpar(fontsize = 8) )#genoWindow()$keylabels))

})

output$SnpFstPlot <- renderPlot({
if (is.null(genoWindow()) ) return(NULL)

  getFST_diploids_fromCodage = function(popnames, SNPDataColumn){  
    # adapted from OutFLANK for bi-allelic snp only coded as 0 ref/ref, 1 ref/alt, 2 alt/alt
    #remove missing data for this locus
    
    popNameTemp=popnames[!is.na(SNPDataColumn)]
    snpDataTemp=SNPDataColumn[!is.na(SNPDataColumn)]
    
    popGenoCounts <- table(popNameTemp,snpDataTemp)
    popGenoCounts[is.na(popGenoCounts)] = 0
    
    #fixed allele or not biallelic
    nb_alleles = dim(popGenoCounts)[2]
    if(nb_alleles!=2){
      return (list(He=NA,FST=NA, MeanPi=NA, MeanDxy=NA))
    }

    pops_sizes = rowSums(popGenoCounts)
    nbar = mean(pops_sizes)
    r = nrow(popGenoCounts) 
    #W&C Fst
    n_c = (r*nbar - sum(pops_sizes^2)/(r*nbar))/(r-1)
    if ("0" %in% colnames(popGenoCounts)) p0=popGenoCounts[,"0"] else p0=rep(0,r)
    if ("1" %in% colnames(popGenoCounts)) p1=popGenoCounts[,"1"] else p1=rep(0,r)
    if ("2" %in% colnames(popGenoCounts)) p2=popGenoCounts[,"2"] else p2=rep(0,r)

    p_freqs = (p0 + p1/2) /pops_sizes
    pbar = sum(pops_sizes*p_freqs)/(nbar*r)
    
    s2 = sum(pops_sizes*(p_freqs - pbar)^2)/((r-1)*nbar)
    if(s2==0){return(0); break}  
    
    h_freqs = p1/pops_sizes
    hbar = sum(pops_sizes*h_freqs)/(nbar*r)
    
    a <- nbar/n_c*(s2 - 1/(nbar-1)*(pbar*(1-pbar)-((r-1)/r)*s2-(1/4)*hbar))
    b <- nbar/(nbar-1)*(pbar*(1-pbar) - (r-1)/r*s2 - (2*nbar - 1)/(4*nbar)*hbar)
    c <- hbar/2
    FST <- a/(a+b+c) 
    
    He <- 1-sum(pbar^2, (1-pbar)^2)
    
    #http://seqanswers.com/forums/showthread.php?t=33639
    #samples are diploid
    al1_counts = (2*p0 + p1) 
    al2_counts = (2*p2 + p1) 
    
    #Pi adapted from popgenome calc_nuc_diversity_within in calc_diversities.R
    N = al1_counts + al2_counts
    intrapops_pi = (al1_counts * al2_counts)/(N*(N-1)/2)

    #Dxy :  adapted from popgenome calc_average_nuc_diversity_between_per_site in site_FST.R for bi-allelic snp only
    Dxy=numeric(r)
    i=0
    for (pop1 in 1:(r-1))
    {
      for (pop2 in (pop1+1):r)
      {
        i = i +1  
        denom=N[pop1]*N[pop2]
        Dxy[i] = (al1_counts[pop1] * al2_counts[pop2]  + al2_counts[pop1] * al1_counts[pop2])/denom
      }
    }

    return(list(He=He, FST=FST, MeanPi=mean(intrapops_pi), MeanDxy=mean(Dxy)))
}



ttt=apply(genoWindow()$geno_matrix, 2, function(snp){getFST_diploids_fromCodage(genoWindow()$pop_code1,snp)} )
He_Fst = data.frame(matrix(unlist(ttt), nrow=length(ttt), byrow=T))
He_Fst = cbind(He_Fst,pos=colnames(genoWindow()$geno_matrix))
colnames(He_Fst) = c("He","FST","MeanPi","MeanDxy","Pos")
my_threshold <- quantile(He_Fst$FST, 0.975, na.rm = T)
# make an outlier column in the data.frame
He_Fst <- He_Fst %>% mutate(outlier = ifelse(He_Fst$FST > my_threshold, "outlier", "background"))
fstplot = ggplot(He_Fst, aes(Pos, FST, colour = outlier)) + geom_point() + scale_color_manual(values=c("black", "red", "grey"))+ labs(title = genoWindow()$titre,
        subtitle = "Weir&Cockerham 84 Fst for bi-allelic snp. (putative outliers are set to be Fst values > the 97.5% quantile of the current window Fsts; see Whitlock and Lotterhos 2015 for an approriate method)")
Piplot = ggplot(He_Fst, aes(Pos, MeanPi, colour = outlier)) + geom_point() + ggtitle("Mean Pi for bi-allelic snp.") + scale_color_manual(values=c("black", "red", "grey"))
Dxyplot = ggplot(He_Fst, aes(Pos, MeanDxy, colour = outlier)) + geom_point() + ggtitle("Mean Dxy for bi-allelic snp.") + scale_color_manual(values=c("black", "red", "grey"))
heplot = ggplot(He_Fst, aes(Pos, He, colour = outlier)) + geom_point() + ggtitle("Expeced H for bi-allelic snp.") + scale_color_manual(values=c("black", "red", "grey"))
hefstplot  = ggplot(He_Fst, aes(He, FST, colour = outlier)) + geom_point() + ggtitle("Fst vs.He for bi-allelic snp.") + scale_color_manual(values=c("black", "red", "grey"))
(fstplot/Piplot/Dxyplot/heplot/hefstplot)
})

# See also Plink : plink --file data --cluster-missing 
# Systematic batch effects that induce missingness in parts of the sample will induce correlation between
# the patterns of missing data that different individuals display. One approach to detecting correlation 
# in these patterns, that might possibly idenity such biases, is to cluster individuals based on their identity-by-missingness (IBM). 
# This approach use exactly the same procedure as the IBS clustering for population stratification, except the distance between 
# two individuals is based not on which (non-missing) allele they have at each site, but rather the proportion of sites for which 
# two individuals are both missing the same genotype. 
# see also the grur package 
# one can also run a ape::pcoa on the pairwise missingness distances and use ape::biplot
output$missingPlot <- renderPlot({
  if (is.null(genoWindow()) ) return(NULL)

  Pairemissingness <- function(i, j){
    #fraction of SNPs with missing data for one of two sample genotypes divided by the total number of SNPs genotyped for at least one of the samples
    #res = mapply ( nbNA, geno_matrix[i,], geno_matrix[j,])
    # 0 if geno_matrix[i,z] and geno_matrix[j,z] != NA; 1 if only one param is NA and 2 if both params are NA 
    res = is.na(genoWindow()$geno_matrix[i,]) + is.na(genoWindow()$geno_matrix[j,])
    nonmissings = sum(res == 0)
    onlyonemissing  = sum(res == 1)
    return(onlyonemissing / (onlyonemissing + nonmissings))
  }
  
  # get pairwise combination of samples index
  combinaisons = combn(nrow(genoWindow()$geno_matrix), 2)

  missingness = mapply(Pairemissingness, combinaisons[1,], combinaisons[2,])
  if (sum(missingness) == 0) return(NULL)
  #here values are i.e for 4 samples [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4) ]
  #to put the vector in a upper diag matrix
  #a = 1:6; b= matrix(0, 4, 4); b[lower.tri(b, diag=FALSE)] <- a; b <- t(b)
  pairwisemissingness= matrix(0,nrow(genoWindow()$geno_matrix), nrow(genoWindow()$geno_matrix))
  pairwisemissingness[lower.tri(pairwisemissingness, diag=FALSE)] <- missingness
  pairwisemissingness = t(pairwisemissingness)
  #lower diag matrix
  pairwisemissingness[lower.tri(pairwisemissingness, diag=FALSE)] <- missingness
  
  ibsDist = c(1 - ibs$ibs[lower.tri(ibs$ibs, diag = F)])
  pearson <- cor.test(missingness, ibsDist, method="pearson")

  #reorder based on ibsDendro
  pairwisemissingness <- pairwisemissingness[ibsDendro$order, ibsDendro$order]
  
  rownames(pairwisemissingness)=rownames(genoWindow()$geno_matrix)[ibsDendro$order]
  colnames(pairwisemissingness)=rownames(pairwisemissingness)


  pops_ha = rowAnnotation(pops= genoWindow()$pop_code1, col=list(pops=genoWindow()$popc))

  Heatmap(pairwisemissingness, name = "Missingness", cluster_rows=ibsDendro, cluster_columns=ibsDendro, left_annotation =pops_ha ,column_title=paste0(genoWindow()$titre, " pearson:", round(pearson$estimate,digits=4), " p=",signif(pearson$p.value, digits=4) ), col= genoWindow()$my_palette, column_names_side="top", row_names_side = "left")

})


output$tsnePlot <- renderPlot({
 if (is.null(tsneRun()) ) return(NULL)
 df = cbind(as.data.frame(tsneRun()$tsneRes$Y), pops=tsneRun()$pops, samples=tsneRun()$samples)
 colnames(df) = c('axe1','axe2', 'pops', 'samples')

 nbpops=length(unique(tsneRun()$pops) )
 #if (nbpops == 1) nbpops = 2

#colorier les points avec info pop et entourer les clusters

  if(nbpops > 1)
      plot_p = ggplot(df, aes_string(x="axe1", y="axe2", color="pops")) + ggtitle(paste0("maxLD=",maxLD," maxMissing=", globalmaxMissing," minMaf=", globalminMAF)) +  geom_point(size=3) +  xlab("") + ylab("") + ggtitle("") +
      theme(axis.text.x=element_blank(), axis.text.y=element_blank(), legend.direction = "horizontal", legend.position  = "bottom", legend.box = "horizontal", panel.background = element_blank()) + scale_colour_manual(values=tsneRun()$popc)
  else
      plot_p=ggplot(df, aes_string(x="axe1", y="axe2")) +  geom_point(size=3) + ggtitle(paste0("maxLD=",maxLD," maxMissing=", globalmaxMissing," minMaf=", globalminMAF))
  if (input$showLabels)   plot_p = plot_p + geom_text(aes(label=samples),hjust=0,vjust=0.2)

  plot_p


})


pcaRun <- eventReactive(input$runpca,{
    if (is.null(input$vcffile))
    {
        fics = parseFilePaths(c(Data="/Data",Results="/Results"),input$servervcffile)
        if (nrow(fics)>0) vcf.fn = fics$datapath[1]
        else return(dummy)
    }
    else
    {
      vcf.fn = input$vcffile$datapath
    }

    if (! file.exists(gds_file))
    {
      seqVCF2GDS(filteredFile,gds_file,  verbose=FALSE,  storage.option="LZ4_RA", parallel = 8)
    }

    genofile <- seqOpen(gds_file)
      
    PCA1 <- snpgdsPCA(genofile,  snp.id=NULL, maf=NaN, missing.rate=1, num.thread=8, verbose=FALSE, autosome.only=FALSE)
    sample.id1 <- PCA1$sample.id
    seqClose(genofile)
            
    popDeco = popDeco()
    
    # make a data.frame
    df =  data.frame(sample.id = sample.id1,
                    pop = factor(popDeco()$pop_code1),
                    Axe1 = PCA1$eigenvect[,1],    # the first eigenvector
                    Axe2 = PCA1$eigenvect[,2],    # the second eigenvector
                    Axe3 = PCA1$eigenvect[,3],    # the 3 eigenvector
                    Axe4 = PCA1$eigenvect[,4],    # the 4 eigenvector
                    Axe5 = PCA1$eigenvect[,5],    # the 5 eigenvector
                    Axe6 = PCA1$eigenvect[,6],    # the 6 eigenvector
                    Axe7 = PCA1$eigenvect[,7],    # the 7 eigenvector                                                            
                    Axe8 = PCA1$eigenvect[,8],    # the 8 eigenvector
                    Axe9 = PCA1$eigenvect[,9],    # the 9 eigenvector                                        
                    Axe10 = PCA1$eigenvect[,10],    # the 10 eigenvector                    
                    couleur = popDeco()$pop_color1,
                    stringsAsFactors = FALSE)
    rownames(df) = sample.id1

    pcaresu <<-list(tableau = df ,
              valpropre = PCA1$eigenval,
              popc = popDeco()$popc,
              minMAF = globalminMAF,
              maxMissing = globalmaxMissing,
              maxLD = maxLD
            )          
})

shared<- reactive ({
 if  (is.null(pcaRun())) return(dummy$tableau)
 df = pcaRun()$tableau[,c(1,2,as.numeric(input$xcol)+2, as.numeric(input$ycol)+2 , 13)]
 firstAxis = paste0("Axe",as.numeric(input$xcol))
 secondAxis = paste0("Axe",as.numeric(input$ycol))
 colnames(df) = c("sample.id","pop",firstAxis,secondAxis,"couleur")
 return(df)
})
shared_dd <-SharedData$new(shared)


output$PCAdownloadData <- downloadHandler(
    filename = function() {
        paste0("filtered-maxLD", maxLD,"-missing",globalmaxMissing,"-maf",globalminMAF,"-PCA-EigenVectors.txt")
    },
    content = function(file) {
      write.table(pcaRun()$tableau, file , append = FALSE, quote = F, sep = "\t",row.names = FALSE, col.names=TRUE)

    }
  )

output$Eigen <- renderPlot({
  if (is.null(pcaRun()) ) return(NULL)
  barplot(pcaRun()$valpropre[1:20], main=paste0("PCA top 20 eigenvalues maxLD=",maxLD," minMAF=", globalminMAF," maxMissing=", globalmaxMissing), col=heat.colors(25),names.arg=paste("axe",seq(1:length(pcaRun()$valpropre[1:20])), sep=""))
})

### Crosstalk entre d3scatter et ggplot ###########
# voir possibilité entre d3scatter et DT https://stackoverflow.com/questions/48238055/selecting-rows-from-a-dt-table-using-crosstalk-in-shiny
# voir éagelement possibilité avec rgl ?rgl::rglShared et https://bwlewis.github.io/rthreejs/crosstalk.html

output$Scatterplot <-  renderD3scatter ({
  if (is.null(shared_dd) ) return(NULL)  
  firstAxis = paste0("Axe",as.numeric(input$xcol))
  secondAxis = paste0("Axe",as.numeric(input$ycol))
  d3scatter(shared_dd, as.formula(paste0("~",firstAxis)), as.formula(paste0("~", secondAxis)), color = ~pop , width = "100%")
})

df <- debounce(reactive(shared_dd$data(withSelection = TRUE)), 250)

output$by_pop <- renderPlot({
    if (is.null(df) ) return(NULL) 
    ggplot(df(), aes(pop, fill = crosstalk::selection_factor(selected_))) +
      geom_bar(stat = "count") +
      ggtitle(paste0("Populations membership of selected individuals: maxLD=",maxLD," minMAF=", globalminMAF," maxMissing=", globalmaxMissing) )+
      crosstalk::scale_fill_selection("#444444", "#9999FF")
  })


output$PcaHclust <- renderPlot({
 if (is.null(pcaRun()) ) return(NULL)
 df = pcaRun()$tableau[,c(as.numeric(input$xcolhc)+2, as.numeric(input$ycolhc)+2 ,1,2, 13)]
 colnames(df) = c("axe1","axe2","samples","pops","couleur")

 fit_cluster_hierarchical=hclust(dist(scale(df[,c(1,2)])))
 nbpops=length(unique(df$pops) )
 
 h_clusters = factor(cutree(fit_cluster_hierarchical, k=input$pcahclustK))
 df = cbind(df, clusters = h_clusters)
 
#colorier les points avec info pop et entourer les clusters
  if(nbpops > 1)
      plot_p=plot_cluster(df, "pops", pcaRun()$popc,  "clusters" )
  else
      plot_p=ggplot(df, aes_string(x="axe1", y="axe2", shape="clusters")) +  geom_point(size=3)
  if (input$pcashowLabels)   plot_p = plot_p + geom_text(aes(label=samples),hjust=0,vjust=0.2)
  if (input$pcashowClusters) plot_p = plot_p + geom_encircle(aes_string(group = "clusters"), s_shape = 1, expand = 0,alpha = 0.2, color = "red", show.legend = FALSE) 
  plot_p + ggtitle(paste0("maxLD=",maxLD," minMAF=", globalminMAF," maxMissing=", globalmaxMissing))

})

output$Scatterplot3d <- renderScatterplotThree({
#output$Scatterplot3d <-renderPlotly({
   
  if (is.null(pcaRun()) ) return(NULL)
  
  x <- pcaRun()$tableau[,as.numeric(input$xcol3d)+2]
  y <- pcaRun()$tableau[,as.numeric(input$ycol3d)+2]
  z <- pcaRun()$tableau[,as.numeric(input$zcol3d)+2]
  axisLabels = paste0('Axis',c(input$xcol3d, input$ycol3d, input$zcol3d) )
  colors =  unique(pcaRun()$tableau$couleur)
  scatterplot3 <<- scatterplot3js(x,y,z, height = '1000px', width='1000px', axisLabels = axisLabels, size=0.5, labels = pcaRun()$tableau[,1], color=pcaRun()$tableau$couleur)
  #plot_ly(pcaRun()$tableau, x = ~Axe1, y = ~Axe2, z=~Axe3, color=~pop, colors=colors, hoverinfo = 'text', text = ~paste(pop, " : ", sample.id)) %>%
  #layout(scene = list(xaxis = list(title = 'Axe1'), yaxis = list(title = 'Axe2'), zaxis = list(title = 'Axe3')))

})

#output$Scatterplotpairs <- renderPairsD3({
output$Scatterplotpairs <- renderPlot({
    if (is.null(pcaRun()) ) return(NULL)
    
    labels = paste("axe", 1:input$nbaxes, sep='')
    labels = paste(labels,format(round(pcaRun()$valpropre[1:input$nbaxes], 2) ), "%")
    pairs(pcaRun()$tableau[,3:(2+input$nbaxes)],lower.panel=NULL, col = pcaRun()$tableau$couleur,  labels = labels, main=paste0("maxLD=",maxLD," minMAF=", globalminMAF," maxMissing=", globalmaxMissing))
    #pairsD3(pcaRun()$tableau[,c( 3,4,5,6 )], group = pcaRun()$tableau$pop, big=TRUE)
})

output$downScatter3d <- downloadHandler(
    filename = function() {
      paste0(basename(pcafilteredFile),"PCA.html")
    },
    content = function(file) {
        #write.table(df, file , append = FALSE, quote = F, sep = "\t",row.names = TRUE, col.names=TRUE)
        savePairs(scatterplot3, file, selfcontained = T)
    }
  )

ibsRun <- eventReactive(input$runibs,{
    if (is.null(input$vcffile))
    {
        fics = parseFilePaths(c(Data="/Data",Results="/Results"),input$servervcffile)
        if (nrow(fics)>0) vcf.fn = fics$datapath[1]
        else return(dummy)
    }
    else
    {
      vcf.fn = input$vcffile$datapath
    }
    
    if (! file.exists(gds_file))
    {
      snpgdsVCF2GDS(filteredFile,gds_file, method="biallelic.only", verbose=FALSE)
    }

    genofile <- seqOpen(gds_file)
    
    ibs  <- snpgdsIBS(genofile,  snp.id=NULL,  num.thread=8, verbose=TRUE, autosome.only=FALSE)
    ibs.hc <- snpgdsHCluster(ibs)

    colnames(ibs$ibs) = ibs$sample.id
    rownames(ibs$ibs) = ibs$sample.id
    seqClose(genofile)
   
    popDeco = popDeco()
    ibsresu=list(
              ibs = ibs ,
              ibs.hc = ibs.hc,
              popc = popDeco()$popc,
              couleur = popDeco()$pop_color1,
              pop_code1 = popDeco()$pop_code1,
              minMAF = input$ibsminMAF,
              maxMissing = input$ibsmaxMissing
            )          
    
    return(ibsresu) 

})

output$ibs <- renderPlot({
  if (is.null(ibsRun()) ) return(NULL)
  
pops_ha = rowAnnotation(pops= ibsRun()$pop_code1, col=list(pops=ibsRun()$popc))
Heatmap(ibsRun()$ibs$ibs, cluster_rows = ibsRun()$ibs.hc$hclust, cluster_columns = ibsRun()$ibs.hc$hclust, name = "IBS", row_names_side = "left", left_annotation =pops_ha , show_column_names=FALSE, column_title=paste0("maxLD=",maxLD," minMAF=",globalminMAF," maxMissing=",globalmaxMissing), na_col = "white", col=terrain.colors(20), row_names_gp = gpar(fontsize = 8) ,column_dend_height = unit(3, "cm"), row_dend_width = unit(3, "cm"))

 })

output$ibsClusters <- renderPlot({
  if (is.null(ibsRun()) ) return(NULL)
  set.seed(1234)

  # Determine groups of individuals automatically
  rv <- snpgdsCutTree(ibsRun()$ibs.hc)
  op <- par(cex=0.7)
  plot(rv$dendrogram, leaflab="perpendicular",  main="Hclustering and automatic grouping")
  par(op)
})

output$distruct <- renderImage({
   if (is.null(plotfastStructpopHelper() )) 
        return(NULL)
   
   width  <- session$clientData$output_distruct_width
   height <- session$clientData$output_distruct_height
   
   # Return a list containing the filename
   list( src = normalizePath(plotfastStructpopHelper() ),
         contentType = 'image/png+xml',
         width = '100%',
         height = '100%',
         alt = paste0("Admixture proportions inferred by fastStructure K=", isolate(input$K)) )
}, deleteFile = FALSE)       


fstRun <- eventReactive(input$runfst,{
    if (is.null(input$vcffile))
    {
        fics = parseFilePaths(c(Data="/Data",Results="/Results"),input$servervcffile)
        if (nrow(fics)>0) vcf.fn = fics$datapath[1]
        else return(dummy)
    }
    else
    {
      vcf.fn = input$vcffile$datapath
    }

    popDeco = popDeco()
    fics = parseFilePaths(c(Data="/Data",Results="/Results"),input$serverpopMapfile)
    if (nrow(fics)==0 | length(popDeco$popc)==1) {
          fstpairs = c()
          return(NULL)
    } else {
    
    if (! file.exists(gds_file))
    {
      snpgdsVCF2GDS(filteredFile,gds_file, method="biallelic.only", verbose=FALSE)
    }
    genofile <- seqOpen(gds_file)

    nbpops = length(popDeco$popc)
    fstpairs = matrix(0,nrow=nbpops,ncol=nbpops,dimnames= list(names(popDeco$popc), names(popDeco$popc) ) )
    
    for (i in 1:(nbpops-1) )
    {
      for (j in (i+1):nbpops )
      {
         pop1 = names(popDeco$popc)[i]
         pop2 = names(popDeco$popc)[j]
         filteredpops = popDeco$pop_code1[popDeco$pop_code1 %in% c(pop1,pop2)]
         indpop1pop2 = names(filteredpops)
         
         fst <- snpgdsFst(genofile, sample.id=indpop1pop2, population=as.factor(filteredpops),  autosome.only=F,  method="W&C84",  remove.monosnp=T)

         fstpairs[i,j] = round(fst$MeanFst, digits=5)
         fstpairs[j,i] = round(fst$Fst , digits=5)
      }
    }

    seqClose(genofile)
    }

    return(list(fstpairs=fstpairs, minMAF=input$fstminMAF, maxMissing=input$fstmaxMissing ) ) 

})

output$Fstpairs = DT::renderDataTable({
  if (is.null(fstRun()) ) return(NULL)

   DT::datatable( as.data.frame(fstRun()$fstpairs), caption = paste0("Data filtered : minMAF=",fstRun()$minMAF," maxMissing=",fstRun()$maxMissing), options = list(pageLength = 30, searching = FALSE) )
})

output$treePlot <- renderPlot({
      forme   = "rectangular" 
	    title=paste0("NJ plot on pairewise Fst. Data filtered : maxLD=",maxLD," minMAF=" , globalminMAF," maxMissing=",globalmaxMissing)
	    if (is.null(fstRun()) | (nrow(fstRun()$fstpairs) < 3) ) return(NULL)
      distance = fstRun()$fstpairs
		  tree = nj(distance)
		  p <- ggtree(tree,layout=forme) + geom_text(aes(label=label), size=3, color="purple", hjust=-0.3) + theme_tree2("gray86") +ggtitle(title)
		  plot(p)
      
})


plotfastStructpopHelper <- eventReactive(input$runfastStructure, {
   if (is.null(input$vcffile))
    {
        fics = parseFilePaths(c(Data="/Data",Results="/Results"),input$servervcffile)
        if (nrow(fics)>0) vcf.fn = fics$datapath[1]
        else return(dummy)
    }
    else
    {
      vcf.fn = input$vcffile$datapath
    }
    
    #bed file is now cereated with SeqArray
    #cmd=  paste0("/opt/biotools/bin/plink --vcf ",filteredFile," --make-bed --double-id --allow-extra-chr -autosome-num 950 --out ",filteredFile)
    #cmd=  paste0("/opt/biotools/bin/plink --vcf ",filteredFile," --make-bed --double-id --allow-extra-chr --chr-set 950 --out ",filteredFile)
    #system(cmd)
    
    files <- Sys.glob(paste0(filteredFile ,'*.meanQ') )
    unlink(files)
    print(filteredFile)

    for (K in isolate(input$K[1]):isolate(input$K[2] ) )
    {
      if (K < isolate(input$K[2]))
       cmd = paste0("python /opt/biotools/bin/fastStructure/structure.py -K ",K, " --prior=",isolate(input$Prior)," --cv ", input$fastcv," --input=",filteredFile," --output=",filteredFile,"-fastStr &")
      else 
       cmd = paste0("python /opt/biotools/bin/fastStructure/structure.py -K ",K, " --prior=",isolate(input$Prior)," --cv ", input$fastcv," --input=",filteredFile," --output=",filteredFile,"-fastStr ")

      system(cmd)
    }

  
    cmd = paste0("python /opt/biotools/bin/fastStructure/chooseK.py --input=",filteredFile,"-fastStr ")
    ModelComplexity=system(cmd, intern=TRUE)
    titre = paste0("maxLD=",maxLD," minMAF=",globalminMAF, " maxMissing=",globalmaxMissing," ", paste0(ModelComplexity, collapse=" ") )

    files <- Sys.glob(paste0(filteredFile,'-fastStr*.meanQ') )
    qlist <- readQ(files)
        
    kelly <- c("#BE0032","#F3C300","#875692","#F38400","#A1CAF1","#C2B280","#848482","#008856","#E68FAC","#0067A5","#F99379","#604E97","#F6A600","#B3446C","#DCD300","#882D17","#8DB600", "#654522","#E25822","#2B3D26")

    ttt = system(paste0('vcf-query -l ', filteredFile), intern=T)
    qlist <- lapply(qlist,"rownames<-", ttt)
    haut = round(40 / length(qlist))
    nopops =T

# TODO check if we have more than one plot nb K to test > 1

    popDeco = popDeco()
    fics = parseFilePaths(c(Data="/Data",Results="/Results"),input$serverpopMapfile)
    if (nrow(fics)>0) {
      popGroups = as.data.frame(popDeco$pop_code1)  
      if (length(popDeco$popc) > 1) {
        nopops = F 
        colnames(popGroups) = 'Pops'

        p = plotQ(qlist,returnplot=T, exportplot=F,basesize=11, grplab=popGroups,grplabsize=3,linesize=0.8,pointsize=3, showindlab=T,useindlab=T, ordergrp=T, indlabangle=90,indlabvjust=1, height=haut, width=120, indlabsize=3, divsize=2, spbgcol="grey", imgoutput="join", panelspacer=0.4, splab=paste0("K=",sapply(qlist,ncol)),splabsize=20 , clustercol=kelly,barbordersize=0.2 , barbordercolour='white', showtitle=T, titlesize=6,titlelab= titre)
        
      } 
      else 
      {
        showModal(modalDialog(title = "Warning", "The list of samples in popMap is different from samples in the VCF. Population info is ignored. !"))
        nopops = T
      }
    }

    if (nrow(fics) == 0 | nopops) 
      p= plotQ(qlist,returnplot=T,exportplot=F,basesize=11,linesize=0.8,pointsize=3,  height=haut, width=140,indlabsize=, divsize=2, spbgcol="grey", imgoutput="join", panelspacer=0.4, splab=paste0("K=",sapply(qlist,ncol)),splabsize=20 , indlabangle=90, indlabvjust=1, showindlab=T,useindlab=T, clustercol=kelly,barbordersize=0.2 , barbordercolour='white', showtitle=T, titlesize=6, titlelab= titre)
      

  png('/tmp/test.png', width = 1800, height = 800)
  gridExtra::grid.arrange(p$plot[[1]])
  dev.off()

output$meanQtable <- DT::renderDataTable({
    fic = paste0(filteredFile,"-fastStr.",input$kfile, ".meanQ")
    
    cmd = paste0("grep  'Marginal Likelihood =' ", filteredFile,"-fastStr.",input$kfile,".log")
    likelihood = system(cmd, intern=T)

    if  ( file.exists(fic) )  {
        df = read.table(fic, header=F)

        ttt = system(paste0('vcf-query -l ', filteredFile), intern=T)
        rownames(df)  = ttt
        
        colnames(df) =  paste0("cluster", 1:input$kfile)
        
        popDeco = popDeco()
        if (length(popDeco$popc) > 1 )  df = cbind(pop=popDeco$pop_code1, df, stringsAsFactors=F)

        DT::datatable(df, caption=likelihood, options= list(pageLength = 30, searching = FALSE))# %>% DT::formatStyle(names(df), backgroundColor =  DT::styleInterval(0.51, c( "white", "lightblue")))
    }
  })
  outfile <- "/tmp/test.png"

  })

observe({
    minK <- input$K[1]
    maxK <- input$K[2]

    x = as.list(minK:maxK)
    names(x) = paste0("k=", minK:maxK)
    
    updateSelectInput(session, "kfile",
      choices  = x,
      selected = paste0("k=", maxK)
    )
  })

 

output$downloadData <- downloadHandler(
    filename = function() {
      paste0(basename(filteredFile),"-fastStr.K.",input$kfile, ".meanQ")
    },
    content = function(file) {
      fileTodownload = paste0(filteredFile,"-fastStr.",input$kfile, ".meanQ")
      
      if (file.exists(fileTodownload)) 
      {
        df = read.table(fileTodownload, header=F)

        ttt = system(paste0('vcf-query -l ', filteredFile), intern=T)
        df = cbind(df, ttt)
        
        colnames(df) =  c("Sample", paste0("cluster", 1:input$kfile) )
        popDeco = popDeco()
        if (length(popDeco$popc) > 1 )  df = cbind( df, pop=popDeco$pop_code1, stringsAsFactors=F)
        write.table(df, file , append = FALSE, quote = F, sep = "\t",row.names = TRUE, col.names=TRUE)
        #file.copy(fileTodownload,file)
      }
    }
  )

# Show modal when set model params button is clicked.
observeEvent(input$setParams, {
       if (is.null(sfsRun())  ) return(NULL)
       showModal(ModelParams(input$Model, models_tooltips[input$Model,2]))
      #  if (input$Model == "IM")  showModal(IMParams())
      #  if (input$Model == "prior_onegrow_mig")  showModal(prior_onegrow_migParams())
      #  if (input$Model == "SI")  showModal(SIParams())
      #  if (input$Model == "SC")  showModal(SCParams())
})

sfsRun <- eventReactive(input$runsfs,{
    if (is.null(input$vcffile))
    {
        fics = parseFilePaths(c(Data="/Data",Results="/Results"),input$servervcffile)
        if (nrow(fics)>0) vcf.fn = fics$datapath[1]
        else return(NULL)
    }
    else
    {
      vcf.fn = input$vcffile$datapath
    }

   popDeco = popDeco()
   if (length(popDeco$popc)>1){
    choixPops = 1:length(popDeco$popc)
    names(choixPops) = names(popDeco$popc)
    
    updateSelectInput(session, "firstPop", label = NULL, choices = choixPops,selected = 1)
    updateSelectInput(session, "secondPop", label = NULL, choices = choixPops,selected = 2)
    updateSelectInput(session, "firstPop1", label = NULL, choices = choixPops,selected = 1)
    updateSelectInput(session, "secondPop1", label = NULL, choices = choixPops,selected = 2)

    fics = parseFilePaths(c(Data="/Data",Results="/Results"),input$serverpopMapfile)
    ficpopmap = fics$datapath[1]

    #la version de moments ne semble pas fonctionner quand il y a des données manquantes ?!
    #dd = momentsMisc$make_data_dict_vcf(filteredFile, ficpopmap, filter=FALSE)
    ddRes = make_data_dict_vcf_Mem(filteredFile, ficpopmap, filter=FALSE, maxMemPercent=isolate(input$maxMemPercent))
    if (py_to_r(ddRes[0])) showModal(modalDialog(title = "Warning", "Mem used >= ",isolate(input$maxMemPercent),". Only the first ", py_to_r(ddRes[1]), "Snps are used"))

    if (py_to_r(ddRes[1])== 0) {
      dd = dadiMisc$make_data_dict_vcf(filteredFile, ficpopmap, filter=FALSE)
      return(dd)
    }else  return(ddRes[2])
   }else{
     showModal(modalDialog(title = "Warning", "You have to provide a population Map file !"))
     return(NULL)
   }
})
    
output$sfsPlot <- renderImage({
  
        if (is.null(sfsRun())  ) return(NULL)
        
        i = as.integer(input$firstPop)
        j = as.integer(input$secondPop)
        if( is.na(i)|| is.na(j)) {return(NULL) }

        popDeco  = popDeco()
        popnames = names(popDeco$popc)
        popsizes = table(popDeco$pop_code1)
        pop_ids  = popnames[c(i,j)]

        ns = popsizes[pop_ids]

        #fs = dadiSpectrum$Spectrum$from_data_dict(sfsRun(), pop_ids, ns, polarized=F)
        fs = momentsSpectrum$Spectrum$from_data_dict(sfsRun(), pop_ids, ns, polarized=input$fold)
        
        #projection sizes
        down1 = min(input$downFirst, ns[1])
        down2 = min(input$downSencond, ns[2])
        # updateNumericInput(session, "downFirst", label = NULL, value = down1,  max = popsizes[pop_ids[1]])
        # updateNumericInput(session, "downSecond", label = NULL, value = down2,  max = popsizes[pop_ids[2]]) 

        fs = fs$project(as.integer(c(down1,down2)) )
        ns = fs$sample_sizes        
        
        plt$figure(figsize=c(8,8)) 
        plt$suptitle(paste0(pop_ids[1],"-", pop_ids[2], " Fst: ", round(py_to_r(fs$Fst()),3)," S: ", round(py_to_r(fs$S()),3)),fontsize=10)
        
        moments$Plotting$plot_single_2d_sfs(fs, vmin=1e-2)
        
        ImgFs = paste0(tempfile(),".png")
        plt$savefig(ImgFs)
        plt$close()
        list(src = ImgFs  ,
         contentType = 'image/png',
         width = "100%",
         height = "100%",
         alt = "2D FS")

}, deleteFile = TRUE)

demoModel <- eventReactive( input$runInference  ,{    
if (is.null(sfsRun())  ){
  showModal(modalDialog(
        title = "Warning", "You have to set parameters for the selected Model !"))
  return(NULL)
}
params = isolate(Modelparams$data)

if (is.null(params) ){showModal(modalDialog(
        title = "Warning", "You have to set parameters for the selected Model !",))
  return(NULL)
} 

if ( isolate(input$Model) !=  params$model){
  showModal(modalDialog(
        title = "Warning", "You have to set parameters for the selected Model !"))
  return(NULL)
}

  i = as.integer(isolate(input$firstPop1))
  j = as.integer(isolate(input$secondPop1))
  if( is.na(i)|| is.na(j)) {return(NULL) }

  popDeco  = popDeco()
  popnames = names(popDeco$popc)
  popsizes = table(popDeco$pop_code1)

  pop_ids = popnames[c(i,j)]
  ns = popsizes[pop_ids]
  #fs = dadiSpectrum$Spectrum$from_data_dict(sfsRun(), pop_ids, ns, polarized=F)
  fs = momentsSpectrum$Spectrum$from_data_dict(sfsRun(), pop_ids, ns, polarized=isolate(input$foldInference))
  
  #projection sizes
  pop1size = min(isolate(input$downFirst1), ns[1])
  pop2size = min(isolate(input$downSencond1), ns[2])
  # updateNumericInput(session, "downFirst1", label = NULL, value = pop1size,  max = popsizes[i])
  # updateNumericInput(session, "downSecond1", label = NULL, value = pop1size,  max = popsizes[j])
  
  fs = fs$project(as.integer(c(pop1size,pop2size)) )
  ns = fs$sample_sizes        
  
  # pts_l = c(40L,50L,60L)
  # if (input$Model == "IM")  func = dadi$Demographics2D$IM
  # if (input$Model == "prior_onegrow_mig")  func = prior_onegrow_mig
  # if (input$Model == "SC")  func = SC
  # if (input$Model == "SI")  func = SI
  testModel = isolate(input$Model)
  # if (testModel == "IM")  func = moments$Demographics2D$IM
  # if (testModel == "prior_onegrow_mig")  func = prior_onegrow_mig
  # if (testModel == "SC")  func = SC
  # if (testModel == "SI")  func = SI
  
  cmd = paste0("func=globals()['",testModel,"']")
  py_run_string(cmd)
  func = py$func
  upper_bound = params$upper_bounds
  lower_bound = params$lower_bounds
  lower_bound[lower_bound==0] = 1e-3
  # This is our initial guess for the parameters, which is somewhat arbitrary.
  p0 = params$p0
  print(p0)
  print(upper_bound)
  print(lower_bound)
  # Make the extrapolating version of our demographic model function.
  #func_ex = dadi$Numerics$make_extrap_log_func(func)
  # Perturb our parameters before optimization. This does so by taking each
  # parameter a up to a factor of two up or down.
  #p0 = dadi$Misc$perturb_params(p0, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)
  p0 = momentsMisc$perturb_params(p0, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)
  logfile = paste0(tempfile(),".log")
  
  verbosity = 1
  if (input$OptiMethod=="optimize_log" )
  {
    # popt = dadi$Inference$optimize_log(p0, fs, func_ex, pts_l, lower_bound=lower_bound,upper_bound=upper_bound,verbose=1, maxiter=params$MaxIter)}
    popt = moments$Inference$optimize_log(p0, fs, r_to_py(func), lower_bound=lower_bound, upper_bound=upper_bound, verbose=verbosity, maxiter=input$MaxIter, output_file=logfile)
  }  
  if (input$OptiMethod=="optimize" )
  {
    popt = moments$Inference$optimize(p0, fs, r_to_py(func), lower_bound=lower_bound, upper_bound=upper_bound, verbose=verbosity, maxiter=input$MaxIter, output_file=logfile)
  }
  if (input$OptiMethod=="optimize_log_lbfgsb" )
  {
    popt = moments$Inference$optimize_log_lbfgsb(p0, fs, r_to_py(func), lower_bound=lower_bound, upper_bound=upper_bound, verbose=verbosity, maxiter=input$MaxIter, output_file=logfile)
  }    
  if ( input$OptiMethod=="dual_anneal")
  {
    ##### Set optimization values
    maxiterGlobal=input$MaxIter #maximum global search iterations - in each iteration, it explores twice the number of parameters
    accept=1 #parameter for acceptance distribution (lower values makes the probability of acceptance smaller)
    visit=1.01 #parameter for visiting distribution (higher values makes the algorithm jump to a more distant region)
    Tini=50 #initial temperature
    no_local_search=FALSE #If set to True, a Generalized Simulated Annealing will be performed with no local search strategy applied
    local_method='L-BFGS-B' #local search method
    maxiterLocal=20 #maximum local search iterations
    popt = myInference$optimize_dual_anneal(p0, fs, r_to_py(func), lower_bound=lower_bound, upper_bound=upper_bound,
    no_local_search=no_local_search, local_method=local_method, local_maxiter=maxiterLocal,maxiter=maxiterGlobal, Tini=Tini, accept=accept, visit=visit, verbose=verbosity, full_output=FALSE, output_file=logfile)
  }

  # Calculate the best-fit model AFS.
  #bestmodel = py_call(func,popt, ns, pts_l)
  bestmodel = py_call(func,popt, ns)
  ll_model = moments$Inference$ll_multinom(bestmodel, fs)
  theta = moments$Inference$optimal_sfs_scaling(bestmodel, fs)

return( list(params = params, pop_ids = pop_ids , testModel = isolate(input$Model) , fs=fs, bestmodel=bestmodel, ll_model = ll_model, theta=theta, popt=popt, func = func, ns=ns, logfile=logfile) ) 

}, ignoreInit = TRUE)

output$demoPlot <- renderImage({
        if (is.null(demoModel()$fs) )return(NULL)
        plt$figure(figsize=c(8,8)) 
       
        moments$Plotting$plot_2d_comp_multinom(demoModel()$bestmodel, demoModel()$fs, vmin=1, pop_ids =c(demoModel()$pop1,demoModel()$pop2), adjust=FALSE, show=FALSE)
        plt$subplots_adjust(top=0.90)
        plt$suptitle(paste0(demoModel()$testModel," Model on ",demoModel()$pop_ids[1],"-", demoModel()$pop_ids[2], " popt:[", toString(round(py_to_r(demoModel()$popt),3)), "] Theta:", round(py_to_r(demoModel()$theta), 2), " logL:", round(py_to_r(demoModel()$ll_model),4) ),fontsize=10) 
        # " Fst:", round(py_to_r(fs$Fst()),3)," S:", round(py_to_r(fs$S()),3),
 
        ImgFs = paste0(tempfile(),".png")
        plt$savefig(ImgFs)
        plt$close()
        list(src = ImgFs  ,
         contentType = 'image/png',
         width = "100%",
         height = "100%",
         alt = "2D FS")

}, deleteFile = TRUE)

output$PlotModel <- renderImage({
  if (is.null(demoModel()$fs) )return(NULL)
  # plotting demographic model
         mu    = input$mutationRate/100
         gtime = input$generationTime
         nref=py_to_r(demoModel()$theta)/(4*mu)
         plot_mod = moments$ModelPlot$generate_model(demoModel()$func, demoModel()$popt, demoModel()$ns)
         moments$ModelPlot$plot_model(plot_mod, pop_labels=demoModel()$pop_ids, nref=nref, draw_scale=FALSE, gen_time=gtime, gen_time_units="KYear", reverse_timeline=TRUE)

        ImgModel = paste0(tempfile(),".png")
        plt$savefig(ImgModel)
        plt$close()
        list(src = ImgModel  ,
         contentType = 'image/png',
         width = "100%",
         height = "100%",
         alt = "Representation of the model")

      }, deleteFile = TRUE)


output$logLikPlot <- renderPlot({
  if (is.null(demoModel()$fs) )return(NULL)
   library("data.table")
  # plotting demographic model
  logL = fread(demoModel()$logfile,select=c(1,2),sep=",")
  logL = cbind(logL, 1:nrow(logL))
  colnames(logL) = c("cpt","logLikelihood", "Step")
  
  ggplot(logL, aes(x=Step, y=logLikelihood)) +  geom_line() 
  })
    
}
)



######################## Refs
# https://indico.io/blog/visualizing-with-t-sne/
# https://www.r-bloggers.com/playing-with-dimensions-from-clustering-pca-t-sne-to-carl-sagan/
# https://distill.pub/2016/misread-tsne/
# L.J.P. van der Maaten and G.E. Hinton. Visualizing High-Dimensional Data Using t-SNE. Journal of Machine Learning Research 9(Nov):2579-2605, 2008.

# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0056883
# https://www.biorxiv.org/content/biorxiv/early/2017/03/08/114884.full.pdf ====>  https://www.ncbi.nlm.nih.gov/pubmed/28718343
# https://www.biostars.org/p/295174/

# https://towardsdatascience.com/dimensionality-reduction-for-data-visualization-pca-vs-tsne-vs-umap-be4aa7b1cb29

# Zheng X, Levine D, Shen J, Gogarten S, Laurie C, Weir B (2012). "A High-performance Computing Toolset for Relatedness and Principal Component Analysis of SNP Data." Bioinformatics, 28(24), 3326-3328. doi: 10.1093/bioinformatics/bts606.

# James R. Whiting, Josephine R. Paris, Mijke J. van der Zee, Paul J. Parsons, Detlef Weigel, Bonnie A. Fraser. https://doi.org/10.1101/2020.10.14.339333 

# Francis, R. M. (2017). POPHELPER: an R package and web app to analyse and visualize population structure. Molecular Ecology Resources, 17(1), 27-32. DOI: 10.1111/1755-0998.12509

# example data : https://datadryad.org/stash/dataset/doi:10.5061/dryad.kp11q

# Whitlock, M. C., and K. E. Lotterhos 2015. Reliable detection of loci responsible for local adaptation: inference of a null model through trimming the distribution of {FST}. Am. Nat. 186:S24-S36. 
  
# Weir, B. S., and C. C. Cockerham 1984. Estimating F-statistics for the analysis of population structure. Evolution 38:1358-1370. 