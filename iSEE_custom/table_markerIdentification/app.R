stopifnot(suppressPackageStartupMessages({
  require(iSEE)
  require(scRNAseq)
  require(Seurat)
  require(scater)
  require(shiny)
  require(shinyjs)
}))

data(allen)

# Example data ----

sce <- as(allen, "SingleCellExperiment")
counts(sce) <- assay(sce, "tophat_counts")
sce <- normalize(sce)

set.seed(1234)
sce <- runPCA(sce, ncomponents=4)
set.seed(1234)
sce <- runTSNE(sce)

rowData(sce)$mean_log <- rowMeans(logcounts(sce))
rowData(sce)$var_log <- apply(logcounts(sce), 1, var)

# Import custom panel + function ----

source("custom.R")

# Import tour steps ----

# tour <- read.delim("tour.txt", sep=";", quote="")

# Configure the app ----

# RedDimPlot1
redDimArgs <- redDimPlotDefaults(sce, 1)
redDimArgs$Type <- 2L
redDimArgs$ColorBy <- "Feature name"
redDimArgs$ColorByFeatNameAssay <- "logcounts"

# Custom_StatTable1
customStatArgs <- customStatTableDefaults(sce, 1)
customStatArgs$Function <- "CUSTOM_PairwiseMarker"
customStatArgs$ColumnSource <- "Reduced dimension plot 1"
customStatArgs$DataBoxOpen <- T

# Feature_AssayPlot1
featAssay <- featAssayPlotDefaults(sce, 1)
featAssay$YAxisFeatName <- 1L
featAssay$XAxis <- "Column data" 
featAssay$XAxisColData <- "Core.Type" 
featAssay$Assay <- "logcounts"


# Setting up links between plots. -----------------------------------------
redDimArgs$ColorByRowTable <- "Custom statistics table 1"
customStatArgs$SelectByPlot <- "Reduced dimension plot 1"
customStatArgs$ColumnSource <- "Reduced dimension plot 1"
featAssay$YAxisRowTable <- "Custom statistics table 1"
featAssay$SelectByPlot <- "Reduced dimension plot 1"
featAssay$SelectEffect <- "Transparent"
featAssay$SelectAlpha <- 0.1
featAssay$YAxisRowTable <- "Custom statistics table 1"

initialPanels <- DataFrame(
  Name=c("Reduced dimension plot 1", "Custom statistics table 1", "Feature assay plot 1"),
  Width=c(4L, 4L, 4L))

useShinyjs()

app <- iSEE(
  se = sce,
  redDimArgs=redDimArgs, featAssayArgs=featAssay, customStatArgs=customStatArgs,
  redDimMax = 1, colDataMax = 0, featAssayMax = 1, rowDataMax = 0, 
  sampAssayMax = 0, rowStatMax = 0, colStatMax = 0, customDataMax = 0, heatMapMax = 0, customStatMax = 1,
  initialPanels=initialPanels,
  customStatFun=list(CUSTOM_PairwiseMarker = CUSTOM_PairwiseMarker, CUSTOM_GlobalMarker = CUSTOM_GlobalMarker, CUSTOM_SelectedCellMarker = CUSTOM_SelectedCellMarker), 
  #tour=tour, 
  appTitle = "Custom table panel: Running FindMarker on the fly")

# launch the app itself ----

app