# Define custom functions -------------------------------------------------
caches_env = new.env()

caches_env$sharedUI = function(id){
  ns <- NS(id)
  tagList(
    selectInput(ns("markerType"), label = "Specificity of markers", choices = c("low (slow)", "medium (fast)", "high (ultra fast)"), selected = "medium (fast)"),
    iSEE::collapseBox(
      ns("advancedSettings"),
      title = "Advanced Settings",
      open = T,
      splitLayout(
        sliderInput(ns("downsample"), label = "Downsample(prop)", min = 0, max = 1, value = 0.5),
        sliderInput(ns("logfc"), label = "Minimum Fold (log)", min = 0, max = 2, value = 0.5, step = 0.01)
      ),
      splitLayout(
        sliderInput(ns("minDiffPct"), label = "Minimum difference in detection(prop)", min = 0, max = 1, value = 0),
        sliderInput(ns("minPct"), label = "Minimum detection rate", min = 0, max = 1, value = 0.3)  
      )
    ),
    fluidRow(
      column(6, align = "left", selectInput(ns("method"), NULL, choices = c("MAST", "wilcox", "t", "roc", "bimod"), selected = "MAST")),
      column(4, align = "left", actionButton(ns("submit"), "Run"))
    )
  )
}

caches_env$findMarkerSetting <-  list(
  "low (slow)" = list(
    "downsample" = 1,
    "logfc" = 0.25,
    "minDiffPct" = 0,
    "minPct" = 0.1
  ),
  "medium (fast)" = list(
    "downsample" = 0.5,
    "logfc" = 0.5,
    "minDiffPct" = 0,
    "minPct" = 0.3
  ),
  "high (ultra fast)" = list(
    "downsample" = 0.2,
    "logfc" = 0.5,
    "minDiffPct" = 0.5,
    "minPct" = 0.3
  )
)

CUSTOM_PairwiseMarker <- list(
  # computation part
  "computation" = function(se, rows, columns, session, input, output, panel_name, funName, caches_env, mincell = 5, latentVars = c("nUMI", "nGene")){
    ns <- NS(paste0(panel_name, "_", funName))
    
    # Calculation -------------------------------------------------------------
    
    # Define return data.frame;
    result <- data.frame(ID=numeric(0))
    if (is.null(rows)){
      rows <-  rownames(se)
    }
    if (is.null(columns)){
      columns <- colnames(se)
    }
    
    
    # Compute the final results -----------------------------------------------
    
    if (is.null(caches_env$seuratObj)){
      if(package_version(x = packageVersion(pkg = "Seurat")) >= package_version(x = "3.0.0")){
        seuratObj <- Seurat::as.Seurat(from = se)
      } else {
        seuratObj <- Seurat::as.seurat(from = se)
      }
      caches_env$seuratObj = seuratObj
    } else {
      seuratObj = caches_env$seuratObj
    }
    
    
    
    # need variable I and variable II to be defined
    if (!any(is.null(input[[ns("compVar")]]), is.null(input[[ns("var1")]]), is.null(input[[ns("var2")]]))) {
      varName <- input[[ns("compVar")]]
      method <- input[[ns("method")]]
      
      cellID1 <- rownames(subset(seuratObj@meta.data, get(varName) %in% input[[ns("var1")]]))
      cellID2 <- rownames(subset(seuratObj@meta.data, get(varName) %in% input[[ns("var2")]]))
      cellID1 <- cellID1[cellID1 %in% columns]
      cellID2 <- cellID2[cellID2 %in% columns]
      mycomp <- paste0(varName, ":", paste(input[[ns("var1")]], collapse = ";"), "_vs_", paste(input[[ns("var2")]], collapse = ";"))
      logfc <- input[[ns("logfc")]]
      minPct <- input[[ns("minPct")]]
      minDiffPct <- input[[ns("minDiffPct")]]
      downSample <- input[[ns("downsample")]]
      downSample = ifelse( downSample < 1, ceiling(downSample * (length(cellID1) + length(cellID2))), Inf) 
      if (length(cellID1) > mincell & length(cellID2) > mincell & !identical(sort(input[[ns("var1")]]), sort(input[[ns("var2")]]))) {
        isolate({
          tmpDF <- tryCatch({
            suppressWarnings(Seurat::FindMarkers(
              seuratObj, ident.1 = cellID1, ident.2 = cellID2,
              test.use = method, logfc.threshold = logfc,
              min.diff.pct = minDiffPct,
              min.pct = minPct, max.cells.per.ident = downSample
              # latent.vars = latentVars, 
            ))
          },
          error = function(x) {
            return(NULL)
          },
          warning = function(w) {}
          )
          
          if (!is.null(tmpDF)) {
            result <- data.frame(
              # "comprison" = mycomp,
              tmpDF,
              row.names = rownames(tmpDF)
            )
          }
        })
      }
    }
    return(result)
  },
  
  # UI + server logic part
  "UI" = function(se, pointsSelected, session, input, output, panel_name, funName, FUN, caches_env, UIonly = T) {
    ns <- NS(paste0(panel_name, "_", funName))
    possibleLevels <- reactiveValues(var = NULL)
    
    dyUI <- tagList(
      fluidRow(
        column(8, align = "left",
               tags$b(paste0(length(pointsSelected$columns), " cells have been selected"))
        )
      ),
      selectInput(ns("compVar"), "Define group from:", choices = names(colData(se))),
      selectizeInput(ns("var1"), "Group I", choices = c("---"), selected = "---", multiple = TRUE),
      selectizeInput(ns("var2"), "Group II", choices = c("---"), selected = "---", multiple = TRUE),
      caches_env$sharedUI(paste0(panel_name, "_", funName))
    )
    
    if (!UIonly) {
      observe({
        markType = ifelse(is.null(input[[ns("markerType")]]), "medium (fast)", input[[ns("markerType")]])
        selectedSetting <- caches_env$findMarkerSetting[[markType]]
        for(myid in names(selectedSetting)){
          updateSliderInput(session, ns(myid), value = selectedSetting[[myid]])
        }
      })
      
      observe({
        if (is.null(input[[ns("compVar")]])) {
          possibleLevels$var <- "---"
        } else {
          if(is.null(pointsSelected$rows)){
            rows <- rownames(se)
          } else {
            rows <- pointsSelected$rows
          }
          if(is.null(pointsSelected$columns)){
            columns <- colnames(se)
          } else {
            columns <- pointsSelected$columns
          }
          possibleLevels$var <- colData(se[rows, columns])[, unique(input[[ns("compVar")]])]
          updateSelectizeInput(session, ns("var1"), choices = possibleLevels$var, selected = possibleLevels$var[1])
          updateSelectizeInput(session, ns("var2"), choices = possibleLevels$var, selected = possibleLevels$var[1])
          shinyjs::hide(paste0(panel_name, "_download"))
        }
      })
      
      observeEvent(input[[ns("submit")]], {
        shinyjs::disable(ns("submit"))
        shinyjs::hide(paste0(panel_name, "_download"))
        updateActionButton(session, ns("submit"), "Running")
        
        tmpDF <- withProgress(message = "Finding marker between two groups...", value = 0, {
          do.call(FUN, c(list(se, pointsSelected$rows, pointsSelected$columns, session, input, output, panel_name, funName, caches_env)))
        })
        
        shinyjs::enable(ns("submit"))
        updateActionButton(session, ns("submit"), "Run")
        
        mycomp <- paste0(input[[ns("compVar")]], ":", input[[ns("var1")]], "_vs_", input[[ns("var2")]])
        
        output[[panel_name]] <- DT::renderDataTable({
          validate(need(dim(tmpDF)[1] > 0, "No marker identified"))
          
          # update download button
          outDF = tmpDF
          shinyjs::show(paste0(panel_name, "_download"))  
          output[[paste0(panel_name, "_download")]] = downloadHandler(
            filename = paste0(paste0(panel_name, "_download"), ".txt"),
            content = function(f){
              write.table(outDF, f, col.names = T, row.names = F, sep = "\t", quote = F);
            }
          )
          
          # the output should be a data.frame with identical # of rows of original rowData
          tmpDF$ID <- rownames(tmpDF)
          origDF <- data.frame("ID" = rownames(se), "SELECTED" = F, "ORDER" = 1:dim(se)[1], row.names = rownames(se))
          mergeDF <- merge(origDF, tmpDF, by = "ID", all.x = T)
          rownames(mergeDF) <- mergeDF$ID
          mergeDF <- mergeDF[order(mergeDF$ORDER), ]
          mergeDF$SELECTED <- ifelse(is.na(mergeDF[, 4]), F, T)
          
          search_col <- rep(list(NULL), ncol(mergeDF))
          search_col[[2]] <- list(search = "[\"true\"]")
          
          DT::datatable(
            mergeDF,
            filter = "top", rownames = TRUE, caption = mycomp,
            selection = list(mode = "single"),
            options = list(
              search = list(smart = FALSE, regex = TRUE, caseInsensitive = FALSE),
              searchCols = c(list(NULL), search_col),
              columnDefs = list(list(visible = FALSE, targets = c(1, 2, 3))), # first column ID, second SELECTED and third column ORDER are invisible
              pageLength = 5,
              stateSave = F,
              scrollX = TRUE
            )
          )
        })
      })
    }
    return(dyUI);
  },
  
  "caches" = caches_env
)

CUSTOM_SelectedCellMarker <- list(
  # computation part
  "computation" = function(se, rows, columns, session, input, output, panel_name, funName, caches_env, mincell = 5, latentVars = c("nUMI", "nGene")){
    ns <- NS(paste0(panel_name, "_", funName))
    
    # Calculation -------------------------------------------------------------
    
    result = result = data.frame(ID=numeric(0));
    if (is.null(rows)){
      rows = rownames(se);
    }
    if (is.null(columns)){
      columns = colnames(se);
    }
    
    if (is.null(caches_env$seuratObj)){
      if(package_version(x = packageVersion(pkg = "Seurat")) >= package_version(x = "3.0.0")){
        seuratObj <- Seurat::as.Seurat(from = se)
      } else {
        seuratObj <- Seurat::as.seurat(from = se)
      }
      caches_env$seuratObj = seuratObj
    } else {
      seuratObj = caches_env$seuratObj
    }
    
    
    if(!any(is.null(columns))){
      allCells = rownames(seuratObj@meta.data)
      cellID1 = allCells[allCells %in% columns]
      cellID2 = allCells[!allCells %in% cellID1]
      mycomp = paste0("Markers for ", length(cellID1), " selected cells");
      logfc <- input[[ns("logfc")]]
      method <- input[[ns("method")]]
      minPct <- input[[ns("minPct")]]
      minDiffPct <- input[[ns("minDiffPct")]]
      downSample <- input[[ns("downsample")]]
      downSample = ifelse( downSample < 1, ceiling(downSample * (length(cellID1) + length(cellID2))), Inf) 
      if(length(cellID1) > mincell & length(cellID2) > mincell){
        # oldAssay = DefaultAssay(scDataClustered);
        # DefaultAssay(scDataClustered) = "RNA"
        
        isolate({
          tmpDF = tryCatch({
            suppressWarnings(Seurat::FindMarkers(
              seuratObj, ident.1 = cellID1, ident.2 = cellID2,
              test.use = method, logfc.threshold = logfc,
              min.diff.pct = minDiffPct,
              min.pct = minPct, max.cells.per.ident = downSample
              # latent.vars = latentVars, 
            ))}, 
            error = function(x){return(NULL)}
          );
          if (!is.null(tmpDF)) {
            result = data.frame(
              # "comprison" = mycomp,
              tmpDF,
              row.names = rownames(tmpDF)
            )
          }  
        })
      }
    }
    return(result)
  },
  
  # UI + server logic part
  "UI" = function(se, pointsSelected, session, input, output, panel_name, funName, FUN, caches_env, UIonly = T) {
    ns <- NS(paste0(panel_name, "_", funName))
    possibleLevels = reactiveValues( var = NULL);
    
    # update function specific UI ---------------------------------------------
    dyUI = tagList(
      fluidRow(
        column(8, align = "left",
               tags$b(paste0(length(pointsSelected$columns), " cells have been selected"))
        )
      ),
      caches_env$sharedUI(paste0(panel_name, "_", funName))
    )
    
    
    
    if(!UIonly){
      observe({
        markType = ifelse(is.null(input[[ns("markerType")]]), "medium (fast)", input[[ns("markerType")]])
        selectedSetting <- caches_env$findMarkerSetting[[markType]]
        for(myid in names(selectedSetting)){
          updateSliderInput(session, ns(myid), value = selectedSetting[[myid]])
        }
      })
      
      observe({
        if(is.null(input[[ns("compVar")]])){
          shinyjs::hide(paste0(panel_name, "_download"))
        }
      })
      
      observeEvent(input[[ns("submit")]], {
        # browser();
        shinyjs::disable(ns("submit"))
        updateActionButton(session, ns("submit"), "Running")
        shinyjs::hide(paste0(panel_name, "_download"))
        
        
        tmpDF = withProgress(message="Finding marker for selected cells...", value=0,{
          do.call(FUN, c(list(se, pointsSelected$rows, pointsSelected$columns, session, input, output, panel_name, funName, caches_env)));
        })
        
        shinyjs::enable(ns("submit"))
        updateActionButton(session, ns("submit"), "Run")
        
        mycomp = paste0("Markers for ", length(pointsSelected$columns), " selected cells");
        
        output[[panel_name]] = DT::renderDataTable({
          validate(need(dim(tmpDF)[1]>0, "No marker identified"));
          
          # update download button
          outDF = tmpDF
          shinyjs::show(paste0(panel_name, "_download"))  
          output[[paste0(panel_name, "_download")]] = downloadHandler(
            filename = paste0(paste0(panel_name, "_download"), ".txt"),
            content = function(f){
              write.table(outDF, f, col.names = T, row.names = F, sep = "\t", quote = F);
            }
          )
          
          tmpDF$ID = rownames(tmpDF);
          origDF = data.frame("ID" = rownames(se), "SELECTED" = F, "ORDER" = 1:dim(se)[1], row.names = rownames(se));
          mergeDF = merge(origDF, tmpDF, by = "ID", all.x = T);
          rownames(mergeDF) = mergeDF$ID;
          mergeDF = mergeDF[order(mergeDF$ORDER), ];
          mergeDF$SELECTED = ifelse(is.na(mergeDF[, 4]), F, T);
          search_col = rep(list(NULL), ncol(mergeDF));
          search_col[[2]] = list(search = "[\"true\"]")
          
          DT::datatable(
            mergeDF, filter="top", rownames=TRUE, caption = mycomp,
            selection = list(mode="single"),
            options=list(
              search=list(smart=FALSE, regex=TRUE, caseInsensitive=FALSE),
              searchCols=c(list(NULL), search_col),
              columnDefs = list(list(visible=FALSE, targets=c(1,2,3))), # first column ID, second SELECTED and third column ORDER are invisible
              pageLength = 5, 
              stateSave = F,
              scrollX=TRUE
            )
          )
        })
      })  
    }
    return(dyUI)
  },
  
  "caches" = caches_env
  
)

CUSTOM_GlobalMarker <- list(
  # computation part
  "computation" = function(se, rows, columns, session, input, output, panel_name, funName, caches_env, mincell = 5, latentVars = c("nUMI", "nGene")){
    ns <- NS(paste0(panel_name, "_", funName))
    
    # Calculation -------------------------------------------------------------
    
    result = result = data.frame(ID=numeric(0));
    if (is.null(rows)){
      rows = rownames(se);
    }
    if (is.null(columns)){
      columns = colnames(se);
    }
    
    if (is.null(caches_env$seuratObj)){
      if(package_version(x = packageVersion(pkg = "Seurat")) >= package_version(x = "3.0.0")){
        seuratObj <- Seurat::as.Seurat(from = se)
      } else {
        seuratObj <- Seurat::as.seurat(from = se)
      }
      caches_env$seuratObj = seuratObj
    } else {
      seuratObj = caches_env$seuratObj
    }
    
    if(!any(is.null(input[[ns("compVar")]]), is.null(input[[ns("var")]]))){
      varName <- input[[ns("compVar")]];
      cellID1 <- rownames(subset(seuratObj@meta.data, get(varName) %in% input[[ns("var")]]))
      cellID2 <- rownames(subset(seuratObj@meta.data, ! get(varName) %in% input[[ns("var")]]))
      cellID1 <- cellID1[cellID1 %in% columns]
      cellID2 <- cellID2[cellID2 %in% columns]
      
      mycomp = paste0(input[[ns("compVar")]], ":", input[[ns("var")]], "_vs_", "other")
      method <- input[[ns("method")]]
      logfc <- input[[ns("logfc")]]
      minPct <- input[[ns("minPct")]]
      minDiffPct <- input[[ns("minDiffPct")]]
      downSample <- input[[ns("downsample")]]
      downSample = ifelse( downSample < 1, ceiling(downSample * (length(cellID1) + length(cellID2))), Inf) 
      if(length(cellID1) > mincell & length(cellID2) > mincell){
        # oldAssay = DefaultAssay(scDataClustered);
        # DefaultAssay(scDataClustered) = "RNA"
        #      browser();
        isolate({
          tmpDF = tryCatch({
            suppressWarnings(Seurat::FindMarkers(
              seuratObj, ident.1 = cellID1, ident.2 = cellID2,
              test.use = method, logfc.threshold = logfc,
              min.diff.pct = minDiffPct,
              min.pct = minPct, max.cells.per.ident = downSample
              # latent.vars = latentVars, 
            ))}, 
            error = function(x){return(NULL)}
          );
          if (!is.null(tmpDF)) {
            result = data.frame(
              # "comprison" = mycomp,
              tmpDF,
              row.names = rownames(tmpDF)
            )
          }  
        })
      }
    }
    return(result)
  },
  
  # UI + server logic part
  "UI" = function(se, pointsSelected, session, input, output, panel_name, funName, FUN, caches_env, UIonly = T) {
    ns <- NS(paste0(panel_name, "_", funName))
    possibleLevels = reactiveValues( var = NULL);
    
    # update function specific UI ---------------------------------------------
    dyUI = tagList(
      fluidRow(
        column(8, align = "left",
               tags$b(paste0(length(pointsSelected$columns), " cells have been selected"))
        )
      ),
      selectInput(ns("compVar"), "Variable of interest", choices = names(colData(se))),
      selectizeInput(ns("var"), "Group", choices = c("---"), selected = "---", multiple = TRUE),
      caches_env$sharedUI(paste0(panel_name, "_", funName))
    )
    
    if(!UIonly){
      
      observe({
        markType = ifelse(is.null(input[[ns("markerType")]]), "medium (fast)", input[[ns("markerType")]])
        selectedSetting <- caches_env$findMarkerSetting[[markType]]
        for(myid in names(selectedSetting)){
          updateSliderInput(session, ns(myid), value = selectedSetting[[myid]])
        }
      })
      
      observe({
        if(is.null(input[[ns("compVar")]])){
          possibleLevels$var = "---"
        } else {
          if(is.null(pointsSelected$rows)){
            rows <- rownames(se)
          } else {
            rows <- pointsSelected$rows
          }
          if(is.null(pointsSelected$columns)){
            columns <- colnames(se)
          } else {
            columns <- pointsSelected$columns
          }
          possibleLevels$var = colData(se[rows, columns])[, unique(input[[ns("compVar")]])];
          updateSelectizeInput(session, ns("var"), choices = possibleLevels$var, selected = possibleLevels$var[1])
          shinyjs::hide(paste0(panel_name, "_download"))
        }
      })
      
      observeEvent(input[[ns("submit")]], {
        shinyjs::disable(ns("submit"))
        updateActionButton(session, ns("submit"), "Running")
        shinyjs::hide(paste0(panel_name, "_download"))
        
        tmpDF = withProgress(message="Finding marker compared to all...", value=0,{
          do.call(FUN, c(list(se, pointsSelected$rows, pointsSelected$columns, session, input, output, panel_name, funName, caches_env)));
        })
        
        shinyjs::enable(ns("submit"))
        updateActionButton(session, ns("submit"), "Run")
        
        mycomp = paste0(input[[ns("compVar")]], ":", input[[ns("var")]], "_vs_", "other")
        
        output[[panel_name]] = DT::renderDataTable({
          validate(need(dim(tmpDF)[1]>0, "No marker identified"));
          # update download button
          outDF = tmpDF
          shinyjs::show(paste0(panel_name, "_download"))  
          output[[paste0(panel_name, "_download")]] = downloadHandler(
            filename = paste0(paste0(panel_name, "_download"), ".txt"),
            content = function(f){
              write.table(tmpDF, f, col.names = T, row.names = F, sep = "\t", quote = F);
            }
          )
          
          tmpDF$ID = rownames(tmpDF);
          origDF = data.frame("ID" = rownames(se), "SELECTED" = F, "ORDER" = 1:dim(se)[1], row.names = rownames(se));
          mergeDF = merge(origDF, tmpDF, by = "ID", all.x = T);
          rownames(mergeDF) = mergeDF$ID;
          mergeDF = mergeDF[order(mergeDF$ORDER), ];
          mergeDF$SELECTED = ifelse(is.na(mergeDF[, 4]), F, T);
          search_col = rep(list(NULL), ncol(mergeDF));
          search_col[[2]] = list(search = "[\"true\"]")
          
          DT::datatable(
            mergeDF, filter="top", rownames=TRUE, caption = mycomp,
            selection = list(mode="single"),
            options=list(
              search=list(smart=FALSE, regex=TRUE, caseInsensitive=FALSE),
              searchCols=c(list(NULL), search_col),
              columnDefs = list(list(visible=FALSE, targets=c(1,2,3))), # first column ID, second SELECTED and third column ORDER are invisible
              pageLength = 5, 
              stateSave = F,
              scrollX=TRUE
            )
          )
        })
      })  
    }
    return(dyUI)
  },
  
  "caches" = caches_env
)