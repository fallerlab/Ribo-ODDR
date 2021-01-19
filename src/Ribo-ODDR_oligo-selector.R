
##################################################
#
#  ---  Ribo-ODDR - Ribo-seq focused Oligo Design pipeline for ---
#  ---     experiment-specific Depletion of Ribosomal RNAs     ---
#  Version 1.0 (See the ChangeLog.md file for changes.)
#
#  Copyright 2020   Faller Lab <fallerlab@gmail.com>
#                   Ferhat Alkan <f.alkan@nki.nl>
#                   William Faller <w.faller@nki.nl>
#  (See the AUTHORS file for other contributors.)
#
#  This file is part of the Ribo-ODDR pipeline.
#
#  Ribo-ODDR is a free software: you can
#  redistribute it and/or modify it under the terms of the
#  GNU General Public License as published by the
#  Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Ribo-ODDR is distributed in the hope that it will
#  be useful, but WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with the Ribo-ODDR Pipeline, see file LICENSE.
#  If not, see <http://www.gnu.org/licenses/>.
#################################################
require(shiny)
require(shinythemes)
require(DT)
require(ggplot2)

## Premade functions
# This one reads the csv file
get_df_from_csv <- function(rawfile){
  print(rawfile)
  raw_df <- read.csv(file = rawfile, stringsAsFactors = FALSE, comment.char = "#")
  tmp <- data.frame(oligoID=raw_df$oligoID, seq=raw_df$sequence, size=raw_df$length,
                    avg_dep_per=round(raw_df$avg_dep_per,1),
                    dep.score=round(as.numeric(raw_df$score),2), 
                    GC=round(as.numeric(raw_df$gc_content),2), 
                    energy=round(as.numeric(raw_df$eng_min),2), 
                    MFE=round(as.numeric(raw_df$mfe),2), 
                    BPper=round(as.numeric(raw_df$bp_per),2), 
                    foldedStructure=raw_df$structure,
                    off=as.numeric(raw_df$no_of_OFFs),
                    target=sapply(raw_df$target,FUN=function(x){return(strsplit(x,":")[[1]][1])}),
                    spos=as.numeric(sapply(raw_df$target,FUN=function(x){return(strsplit(strsplit(x,":")[[1]][2], "-")[[1]][1])})),
                    epos=as.numeric(sapply(raw_df$target,FUN=function(x){return(strsplit(strsplit(x,":")[[1]][2], "-")[[1]][2])})),
                    stringsAsFactors = FALSE)
  rownames(tmp)<-tmp$oligoID
  samples <- colnames(raw_df)[14:dim(raw_df)[2]]
  
  reads <- NULL
  rperc <- NULL
  tperc <- NULL
  for (i in 1:(dim(raw_df)[1])){
    reads <- rbind(reads, sapply(raw_df[i,samples], FUN = function(x){return(as.numeric(strsplit(x = strsplit(x = x, split = "_",fixed = TRUE)[[1]][1], split = '-', fixed = TRUE)[[1]][1]))}))
    rperc <- rbind(rperc, sapply(raw_df[i,samples], FUN = function(x){return(as.numeric(strsplit(x = strsplit(x = x, split = "_",fixed = TRUE)[[1]][2], split = '-', fixed = TRUE)[[1]][1]))}))
  }
  tmp$dep_reads <- rowMeans(reads)
  #tmp$avg_dep_per <- round(rowMeans(rperc),1)
  if (length(samples)==1){
    tmp[[paste("S1_", samples, "_reads_info", sep="")]] <- reads[,1]
  } else {
    for (i in 1:length(samples)){
      sample<-samples[i]
      tmp[[paste("S", as.character(i), "_", sample, "_reads_info", sep="")]] <- reads[,sample]
      tmp[[paste("S", as.character(i), "_", sample, "_rrna_per_info", sep="")]] <- rperc[,sample]
    }
  }
  return(tmp)
}

# Defining UI ----
ui <- fluidPage(
  tags$head(
  tags$style(
    HTML(".shiny-notification {
             position:fixed;
             top: calc(60%);
             left: calc(30%);
             }
             "
      )
    )
  ),
  
  sidebarPanel(
    # Give the version number here
    titlePanel(title = HTML('<i>Ribo-ODDR</i> <small>v1.0</small></br>
                            <small><small>Ribo-seq focused Oligo Design tool for Depleting Ribosomal RNAs</small></small></br>')),

    conditionalPanel(
      condition = "input.tabchoice == 'Welcome' | input.tabchoice == 'About' ",
      h4(HTML('<hr><small>Created by <a href="https://www.fallerlab.com/" target="_blank">Faller Lab</a></small>'))
    ),
    
    conditionalPanel(
      condition = "input.tabchoice == 'Instructions' | input.tabchoice == 'Select Oligos'",
      hr(),
      h4(HTML('<i><b>Ribo-ODDR:oligo-selector</b></i> - shiny app for choosing rRNA depletion oligos')),
      hr(),
      h4("Upload your Ribo-ODDR results"),
      h6("(restarts the session)"),
      fileInput("thecsv", label = HTML("Please choose your Ribo-ODDR <i>CSV</i> output"),
                accept = c("text/csv", "text/csv,text/plain", ".csv")),
      actionButton(inputId = "clearall", label = "clear and reload")
    ),
    
    conditionalPanel(
      condition = "input.tabchoice == 'Instructions'",
      h4(HTML('<hr><small>Created by <a href="https://www.fallerlab.com/" target="_blank">Faller Lab</a></small>'))
    ),
    
    conditionalPanel(
      condition = "input.tabchoice == 'Select Oligos'",
      hr(),
      h4("Filter designed oligos by"),
      sliderInput("o_l", label = h5("Oligo length (size)"), min = 0, max = 1, value = c(0,1)),
      sliderInput("o_tdp", label = h5("Average depleting potential"), min = 0, max = 1, value = c(0,1)),
      selectInput(inputId = "o_tar", choices = c("ALL"), label = h5("Target (ribosomal) RNA to deplete")),
      sliderInput("o_gc", label = h5("GC ratio"), min = 0, max = 1, value = c(0,1)),
      # sliderInput("o_be", label = h5("Oligo binding energy"), min = 0, max = 1, value = c(0,1)),
      sliderInput("o_off", label = h5("Number of potential off-targets"), min = 0, max = 1, value = c(0,1)),
      sliderInput("o_sc", label = h5("Depletion score"), min = 0, max = 1, value = c(0,1)),
      # sliderInput("o_eng", label = h5("Oligo self-fold MFE"), min = 0, max = 1, value = c(0,1)),
      hr(),
      h4(HTML('<hr><small>Created by <a href="https://www.fallerlab.com/" target="_blank">Faller Lab</a></small>'))
    ),
    width = 3
  ),
  
  mainPanel(
    tabsetPanel(id = "tabchoice",
      tabPanel(title = "Welcome",
               titlePanel(title = HTML('<i><b>Ribo-ODDR:oligo-selector</b></i> - shiny app for choosing rRNA depletion oligos'),
                          windowTitle = "Ribo-ODDR:oligo-selector"),
               hr(),
               h4(HTML('Welcome! <br><br>
                        This app has been designed as an extension to the <i>Ribo-ODDR</i> oligo design pipeline
                        with the hope that it will be useful for selecting your rRNA depletion oligos 
                        from the list of oligos reported by <i>Ribo-ODDR</i>. 
                        <br><br>
                        From the tabs above, please click "<i>Instructions</i>" and read the manual on how to use this app.
                        <br><br>
                        If you are already familiar with the app, continue by clicking "<i>Select Oligos</i>".
                        <br><br>
                        To learn more about <i>Ribo-ODDR</i>, please visit the "<i>About</i>" tab.
                        <br><br>
                        If you are using this tool, do not forget to cite our publication below.
                        <hr>
                        <i>Publication Citation comes here.</i>
                        <br>
                        <hr>Click <a href="https://github.com/ferhatalkan/Ribo-ODDR" target="_blank">here</a> 
                            to access the source code of <i>Ribo-ODDR</i> pipeline.
                       '))
               ),
      tabPanel(title = "Instructions",
               titlePanel(title = HTML('How to use <i>Ribo-ODDR:oligo-selector</i>?'),
                          windowTitle = "Ribo-ODDR:oligo-selector"),
               hr(),
               h4(HTML('<b>What is this app?</b>: It is an extension to the <i>Ribo-ODDR</i> oligo design pipeline, therefore, 
                        prior to using this app, you need to run <i>Ribo-ODDR.py</i> and generate your CSV output with oligo designs.
                        <hr>
                        <b>Aim of the app?</b>: To help you analyze <i>Ribo-ODDR</i> designed oligos and choose the high potential ones for maximum rRNA depletion.
                        <hr>
                        <b>How to use?</b><br>
                              <ul>
                                <li>Upload your <i>Ribo-ODDR</i> generated CSV file in the left panel.</li>
                                <li>Click "Select Oligos" tab above.</li>
                                <li>Browse through the "All oligo designs - Overview", the bottom-left table, and, if needed, use the filter options on the left panel to narrow down this list.</li>
                                <li>If you are in doubt about what these features are, please consult to <i>Ribo-ODDR</i> paper.</li>
                                <li>To see sample-specific depleting potentials of an oligo design, please click the oligo row in the bottom-left table.</li>
                                <li>When you click an oligo row in the left table, other oligo designs that are highly similar to it will appear in the bottom-right table.</li>
                                <li>If you are happy with the oligo chosen in the bottom-left list, click the selection button to add it to the list above. 
                                        This will automatically remove the highly similar oligos from the bottom-left list.</li>
                                <li>Continue your selection by browsing, inspecting and adding the oligos from the bottom-left list to the upper selection list 
                                until you are satisfied with your total depleting potential, shown as a plot next to the "Selected Oligos" list..</li>
                              </ul> 
                        <hr>
                        <b>Note that:</b>
                              <ul>
                                <li>Oligo features are presented to the user with no threshold recommendation. </li>
                                <li>However, you should always consider selecting high depleting potential oligos with mid-range GC content and minimum off-targets.</li>
                              </ul> 
                        <hr>
                        If you are using <i>Ribo-ODDR</i>, do not forget to cite the publication below.
                        <hr>
                        <i>Publication Citation comes here.</i>
                        <br>
                        <hr>Click <a href="https://github.com/ferhatalkan/Ribo-ODDR" target="_blank">here</a> 
                            to access the source code of <i>Ribo-ODDR</i> pipeline.
                       '))
               ),
      tabPanel(title = "Select Oligos",
               h3('Selected oligos'),
               fluidRow(
                 column(3, 
                        plotOutput('selectstats', height = 300)
                 ),
                 column(9, 
                        div(dataTableOutput(outputId = 'selection'), style = "font-size:75%"),
                        h5("Add your selections from the list below and see your overall depleting potential here.
                           When selection is done, copy your sequences and order :)")
                 )
               ),
               hr(),
               h3('All oligo designs - Overview'),
               
               fluidRow(
                 column(4,
                        actionButton(inputId = "addtoselection", width = "100%", icon("paper-plane"), 
                                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                                     label = "Selection Button")
                 ),
                 column(8,
                        plotOutput('selected_pots', height = 200)
                 )
               ),
               fluidRow(
                 column(8, h4('Oligos to choose from'),
                        verbatimTextOutput(outputId = "design_stats",placeholder = TRUE)),
                 column(4, h5('Below is the list of oligos that are highly similar to the oligo selected in the left table (they will be removed from the left table upon adding the left oligo to the selection)'),
                        verbatimTextOutput(outputId = "overlap_stats",placeholder = TRUE))
                 ),
               hr(),
               fluidRow(
                 column(8, div(DT::dataTableOutput('designs'), style = "font-size:75%")),
                 column(4, div(DT::dataTableOutput('overlap'), style = "font-size:60%"))
               )
               ),
      tabPanel(title = "About",
               titlePanel(title = HTML('<b><i>Ribo-ODDR</i></b>'),
                          windowTitle = "Ribo-ODDR:About"),
               hr(),
               h4(HTML('The most common limitation in Ribo-Seq experiments is the overabundance of ribosomal RNA
                (rRNA) fragments. <br><br>
                Various strategies are employed
                to tackle this issue, including the use of commercial
                rRNA depletion kits. Unfortunately, as these have largely
                been designed with RNASeq in mind, they may perform
                suboptimally. <br><br>
                In <i>Ribo-ODDR</i> publication, we show that the rRNA fragments
                generated via Ribo-seq vary significantly with differing
                experimental conditions, suggesting that a “one size fits all”
                approach may result in inefficient rRNA depletion. <br><br>
                In order to overcome this issue it is possible to use custom designed
                biotinylated oligos complementary to the most abundant
                rRNA fragments generated under specific experimental
                conditions, however currently no computational framework
                exists to aid the design of optimal oligos. <br><br>
                To meet these demands, we have developed Ribo-ODDR, an oligo
                design pipeline integrated with a user-friendly interface that
                assists in oligo selection for efficient experiment-specific
                rRNA depletion. <br><br>
                One can use Ribo-ODDR with preliminary or previously
                published data to identify the most abundant rRNA
                fragments, and calculate the rRNA depleting potential
                of potential oligos. <br><br>
                Selecting oligos with high potential, computed by Ribo-ODDR, lead to a significant increase in rRNA depletion, and
                increased sequencing depth as a result. 
                <hr>
                We hope that you find it useful. If you are interested in <i>Ribo-ODDR</i> publication, please see below.
                <hr>
                <i>Publication Citation comes here.</i>
                <hr>Click <a href="https://github.com/ferhatalkan/Ribo-ODDR" target="_blank">here</a> 
                            to access the source code of <i>Ribo-ODDR</i> pipeline.
                <hr><small>Created by <a href="https://www.fallerlab.com/" target="_blank">Faller Lab</a></small>
                '))
               )
    )
  )
)


# Defining server logic ----
server <- function(input, output, session) {
  design_ol_DF <- reactiveValues(dfWorking = data.frame()) # ALL designs
  filter_ol_DF <- reactiveValues(dfWorking = data.frame()) # Filtered based on left panel slides
  select_ol_DF <- reactiveValues(dfWorking = data.frame()) # Selected oligos
  overlap_ol_DF <- reactiveValues(dfWorking = data.frame()) # Oligos to show in overlap panel
  
  observeEvent(input$addtoselection, {
    showNotification('Selected row have been added to the upper selection table.', duration = 3)
    s = input$designs_rows_selected
    if (length(s)) {
      selectdf <- filter_ol_DF$dfWorking[s,]
      designdf <- design_ol_DF$dfWorking

      #Get not similar oligos (LeftEnd and RightEnd are both very dissimilar)
      ssize <- selectdf$size[1]
      slim <- max(10,round(ssize/3))
      sseq <- selectdf$seq[1]
      similar <- c()
      for (i in 1:(ssize-slim)){
        sub <- substr(sseq,i,i+slim)
        similar <- unique(c(similar,which(grepl(sub, designdf$seq, fixed = TRUE))))
      }

      designdf <- designdf[-similar,] # Update the designs

      # Update the globals
      overlap_ol_DF$dfWorking <-data.frame()
      select_ol_DF$dfWorking <- rbind(select_ol_DF$dfWorking, selectdf)
      design_ol_DF$dfWorking <- designdf
      filter_ol_DF$dfWorking <- designdf
    }
  })
  
  observeEvent(input$clearselection, {
    showNotification('Row selection cleared.', duration = 3)
  })
 
  observeEvent(input$clearall, {
    showModal(modalDialog(title = "Important message",
              "Selected oligos (top table) have been cleared and designed oligo list (bottom-left table) has been refreshed!"))
    inFile <- input$thecsv
    if (!is.null(inFile)) {
      df <- get_df_from_csv(inFile$datapath)
      design_ol_DF$dfWorking <- df
      filter_ol_DF$dfWorking <- design_ol_DF$dfWorking 
      select_ol_DF$dfWorking <- data.frame()
      overlap_ol_DF$dfWorking <- data.frame()
      
      updateSliderInput(session = session, inputId = "o_l", min = min(df$size), max = max(df$size), 
                        value = c(min(df$size), max(df$size)))
      updateSliderInput(session = session, inputId = "o_tdp", min = min(df$avg_dep_per), max = max(df$avg_dep_per), 
                        value = c(min(df$avg_dep_per), max(df$avg_dep_per)))
      updateSliderInput(session = session, inputId = "o_off", min = min(df$off), max = max(df$off), 
                        value = c(min(df$off), max(df$off)))
      # updateSliderInput(session = session, inputId = "o_be", min = min(df$energy), max = max(df$energy), 
      #                   value = c(min(df$energy), max(df$energy)))
      # updateSliderInput(session = session, inputId = "o_eng", min = min(df$MFE), max = max(df$MFE), 
      #                   value = c(min(df$MFE), max(df$MFE)))
      updateSliderInput(session = session, inputId = "o_gc", min = min(df$GC), max = max(df$GC), 
                        value = c(min(df$GC), max(df$GC)))
      updateSliderInput(session = session, inputId = "o_sc", min = min(df$dep.score), max = max(df$dep.score), 
                        value = c(min(df$dep.score), max(df$dep.score)))
      updateSelectInput(session = session, inputId = "o_tar", choices = c("ALL",unique(as.character(df$target))), selected = "ALL")
      
    }
  }) 
  
  observeEvent(input$thecsv, {
    inFile <- input$thecsv
    if (!is.null(inFile)) {
      updateActionButton(session = session, inputId = "addtoselection",
                         label = HTML("Selection Button<br><small>Click to add the selected oligo<br>to the upper selection list &<br>discard similar oligos from below</small>"))
      HTML("<small>Click to add the selected oligo<br>to the upper selection list &<br>discard similar oligos from below</small>")
      df <- get_df_from_csv(inFile$datapath)
      design_ol_DF$dfWorking <- df
      filter_ol_DF$dfWorking <- design_ol_DF$dfWorking 
      select_ol_DF$dfWorking <- data.frame()
      overlap_ol_DF$dfWorking <- data.frame()
      updateSliderInput(session = session, inputId = "o_l", min = min(df$size), max = max(df$size), 
                        value = c(min(df$size), max(df$size)))
      updateSliderInput(session = session, inputId = "o_tdp", min = min(df$avg_dep_per), max = max(df$avg_dep_per), 
                        value = c(min(df$avg_dep_per), max(df$avg_dep_per)))
      updateSliderInput(session = session, inputId = "o_off", min = min(df$off), max = max(df$off), 
                        value = c(min(df$off), max(df$off)))
      # updateSliderInput(session = session, inputId = "o_be", min = min(df$energy), max = max(df$energy), 
      #                   value = c(min(df$energy), max(df$energy)))
      # updateSliderInput(session = session, inputId = "o_eng", min = min(df$MFE), max = max(df$MFE), 
      #                   value = c(min(df$MFE), max(df$MFE)))
      updateSliderInput(session = session, inputId = "o_gc", min = min(df$GC), max = max(df$GC), 
                        value = c(min(df$GC), max(df$GC)))
      updateSliderInput(session = session, inputId = "o_sc", min = min(df$dep.score), max = max(df$dep.score), 
                        value = c(min(df$dep.score), max(df$dep.score)))
      updateSelectInput(session = session, inputId = "o_tar", choices = c("ALL",unique(as.character(df$target))), selected = "ALL")
    }
  })
  
  output$designs <-  DT::renderDataTable({
    df <- design_ol_DF$dfWorking
    df <- df[df$size >= input$o_l[1] & df$size <= input$o_l[2],]
    df <- df[df$avg_dep_per >= input$o_tdp[1] & df$avg_dep_per <= input$o_tdp[2],]
    df <- df[df$GC >= input$o_gc[1] & df$GC <= input$o_gc[2],]
    df <- df[df$dep.score >= input$o_sc[1] & df$dep.score <= input$o_sc[2],]
    if(dim(df)[1]>0) {
      if (sum(is.na(df$off))==0)
        df <- df[(df$off >= input$o_off[1] & df$off <= input$o_off[2]),]
      # if (sum(is.na(df$energy))==0)
      #   df <- df[(df$energy >= input$o_be[1] & df$energy <= input$o_be[2]),]
      # if (sum(is.na(df$MFE))==0)
      #   df <- df[(df$MFE >= input$o_eng[1] & df$MFE <= input$o_eng[2]),]
      if(input$o_tar!="ALL")
        df <- df[df$target==input$o_tar,]
    }
    filter_ol_DF$dfWorking <- df
    
    sketch = htmltools::withTags(table(
      class = 'display',
      thead(
        tr(
          th('', title = 'Oligo IDs as given in the uploaded CSV file'),
          th('avg_dep_per', title = 'Average depleting potential of the oligo. Click the oligo to see details.'),
          th('seq', title = 'Oligo sequence'),
          th('size', title = 'Oligo length'),
          th('GC', title = 'GC percentage for the oligo'),
          th('target', title = 'Target rRNA for the oligo'),
          # th('spos', title = 'Start position of the oligo target region'),
          # th('epos', title = 'End position of the oligo target region'),
          # th('energy', title = 'Binding energy of the perfect complementary oligo binding'),
          # th('MFE', title = 'Minimum Free Energy of the oligo self-folding at 37C'),
          th('off', title = 'Number of predicted off-targets for this oligo. To see what these off-targets are, please take a look at the OFF file generated by Ribo-ODDR pipeline.')
        )
      )
    ))
    if (dim(df)[1]>0)
      # datatable(filter_ol_DF$dfWorking[,c("avg_dep_per","seq","size","GC","target","spos","epos","energy","MFE","off")], 
      datatable(filter_ol_DF$dfWorking[,c("avg_dep_per","seq","size","GC","target","off")], 
                selection = "single", container = sketch,
                options = list(
                  lengthMenu = c(25, 50, 100, 150, 200),
                  order = list(list(1, 'desc'),list(6, 'desc')),
                  pageLength = 25
                ))
  })
  
  output$overlap <-  DT::renderDataTable({
    s = input$designs_rows_selected
    if (length(s)) {
      selectdf <- filter_ol_DF$dfWorking[s,]
      designdf <- design_ol_DF$dfWorking
      designdf <- designdf[designdf$oligoID!=selectdf$oligoID,]
      
      #Get not similar oligos (LeftEnd and RightEnd are both very dissimilar)
      ssize <- selectdf$size[1]
      slim <- max(10,round(ssize/3))
      sseq <- selectdf$seq[1]
      similar <- c()
      for (i in 1:(ssize-slim)){
        sub <- substr(sseq,i,i+slim)
        similar <- unique(c(similar,which(grepl(sub, designdf$seq, fixed = TRUE))))
      }
      
      df <- designdf[similar,]
      overlap_ol_DF$dfWorking <- df
      
      sketch = htmltools::withTags(table(
        class = 'display',
        thead(
          tr(
            th('', title = 'Oligo IDs as given in the uploaded CSV file'),
            th('avg_dep_per', title = 'Average depleting potential of the oligo. Click the oligo to see details.'),
            th('seq', title = 'Oligo sequence'),
            th('size', title = 'Oligo length'),
            th('GC', title = 'GC percentage for the oligo'),
            # th('target', title = 'Target rRNA for the oligo'),
            # th('spos', title = 'Start position of the oligo target region'),
            # th('epos', title = 'End position of the oligo target region'),
            # th('energy', title = 'Binding energy of the perfect complementary oligo binding'),
            # th('MFE', title = 'Minimum Free Energy of the oligo self-folding at 37C'),
            th('off', title = 'Number of predicted off-targets for this oligo. To see what these off-targets are, please take a look at the OFF file generated by Ribo-ODDR pipeline.')
          )
        )
      ))
      if (dim(df)[1]>0)
        datatable(df[,c("avg_dep_per","seq","size","GC","off")], 
                  selection = "single", container = sketch,
                  options = list(
                    lengthMenu = c(25, 50, 100, 150, 200),
                    order = list(list(1, 'desc'),list(5, 'desc')),
                    pageLength = 25
                  ))
    }
  })
  
  output$selection <-  DT::renderDataTable({
    sketch = htmltools::withTags(table(
      class = 'display',
      thead(
        tr(
          th('', title = 'Oligo IDs as given in the uploaded CSV file'),
          th('avg_dep_per', title = 'Average depleting potential of the oligo. Click the oligo to see details.'),
          th('seq', title = 'Oligo sequence'),
          th('size', title = 'Oligo length'),
          th('GC', title = 'GC percentage for the oligo'),
          th('target', title = 'Target rRNA for the oligo'),
          th('spos', title = 'Start position of the oligo target region'),
          th('epos', title = 'End position of the oligo target region'),
          th('energy', title = 'Binding energy of the perfect complementary oligo binding'),
          th('MFE', title = 'Minimum Free Energy of the oligo self-folding at 37C'),
          th('off', title = 'Number of predicted off-targets for this oligo. To see what these off-targets are, please take a look at the OFF file generated by Ribo-ODDR pipeline.')
        )
      )
    ))
    if (dim(select_ol_DF$dfWorking)[1]>0)
      datatable(select_ol_DF$dfWorking[,c('avg_dep_per',"seq","size","GC","target","spos","epos","energy","MFE","off")],
      # datatable(select_ol_DF$dfWorking[,c('avg_dep_per',"seq","size","GC","target","off")],
                selection = "none", container = sketch,
                options = list(paging = FALSE, searching = FALSE))
  })
  
  output$design_stats <- renderText({ 
    s = input$designs_rows_selected
    if (length(s)){
      one_df <- filter_ol_DF$dfWorking[s,]
      paste("Selected:", one_df$oligoID, "| Sequence:", one_df$seq,
            "| Prim.Target:", one_df$target, " Spos:", one_df$spos, " Epos:", one_df$epos,
            "| On average ", one_df$avg_dep_per, "% rRNA depletion (details plotted)", 
            "| No of off-targets:", as.character(one_df$off),
            "| Self-fold:", as.character(one_df$foldedStructure)
            )
    } else
      "Click an oligo below to summarize its detailed information here"
  })
  
  output$overlap_stats <- renderText({ 
    s = input$overlap_rows_selected
    if (length(s)){
      one_df <- overlap_ol_DF$dfWorking[s,]
      paste("Selected:", one_df$oligoID, "| Sequence:", one_df$seq,
            "| Prim.Target:", one_df$target, " Spos:", one_df$spos, " Epos:", one_df$epos,
            "| On average ", one_df$avg_dep_per, "% rRNA depletion (details plotted)", 
            "| No of off-targets:", as.character(one_df$off),
            "| Self-fold:", as.character(one_df$foldedStructure)
      )
    } else
      "Click an oligo below to summarize its detailed information here"
  })
  
  output$selected_pots <- renderPlot({
    s = input$designs_rows_selected
    if (length(s)) {
      melted <- NULL
      one_df <- filter_ol_DF$dfWorking[s,]
      for(tag in colnames(one_df)){
        if(endsWith(x = tag, suffix = "_info")){
          if(endsWith(x = tag, suffix = "_reads_info"))
            melted<-rbind(melted, data.frame(sample=gsub(pattern = "_reads_info",replacement = "",fixed = TRUE,x = tag), tag="reads", value=one_df[,tag]))
          else if(endsWith(x = tag, suffix = "_rrna_per_info"))
            melted<-rbind(melted, data.frame(sample=gsub(pattern = "_rrna_per_info",replacement = "",fixed = TRUE,x = tag), tag="rrna_per", value=one_df[,tag]))
        }
      }
      if(!is.null(melted))
        ggplot(melted[melted$tag=="rrna_per",],aes(x=sample,y=value,fill=sample))+
        geom_bar(stat = "identity")+ylab("Oligo depleting potential")+
        theme(axis.text.x = element_blank(), legend.title = element_blank(),legend.direction = "vertical",legend.position = "right") 
    }
  })
  
  output$selectstats <- renderPlot({
    melted <- NULL
    for(tag in colnames(select_ol_DF$dfWorking)){
      if(endsWith(x = tag, suffix = "_info")){
        if(endsWith(x = tag, suffix = "_reads_info"))
          melted<-rbind(melted, data.frame(sample=gsub(pattern = "_reads_info",replacement = "",fixed = TRUE,x = tag), tag="reads", value=select_ol_DF$dfWorking[,tag]))
        else if(endsWith(x = tag, suffix = "_rrna_per_info"))
          melted<-rbind(melted, data.frame(sample=gsub(pattern = "_rrna_per_info",replacement = "",fixed = TRUE,x = tag), tag="rrna_per", value=select_ol_DF$dfWorking[,tag]))
      }
    }
    if(!is.null(melted))
      ggplot(melted[melted$tag=="rrna_per",],aes(x=sample,y=value,fill=sample))+
      geom_bar(stat = "identity")+ylab("Total depleting potential")+
      ggtitle("Total Depletion")+ylim(0,100)+
      theme(axis.text.x = element_blank(), legend.direction = "vertical",legend.title = element_blank(),legend.position = "bottom")
  })
}

# Running the app ----
shinyApp(ui = ui, server = server)

