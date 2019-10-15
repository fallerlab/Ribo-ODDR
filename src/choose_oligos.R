
##################################################
#
#  --- Ribo-ODDR - Oligo Design for Depleting rRNAs ---
#  ---      in Ribosome profiling experiments       ---
#  Version 1.0 (See the ChangeLog.md file for changes.)
#
#  Copyright 2019   Ferhat Alkan <f.alkan@nki.nl>
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
require(rtracklayer)
require(ggplot2)

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
    titlePanel(title = HTML('<i>Ribo-ODDR</i> <small>v1.0</small></br>
                            <small><small>Ribo-seq focused Oligo Design tool for Depleting Ribosomal RNAs</small></small></br>
                            <hr>
                            <small><small>Created by <a href="https://www.fallerlab.com/" target="_blank">Faller lab</a></small></small>')),
    conditionalPanel(
      condition = "input.tabchoice == 'Select Oligos' | input.tabchoice == 'Instructions' ",
      hr(),
      h4(HTML('<i><b>Ribo-ODDR:oligo-selector</b></i> - shiny app for choosing rRNA depletion oligos')),
      hr(),
      h4("Upload your Ribo-ODDR results"),
      h6("(restarts the session)"),
      fileInput("thegff", label = HTML("Please choose the Ribo-ODDR <i>gff3</i> output"),
                accept = c("text/gff3", "text/gff3,text/plain", ".gff3")),
      actionButton(inputId = "clearall", label = "clear and reload")
    ),
    conditionalPanel(
      condition = "input.tabchoice == 'Select Oligos'",
      hr(),
      h4("Filter designed oligos by"),
      sliderInput("o_l", label = h5("size"), min = 0, max = 1, value = c(0,1)),
      sliderInput("o_tdp", label = h5("average depletion percentage"), min = 0, max = 1, value = c(0,1)),
      selectInput(inputId = "o_tar", choices = c("ALL"), label = h5("target ribosomal RNA")),
      sliderInput("o_gc", label = h5("GC ratio"), min = 0, max = 1, value = c(0,1)),
      sliderInput("o_be", label = h5("binding energy"), min = 0, max = 1, value = c(0,1)),
      sliderInput("o_off", label = h5("off-targets"), min = 0, max = 1, value = c(0,1)),
      sliderInput("o_sc", label = h5("depletion score"), min = 0, max = 1, value = c(0,1)),
      sliderInput("o_eng", label = h5("self-fold MFE"), min = 0, max = 1, value = c(0,1))
    ),
    width = 3
  ),
  mainPanel(
    tabsetPanel(id = "tabchoice",
      tabPanel(title = "Welcome",
               titlePanel(title = HTML('<i><b>Ribo-ODDR:oligo-selector</b></i> - shiny app for choosing rRNA depletion oligos'),
                          windowTitle = "Ribo-ODDR:oligo-selector"),
               hr(),
               h4("Some texthere")
               ),
      tabPanel(title = "Instructions"),
      tabPanel(title = "Select Oligos",
               h3('Selected oligos overview'),
               fluidRow(
                 column(3, 
                        plotOutput('selectstats', height = 300)
                 ),
                 column(9, 
                        div(tableOutput(outputId = 'selection'), style = "font-size:75%")
                 )
               ),
               hr(),
               h3('Ribo-ODDR designed oligos overview'),
               fluidRow(
                 column(6,
                        actionButton(inputId = "addtoselection", 
                                     label = HTML("Click here to add the selected oligo<br>to the upper selection list<br>& discard the overlapping oligos"))
                        ),
                 column(6,
                        verbatimTextOutput(outputId = "design_stats",placeholder = TRUE)
                        )
               ),
               hr(),
               div(DT::dataTableOutput('designs'), style = "font-size:75%")
               ),
      tabPanel(title = "About",
               titlePanel(title = HTML('<i><b>Ribo-ODDR:oligo-selector</b></i> - shiny app for choosing rRNA depletion oligos'),
                          windowTitle = "Ribo-ODDR:oligo-selector"),
               hr(),
               h4("Some texthere")
               )
    )
  )
)

get_df_from_gff <- function(raw_df){
  tmp <- data.frame(oligoID=raw_df$ID, seq=raw_df$seq, size=nchar(raw_df$seq),
                    dep.score=round(as.numeric(raw_df$score.1),2), GC=round(as.numeric(raw_df$gc_content),2), 
                    energy=round(as.numeric(raw_df$eng_min),2), 
                    MFE=round(as.numeric(raw_df$mfe),2), BPper=round(as.numeric(raw_df$bp_per),2), 
                    foldedStructure=raw_df$structure,
                    off=as.numeric(raw_df$no_of_OFFS),
                    target=raw_df$seqid, spos=as.numeric(raw_df$start), epos=as.numeric(raw_df$end))
  
  samples <- colnames(raw_df)[18:dim(raw_df)[2]]

  reads <- NULL
  rperc <- NULL
  tperc <- NULL
  for (i in 1:(dim(raw_df)[1])){
    reads <- rbind(reads, sapply(raw_df[i,samples], FUN = function(x){return(as.numeric(strsplit(x = strsplit(x = x, split = "-",fixed = TRUE)[[1]][1], split = '_', fixed = TRUE)[[1]][1]))}))
    rperc <- rbind(rperc, sapply(raw_df[i,samples], FUN = function(x){return(as.numeric(strsplit(x = strsplit(x = x, split = "-",fixed = TRUE)[[1]][2], split = '_', fixed = TRUE)[[1]][1]))}))
    tperc <- rbind(tperc, sapply(raw_df[i,samples], FUN = function(x){return(as.numeric(strsplit(x = strsplit(x = x, split = "-",fixed = TRUE)[[1]][3], split = '_', fixed = TRUE)[[1]][1]))}))
  }
  tmp$dep_reads <- rowMeans(reads)
  tmp$avg_dep_per <- round(rowMeans(rperc),0)
  tmp$total_dep_per <- round(rowMeans(tperc),0)
  for (i in 1:length(samples)){
     sample<-samples[i]
     tmp[[paste("S", as.character(i), "_", sample, "_reads_info", sep="")]] <- reads[,sample]
     tmp[[paste("S", as.character(i), "_", sample, "_rrna_per_info", sep="")]] <- rperc[,sample]
     tmp[[paste("S", as.character(i), "_", sample, "_tot_per_info", sep="")]] <- tperc[,sample]
  }
  return(tmp)
}

# Defining server logic ----
server <- function(input, output, session) {
  design_ol_DF <- reactiveValues(dfWorking = data.frame())
  filter_ol_DF <- reactiveValues(dfWorking = data.frame())
  select_ol_DF <- reactiveValues(dfWorking = data.frame())
  
  observeEvent(input$addtoselection, {
    showNotification('Selected row have been added to the upper selection table.', duration = 3)
    s = input$designs_rows_selected
    if (length(s)) {
      selectdf <- filter_ol_DF$dfWorking[s,]

      designdf <- design_ol_DF$dfWorking
      designdf <- designdf[designdf$oligoID!=selectdf$oligoID,]
      
      sametarget <- which(designdf$target==selectdf$target)
      overlap <- sapply(sametarget, FUN = function(x){
        if(length(intersect( ((designdf$spos[x]):(designdf$epos[x])), (selectdf$spos:selectdf$epos))) > (designdf$size[x]*0.5))
          return(TRUE)
        else
          return(FALSE)
      })
      designdf <- designdf[-sametarget[overlap],]
      
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
              "Selected oligos (top table) have been cleared and designed oligo list (bottom table) has been refreshed!"))
    inFile <- input$thegff
    if (!is.null(inFile)) {
      df <- get_df_from_gff(data.frame(rtracklayer::readGFF(inFile$datapath)))
      design_ol_DF$dfWorking <- df
      filter_ol_DF$dfWorking <- design_ol_DF$dfWorking 
      select_ol_DF$dfWorking <- data.frame()
      updateSliderInput(session = session, inputId = "o_l", min = min(df$size), max = max(df$size), 
                        value = c(min(df$size), max(df$size)))
      updateSliderInput(session = session, inputId = "o_tdp", min = min(df$avg_dep_per), max = max(df$avg_dep_per), 
                        value = c(min(df$avg_dep_per), max(df$avg_dep_per)))
      updateSliderInput(session = session, inputId = "o_off", min = min(df$off), max = max(df$off), 
                        value = c(min(df$off), max(df$off)))
      updateSliderInput(session = session, inputId = "o_be", min = min(df$energy), max = max(df$energy), 
                        value = c(min(df$energy), max(df$energy)))
      updateSliderInput(session = session, inputId = "o_eng", min = min(df$MFE), max = max(df$MFE), 
                        value = c(min(df$MFE), max(df$MFE)))
      updateSliderInput(session = session, inputId = "o_gc", min = min(df$GC), max = max(df$GC), 
                        value = c(min(df$GC), max(df$GC)))
      updateSliderInput(session = session, inputId = "o_sc", min = min(df$dep.score), max = max(df$dep.score), 
                        value = c(min(df$dep.score), max(df$dep.score)))
      updateSelectInput(session = session, inputId = "o_tar", choices = c("ALL",unique(as.character(df$target))), selected = "ALL")
      
    }
  }) 
  
  observeEvent(input$thegff, {
    inFile <- input$thegff
    if (!is.null(inFile)) {
      df <- get_df_from_gff(data.frame(rtracklayer::readGFF(inFile$datapath)))
      design_ol_DF$dfWorking <- df
      filter_ol_DF$dfWorking <-  design_ol_DF$dfWorking 
      select_ol_DF$dfWorking <- data.frame()
      updateSliderInput(session = session, inputId = "o_l", min = min(df$size), max = max(df$size), 
                        value = c(min(df$size), max(df$size)))
      updateSliderInput(session = session, inputId = "o_tdp", min = min(df$avg_dep_per), max = max(df$avg_dep_per), 
                        value = c(min(df$avg_dep_per), max(df$avg_dep_per)))
      updateSliderInput(session = session, inputId = "o_off", min = min(df$off), max = max(df$off), 
                        value = c(min(df$off), max(df$off)))
      updateSliderInput(session = session, inputId = "o_be", min = min(df$energy), max = max(df$energy), 
                        value = c(min(df$energy), max(df$energy)))
      updateSliderInput(session = session, inputId = "o_eng", min = min(df$MFE), max = max(df$MFE), 
                        value = c(min(df$MFE), max(df$MFE)))
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
    df <- df[df$off >= input$o_off[1] & df$off <= input$o_off[2],]
    df <- df[df$energy >= input$o_be[1] & df$energy <= input$o_be[2],]
    df <- df[df$MFE >= input$o_eng[1] & df$MFE <= input$o_eng[2],]
    df <- df[df$GC >= input$o_gc[1] & df$GC <= input$o_gc[2],]
    df <- df[df$dep.score >= input$o_sc[1] & df$dep.score <= input$o_sc[2],]
    if(input$o_tar!="ALL")
      df <- df[df$target==input$o_tar,]
    filter_ol_DF$dfWorking <- df
    if (dim(df)[1]>0)
      filter_ol_DF$dfWorking[,c("avg_dep_per","seq","size","GC","target","spos","epos","dep.score","energy","MFE","BPper","off")]
  }, selection = "single", options = list(
    lengthMenu = c(25, 50, 100, 150, 200),
    order = list(list(1, 'desc'),list(8, 'desc')),
    pageLength = 25
  ))
  
  output$design_stats <- renderText({ 
    s = input$designs_rows_selected
    if (length(s)) 
      as.character(filter_ol_DF$dfWorking[s,"foldedStructure"])
    else
      "Selected oligo information will be summarized here"
  })
  
  output$selection <-  renderTable({
    if (dim(select_ol_DF$dfWorking)[1]>0)
      select_ol_DF$dfWorking[,c("seq","size","GC","target","spos","epos","dep.score","energy","MFE","BPper","off")]
  },spacing = "s",width = "auto")
  
  output$selectstats <- renderPlot({
    melted <- NULL
    for(tag in colnames(select_ol_DF$dfWorking)){
      if(endsWith(x = tag, suffix = "_info")){
        if(endsWith(x = tag, suffix = "_reads_info"))
          melted<-rbind(melted, data.frame(sample=gsub(pattern = "_reads_info",replacement = "",fixed = TRUE,x = tag), tag="reads", value=select_ol_DF$dfWorking[,tag]))
        else if(endsWith(x = tag, suffix = "_rrna_per_info"))
          melted<-rbind(melted, data.frame(sample=gsub(pattern = "_rrna_per_info",replacement = "",fixed = TRUE,x = tag), tag="rrna_per", value=select_ol_DF$dfWorking[,tag]))
        else if(endsWith(x = tag, suffix = "_tot_per_info"))
          melted<-rbind(melted, data.frame(sample=gsub(pattern = "_tot_per_info",replacement = "",fixed = TRUE,x = tag), tag="tot_per", value=select_ol_DF$dfWorking[,tag]))
      }
    }
    if(!is.null(melted))
      ggplot(melted[melted$tag=="rrna_per",],aes(x=sample,y=value,fill=sample))+
      geom_bar(stat = "identity")+ylab("rRNA depletion percentage")+
      theme(axis.text.x = element_blank(), legend.direction = "vertical",legend.position = "bottom")
  })
}

# Running the app ----
shinyApp(ui = ui, server = server)

