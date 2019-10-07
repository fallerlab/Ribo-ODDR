# Load packages
library(shiny)
library(shinythemes)
library(DT)

# Defining UI ----
ui <- fluidPage(
  tags$head(
  tags$style(
    HTML(".shiny-notification {
             position:fixed;
             top: calc(50%);
             left: calc(30%);
             }
             "
      )
    )
  ),
  titlePanel(title = HTML("<i><b>Ribo-ODDR:oligo-selector</b></i> - shiny app for choosing the rRNA depletion oligos designed by <i>Ribo-ODDR</i>"),
             windowTitle = "Ribo-ODDR:oligo-selector"),
  fluidRow(
    column(3, 
           h4("Upload your Ribo-ODDR results"),
           h6("(restarts the session)"),
           fileInput("thegff", label = HTML("Please choose the Ribo-ODDR <i>gff3</i> output"),
                     accept = c("text/gff3", "text/gff3,text/plain", ".gff3")),
           actionButton(inputId = "clearall", label = "clear the oligo selection")
          ),
    column(6, 
           h4('Selected Oligos'),
           DT::dataTableOutput('selection')
          ),
    column(3, 
           plotOutput('selectstats', height = 300)
          )
  ),
  hr(),
  h3('Ribo-ODDR oligo design overview'),
  fluidRow(
    column(3, 
           h4("Filter oligos by"),
           sliderInput("o_l", label = h5("length"), min = 0, max = 1, value = c(0,1)),
           sliderInput("o_gc", label = h5("GC ratio"), min = 0, max = 1, value = c(0,1)),
           sliderInput("o_be", label = h5("binding energy"), min = 0, max = 1, value = c(0,1)),
           sliderInput("o_off", label = h5("off-targets"), min = 0, max = 1, value = c(0,1)),
           sliderInput("o_sc", label = h5("depletion score"), min = 0, max = 1, value = c(0,1)),
           sliderInput("o_eng", label = h5("self-fold MFE"), min = 0, max = 1, value = c(0,1)),
           fluidRow(column(6, actionButton(inputId = "addtoselection", label = "add selected")),
                    column(6, actionButton(inputId = "clearselection", label = "unselect rows")))
          ),
    column(7, 
           DT::dataTableOutput('designs')
          ),
    column(2, 
           verbatimTextOutput('designstats')
          ))
)

get_df_from_gff <- function(raw_df){
  tmp <- data.frame(oligoID=raw_df$ID, seq=raw_df$seq, length=nchar(raw_df$seq),
                    dep.score=round(as.numeric(raw_df$score.1),2), GC=round(as.numeric(raw_df$GC_content),2), 
                    energy=round(as.numeric(raw_df$Emin),2), 
                    MFE=round(as.numeric(raw_df$MFE),2), BPper=round(as.numeric(raw_df$BPper),2), 
                    foldedStructure=raw_df$structure,
                    off=as.numeric(raw_df$NOofOFFS),
                    target=raw_df$seqid, spos=as.numeric(raw_df$start), epos=as.numeric(raw_df$end))
  return(tmp)
}

# Defining server logic ----
server <- function(input, output, session) {
  design_ol_DF <- reactiveValues(dfWorking = data.frame())
  filter_ol_df <- reactiveValues(dfWorking = data.frame())
  select_ol_DF <- reactiveValues(dfWorking = data.frame())
  
  observeEvent(input$addtoselection, {
    showNotification('Selected rows have been added to the selection table.', duration = 3)
    s = input$designs_rows_selected
    if (length(s)) {
      selectdf <- filter_ol_df$dfWorking[s,]

      designdf <- design_ol_DF$dfWorking
      designdf <- designdf[designdf$oligoID!=selectdf$oligoID,]
      
      sametarget <- which(designdf$target==selectdf$target)
      overlap <- sapply(sametarget, FUN = function(x){
        if(length(intersect( ((designdf$spos[x]):(designdf$epos[x])), (selectdf$spos:selectdf$epos))) > (designdf$length[x]*0.5))
          return(TRUE)
        else
          return(FALSE)
      })
      designdf <- designdf[-sametarget[overlap],]
      
      select_ol_DF$dfWorking <- rbind(select_ol_DF$dfWorking, selectdf)
      design_ol_DF$dfWorking <- designdf
      filter_ol_df$dfWorking <- designdf
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
      filter_ol_df$dfWorking <- design_ol_DF$dfWorking 
      select_ol_DF$dfWorking <- data.frame()
      updateSliderInput(session = session, inputId = "o_l", min = min(df$length), max = max(df$length), 
                        value = c(min(df$length), max(df$length)))
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
    }
  }) 
  
  observeEvent(input$thegff, {
    inFile <- input$thegff
    if (!is.null(inFile)) {
      df <- get_df_from_gff(data.frame(rtracklayer::readGFF(inFile$datapath)))
      design_ol_DF$dfWorking <- df
      filter_ol_df$dfWorking <-  design_ol_DF$dfWorking 
      select_ol_DF$dfWorking <- data.frame()
      updateSliderInput(session = session, inputId = "o_l", min = min(df$length), max = max(df$length), 
                        value = c(min(df$length), max(df$length)))
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
    }
  })
  
  output$selection <-  DT::renderDataTable({
    select_ol_DF$dfWorking
  }, selection = "single")
  
  output$designs <-  DT::renderDataTable({
    df <- design_ol_DF$dfWorking
    df <- df[df$length >= input$o_l[1] & df$length <= input$o_l[2],]
    df <- df[df$off >= input$o_off[1] & df$off <= input$o_off[2],]
    df <- df[df$energy >= input$o_be[1] & df$energy <= input$o_be[2],]
    df <- df[df$MFE >= input$o_eng[1] & df$MFE <= input$o_eng[2],]
    df <- df[df$GC >= input$o_gc[1] & df$GC <= input$o_gc[2],]
    df <- df[df$dep.score >= input$o_sc[1] & df$dep.score <= input$o_sc[2],]
    filter_ol_df$dfWorking <- df
    filter_ol_df$dfWorking
  }, selection = "single")
}

#df<-data.frame(rtracklayer::readGFF("../../../design_rounds/Round2_Jul19/Ribo9th_and_10th_mix/oligos.gff3"))

# Running the app ----
shinyApp(ui = ui, server = server)

