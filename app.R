library(shiny)
library(pathview)
library(utils)

ui<- fluidPage(
  titlePanel("File upload"),
  sidebarLayout(
    sidebarPanel(
      fileInput("csv",
                label="Upload CSV here",
                multiple = FALSE),
      textInput("pathwayID", label = "Write pathway ID here"),
      textInput("species", label = "Write species here"),
      textInput("IDtype", label = "Write gene ID type here"),
      textInput("outname", label = "Write output name here"),
      numericInput("limit", label = "Write legend limit here",value = 10),
      checkboxGroupInput("ifmultiple", label = "Choice single sample or multiple samples",
                         c("Single sample"="png",
                           "Multiple samples"="multi.png")),
      submitButton("Update View", icon("refresh")),
      
    ),
    mainPanel(
      imageOutput("PathviewGraph")
    )
  )
)

server<-function(input,output){
  
  output$PathviewGraph<-renderImage({
    
    if (is.null(input$csv)) {
      return(NULL) } else {
      infile<-input$csv
      inputdata<-read.csv(infile$datapath, row.names = 1, header = 1)
      pathview(gene.data = inputdata, pathway.id = input$pathwayID, species = input$species, gene.idtype = input$IDtype,limit = input$limit ,kegg.native = T,out.suffix = input$outname, both.dirs = list(gene=F,cpd=F))
      pathway<-paste(input$species,input$pathwayID,sep = "")
      outfile<-paste(pathway,input$outname,input$ifmultiple, sep=".")
      return(list(src = outfile)) }
    
    
    
    
   },deleteFile = FALSE)
}
  

shinyApp(ui = ui, server = server)