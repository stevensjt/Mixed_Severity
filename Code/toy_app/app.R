##mclapply patch example
library(shiny) 
library(parallel) 
library(rgdal)
library(rgeos)
library(raster)
#library(rsconnect)
#rsconnect::deployApp('Code/toy_app') #If you update the code, publish to web here.

ui <- fluidPage(
  fileInput('inputdata', 'Input shapefile and accompanying 
            \nextensions (minimum .shp, .dbf, .prj and .shx)',
            accept=c('.shp','.dbf','.sbn','.sbx','.shx',".prj"), multiple=TRUE),
  numericInput("buf_inc", "Enter buffer interval", 100),
  actionButton(inputId = "go",label = "Run"),
  tableOutput("decay_table")
)

server <- function(input, output) {
  ##Functions up front
  ##load shapefile
  shape <- reactive({
    #Read in the shapefile, store as output$shape (a reactive object)
    myshape<- input$inputdata
    if (is.null(myshape)) 
      return(NULL)       
    
    dir<-dirname(myshape[1,4]) #Get server path to uploaded files
    
    for ( i in 1:nrow(myshape)) {
      #Write full file path by combining server path and filename
      file.rename(myshape[i,4], paste0(dir,"/",myshape[i,1]))
    }
    getshp <- list.files(dir, pattern="*.shp", full.names=TRUE)
    readOGR(getshp)
  })
  
  ##core operation
  sdc <- eventReactive(input$go, {
    ##Implement internal buffering on high severity patches
    decay=function(hs_fire,buf_max=1000,buf_inc=input$buf_inc){
      #Set up long data frame:
      dist.table.sub <-
        data.frame(width=seq(0,buf_max,by=buf_inc),
                   area_ha=NA)
      #Core operation:
      buf.list <- #Apply the buffer at every width in X; returns list of length X.
        mclapply(X=-dist.table.sub$width,
                 FUN=gBuffer,
                 spgeom=hs_fire, 
                 byid = FALSE, id = NULL, quadsegs = 5, 
                 capStyle = "ROUND", joinStyle = "ROUND", mitreLimit = 1,
                 mc.cores = 4) 
      #Should be 16 cores for web version, 4 for personal computer.
      #Other arguments are passed to gBuffer
      
      buf.list <- #Remove any NULL values (when buffer eliminated all HS areas)
        Filter(length,buf.list) 
      dist.table.sub$area_ha <- #Calculate area remaining at each buffer distance
        #0.0001 converts m2 to ha
        (as.vector(sapply(buf.list,raster::area))*0.0001)[1:nrow(dist.table.sub)] 

      #Return the long data frame for the fire in question
      return(dist.table.sub)
    }
    
    decay.table <- decay(hs_fire = shape())
    decay.table
})
  
  ##pass resulting table from sdc() to output
  output$decay_table <- renderTable({
    sdc()
  })

}

shinyApp(ui = ui, server = server)