##mclapply toy example
library(shiny) #originally 1.0.3
library(parallel) #for "mclapply" version 3.3.2
library(rgdal)
library(rgeos)
#library(rsconnect)
#rsconnect::deployApp('Code/toy_app') #If you update the code, publish here.

times_ten <- function(x){x*10}

ui <- fluidPage(
  fileInput('inputdata', 'Input shapefile and accompanying 
            \nextensions (minimum .shp, .dbf, .prj and .shx)',
            accept=c('.shp','.dbf','.sbn','.sbx','.shx',".prj"), multiple=TRUE),
  numericInput("buf_inc", "Enter buffer interval", 100),
  actionButton(inputId = "go",label = "Run"),
  tableOutput("decay_table"),
  tableOutput("toy_table")
)

server <- function(input, output) {
  ##Functions up front
  ##toy function
  reactive_fun <- eventReactive(input$go, {
    get_table = function(x = input$buf_inc){
      df <- data.frame(initial = seq(1:x))
      result.list <- mclapply(X = df$initial,
                             FUN = times_ten,
                             mc.cores = 16)
      df$output <-sapply(result.list,identity)
      return(df)
    }
    t <- get_table()
    t
  })
  
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
    ##Decay profile: 
    ##Implement internal buffering on high severity patches, create "Long" dataset
    decay=function(hs_fire,buf_max=1000,buf_inc=input$buf_inc,name="your_fire"){
      #Set up long data frame:
      dist.table.sub <-
        data.frame(name=rep(name,(buf_max/buf_inc)+1),
                   width=seq(0,buf_max,by=buf_inc),
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
      
      #Post-processing of long data frame:
      buf.list <- #Remove any NULL values (when buffer eliminated all HS areas)
        Filter(length,buf.list) 
      dist.table.sub$area_ha <- #Calculate area remaining at each buffer distance
        #0.0001 converts m2 to ha
        (as.vector(sapply(buf.list,raster::area))*0.0001)[1:nrow(dist.table.sub)] 
#      if(any(is.na(dist.table.sub$area_ha))){
        #If there are NULL values for the area because the buffer got too wide, 
        #remove them. Set the first null value to 0.
#        dist.table.sub$area_ha[which(is.na(dist.table.sub$area_ha))[1]] <- 0 
#        if(any(is.na(dist.table.sub$area_ha))){
          #If there are STILL NULL values, delete the rest
#          dist.table.sub <- #Delete the rest of the table rows with NA in area.
#            dist.table.sub[-which(is.na(dist.table.sub$area_ha)),] 
#        } 
#      }

      #Return the long data frame for the fire in question
      return(dist.table.sub)
    }
    
    decay.table <- decay(hs_fire = shape())
    decay.table
})
  
  ##Begin output definitions
  
  output$decay_table <- renderTable({
    sdc()
  })
  
  output$toy_table <- renderTable({
    reactive_fun()
  })
  
}

shinyApp(ui = ui, server = server)