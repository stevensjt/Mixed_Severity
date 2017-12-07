#This app calculates SDC for a high-severity polygon layer supplied by the user.
#Available at https://stevensjt.shinyapps.io/sdc_app/

#https://www.google.com/search?num=20&source=hp&q=mclapply+in+shiny&oq=mclapply+in+shiny&gs_l=psy-ab.3...2500141.2503852.0.2504036.39.20.3.0.0.0.194.1875.10j9.19.0....0...1.1.64.psy-ab..17.21.1823.0..0j0i131k1j0i10k1j0i22i30k1j0i22i10i30k1j33i160k1.0.lJ16CZUlx6c
#Good stuff here, maybe follow Robin's example:
#https://groups.google.com/forum/#!topic/shiny-discuss/EHXP2OpKLjk
#Good stuff here too, might need to split up functions and/or move them outside of server.R:
#https://shiny.rstudio.com/articles/scoping.html
 
require(devtools) #for session_info and install_version
#install_version("shiny", version = "1.0.5", repos = "http://cran.us.r-project.org"); install_version("rgdal", version = "1.2.5", repos = "http://cran.us.r-project.org"); install_version("parallel", version = "3.3.2", repos = "http://cran.us.r-project.org"); install_version("rgeos", version = "0.3.21", repos = "http://cran.us.r-project.org"); install_version("raster", version = "2.5.8", repos = "http://cran.us.r-project.org"); install_version("ggplot2", version = "2.2.1", repos = "http://cran.us.r-project.org"); install_version("RCurl", version = "1.95.4.8", repos = "http://cran.us.r-project.org")

#If running locally on Jens machine:
library(shiny) #originally 1.0.3
library(rgdal) # 1.2.5
library(parallel) #for "mclapply" version 3.3.2 #Comment out if functions placed outside app.R
library(rgeos) #for gBuffer version 0.3.21
library(raster) #for "area" version 2.5.8
library(ggplot2) #version 2.2.1
library(RCurl) #version 1.95.4.8
#library(rsconnect)
#rsconnect::deployApp('Code/SDC_app') #If you update the code, publish here.

#https://www.shinyapps.io/admin/#/dashboard
#http://shiny.rstudio.com/tutorial/
#Tutorial video:
#https://vimeo.com/rstudioinc/review/131218530/212d8a5a7a/#t=0m0s 
#Observe 1:11:37 or Delay 1:18:45

ui <- fluidPage(
  titlePanel("Stand-replacing Decay Coefficient (SDC) App"),
  helpText("The current application is a beta release and may not work for large files. The user must generate a multipart shapefile with a set of polygons representing the stand-replacing or high-severity burned area of interest. The shapefile must be in a metric projection. Uploading the shapefile will generate a plot of the area, and clicking 'run' will calculate SDC and display the SDC of the shapefile in question against 477 fires in California that burned from 1984-2015."),
  helpText("Citations for this application:"),
  helpText("Collins, B.M., Stevens, J.T., Miller, J.D., Stephens, S.L., Brown, P.M., North, M.P., 2017. Alternative characterization of forest fire regimes: incorporating spatial patterns. Landscape Ecology 32, 1543-1552."),
  helpText("Stevens, J. T., B. M. Collins, J. D. Miller, M. P. North, and S. L. Stephens. 2017. Changing spatial patterns of stand-replacing fire in California conifer forests. Forest Ecology and Management 406:28-36."),
  helpText("Contact Jens Stevens [stevensjt <at> berkeley.edu] with questions"),
  fileInput('inputdata', 'Input shapefile and accompanying 
            \nextensions (minimum .shp, .dbf, .prj and .shx)',
            accept=c('.shp','.dbf','.sbn','.sbx','.shx',".prj"), multiple=TRUE),
  plotOutput("perims"),
  numericInput("buf_inc", "Internal buffer distance interval, in m", 10),
  actionButton(inputId = "go",label = "Run"), 
  helpText("Only click -Run- once; may take some time (5-10 minutes) for large fires. Setting larger buffer distance will cut down on time."),
  textOutput("sdc.name"), #Checkme, commented out for troubleshooting
  #tableOutput("sdc.table"),
  plotOutput("sdc_hist")
  #htmlOutput("checks")
  
  )

server <- function(input, output) {
  
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


  output$perims <-  renderPlot({
    #Code modified from https://groups.google.com/forum/#!topic/shiny-discuss/FtI76rdvoHI
    if(!is.null(shape())){ #might need to replace shape() with input$inputdata
      ggplot(shape(),aes(x=long,y=lat, group = group)) + 
        geom_path(col = "darkred") +
        coord_equal() +
        labs(title = "your fire") + 
        theme_bw()
    }
    
  })
  
  sdc <- eventReactive(input$go, {
    ##Fill holes:
    ##Fill in any internal holes that are <0.81 ha
    
    fill_holes=function(hs_fire){
      hs_fire_p=slot(hs_fire, "polygons")
      holes <- lapply(hs_fire_p, function(x) sapply(slot(x, "Polygons"), slot, "hole"))
      areas <- lapply(hs_fire_p, function(x) sapply(slot(x, "Polygons"), slot, "area"))
      res <- lapply(1:length(hs_fire_p), 
                    function(i) slot(hs_fire_p[[i]], "Polygons")[!(holes[[1]]&areas[[1]]<8100)]
                    )
      #Select the polygons that are not holes. The "i" here is an artifact of the example code; the fires here only have one polygon ID so i=1 (it's a multipart polygon).
      IDs <- row.names(hs_fire)
      hs_fire_fill <- SpatialPolygons(lapply(1:length(res), function(i)
        Polygons(res[[i]], ID=IDs[i])), proj4string=CRS(proj4string(hs_fire)))
      return(hs_fire_fill)
      ##One consequence of this is that it's no longer a sp data frame, and there are some warnings. 
    }
    
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
                 mc.cores = 1) 
      #Setting mc.cores is key
      #Should on the server there are 16 cores available, 4 for personal computer.
      #This is why the app takes so long.
      #
      #Could mess with mc.preschedule = TRUE
      
      #Post-processing of long data frame:
      buf.list <- #Remove any NULL values (when buffer eliminated all HS areas)
        Filter(length,buf.list) 
      dist.table.sub$area_ha <- #Calculate area remaining at each buffer distance
        #0.0001 converts m2 to ha
        (as.vector(sapply(buf.list,raster::area))*0.0001)[1:nrow(dist.table.sub)] 
      if(any(is.na(dist.table.sub$area_ha))){
        #If there are NULL values for the area because the buffer got too wide, 
        #remove them. Set the first null value to 0.
        dist.table.sub$area_ha[which(is.na(dist.table.sub$area_ha))[1]] <- 0 
        if(any(is.na(dist.table.sub$area_ha))){
          #If there are STILL NULL values, delete the rest
          dist.table.sub <- #Delete the rest of the table rows with NA in area.
            dist.table.sub[-which(is.na(dist.table.sub$area_ha)),] 
        } 
      }
      dist.table.sub$prop.hs <- 
        #Calculate the proportion of the original HS area remaining,
        #at each buffer distance
        dist.table.sub$area_ha/dist.table.sub$area_ha[1]
      
      #Return the long data frame for the fire in question
      return(dist.table.sub)
    }
    
    ##Calculate stand-replacing decay coefficient (sdc)
    calculate.sdc=function(decay.table){
      #https://stat.ethz.ch/R-manual/R-devel/library/base/html/by.html
      m.list=with(decay.table,
                  by(decay.table,name,
                     function(x)
                       nls(prop.hs~1/(10^(sdc*width)),data=x,start=list(sdc=0.01))
                  )
      )
      sdc.table=data.frame(name=unique(as.character(decay.table$name)),
                           sdc=sapply(m.list,coef))
      out.table=merge(decay.table,sdc.table,sort=F)
      out.table$sdc.name=as.character(format(round(out.table$sdc,4),scientific=F))
      return(out.table)
    }
    
    hs_fire_filled <- fill_holes(shape())
    decay.table <- decay(hs_fire = hs_fire_filled)
    sdc.table <- calculate.sdc(decay.table) 
    sdc.name <- unique(sdc.table$sdc.name)
    #rm(list = c("hs_fire_filled","decay.table","sdc.table")) Changed, should be uncommented to remove objects that arent needed
    sdc.name #Changed to troubleshoot, should be sdc.name (trying out decay.table and sdc.table)
    #print(paste0("your sdc = ",unique(sdc$sdc.name)))
  })
  
  d <- eventReactive(input$go, {
    d<-read.csv(text = getURL("https://raw.githubusercontent.com/stevensjt/Mixed_Severity/master/Data/all_fires_ForAnalysis_weather.csv"))
  })
  
  output$sdc.name <- renderText({
    print(paste0("your sdc = ",sdc(), 
                 " ln(sdc) = ", round(log(as.numeric(sdc() ) ),4 ) ) )
  })
  
  output$sdc.table <- renderTable({
    sdc()
  })
  
  output$class_pct_plot<- renderPlot({
    ggplot(na.omit(d()[,c("BA90_pct","sdc","class")]),
           aes(x=BA90_pct,y=log(sdc),col=class))+
      geom_point()+
      geom_point(aes(x=0.5,y=log(as.numeric(sdc() ) ) ), col = "black" )+
      geom_smooth(method="lm")+
      labs(y= "ln(sdc)",x=" ", title = "a")+
      annotate("text", x=60, y=-3, label = paste("R^2 == ", 0.67), parse = TRUE) +
      theme_bw()+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=13),
            legend.position = 'none')
  })
  
  output$class_ha_plot<- renderPlot({
    ggplot(na.omit(d()[,c("firesize_ha","sdc","class")]),
           aes(x=log(firesize_ha),y=log(sdc),col=class))+
      geom_point()+
      geom_smooth(method="lm")+
      labs(y= " ", x=" ", title = "b")+
      annotate("text", x=10, y=-3, label = paste("R^2 == ", 0.22), parse = TRUE) +
      theme_bw()+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=13))
  })
  
  output$sdc_hist <- renderPlot({
    ggplot(d(),aes(x = log(sdc)))+
      geom_histogram(binwidth=0.1,fill="white",col="black")+
#      geom_vline(xintercept=log(c(0.0219,0.0068,0.0020,0.0006)),
#                 col=c("#377eb8","#4daf4a","#984ea3","#e41a1c"),
#                 size=1.2)+
      geom_vline(xintercept = log(as.numeric(sdc() ) ), size = 1.2, 
                 col = "black", lty=2)+
      annotate("text", x=log(as.numeric(sdc() ) )-0.05, y=15, label = "your fire",
               angle = 90)+
      labs(x="ln(SDC)",y="Number of fires",title=" ")+
      theme_bw()+
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=15) )
  })
  
  output$checks <- renderPrint({ #For troubleshooting
    capture.output(detectCores())
  })
  
  
}

shinyApp(ui = ui, server = server)