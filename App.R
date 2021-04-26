library(raster)
library(EBImage)
library(raster)

setwd(paste0(getwd()))


# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Torres Lab Cell Segmentation Analyzer"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Text:
      
      h4('Upload images to analyze'),
      
      # Input: Image
      fileInput("upload", "Upload", accept = "image/png/tif", multiple = TRUE),
      
      # Text:
      h4('Dimensions in pixels to consider an object in segmentation'),
      
      # Input: Box pixels to create a mask by
      numericInput("x_mask", 'X pixels', value = 40, min = 1),
      
      # Input: Box pixels to create a mask by
      numericInput("y_mask", 'Y pixels', value = 40, min = 1),
      
      #Text:
      h4('Segement Current Image'),
      
      # Input: Button that starts analysis of chosen image
      actionButton('segc1', label = 'GO'),
      
      # Text:
      h4('Download Calculated Features for All Uploaded Images'),
      
      # Output: Downloadbutton that does cell segmentation for all images uploaded and downloads information
      downloadButton('downloadData', 'Download')
      
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel( tabsetPanel(
      
      tabPanel("Nuclear Cell Segmentation",
               
               h4('Select which image to test cell segmentation which. Remember that the settings applied to this image will be applied to all images when the data is downloaded. '),
               
               # Input: Select an Image
               selectInput('RenderRaster', label = 'Select Image', choices = 'No choices here yet'),
               
               # Text:
               h3('Raw Image:'),
               
               # Output: Render Current Image
               plotOutput('image'),
               
               # Text:
               h3('Segmented Image:'),
               
               # Output: Render the cell segmentation of current image
               plotOutput('seg.im')),
      
      tabPanel("Calculated Features",
               
               # Output: Render datatable of features
               dataTableOutput('features'))
      
    )
  )
))

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  
  
  # Reactive Object: Upload an image when uploaded
  observeEvent(input$upload, {
    
    updateSelectInput(session, "RenderRaster", 
                      label = "Select", 
                      choices = input$upload$datapath)
    
  })
  
  # Reactive Object: Upload an image
  img <- reactive({
    
    if(is.null(input$upload$datapath) == TRUE){
      return(NULL)
    }else{
      
      raster.l <- lapply(input$upload$datapath, readImage)
      names(raster.l) <- input$upload$datapath
      
      return(raster.l)
    }
    
    # imported_raster = pbapply::pblapply(files, raster) #create a raster object
    
  })
  
  # Reactive Object: Current selected image
  im.c <- reactive({
    
    if (is.null(img()) == TRUE){
      return(NULL)
    }else if (input$RenderRaster != 'No choices here yet'){
      
      im <- img()
      im.c <- im[names(im) == input$RenderRaster]
      
      return(im.c)
      
    }
    
    
  })
  
  
  # Output: plot the raster of the image selected in the input
  output$image <- renderPlot({
    
    if (is.null(img()) == TRUE){
      return(NULL)
    }else if (input$RenderRaster != 'No choices here yet'){
      
      y <- EBI.im()
      plot(y, all = TRUE)

    }
    
  })
  
  # Reactive Object: Create an EBI image
  EBI.im <- reactive({
    
    # Normalize by the maximum intensity 
    im <- im.c()
    im <- im[[1]]
    #im <- raster::as.matrix(im)
    im <- im/max(im)
    
    # working with EBIimage
    im <- EBImage::Image(data = im, dim = dim(im)) # make it an EBI object missing color greadient
    
  })
  
  # Reactive Object: Create object segmentation mask
  EBI.mask.im <- reactive({
    
    # Cell segmentation
    # First, we segment the nuclei using thresh, fillHull, bwlabel and opening.
    y <- EBI.im()
    
    x = thresh(y, input$x_mask, input$y_mask, 0.05)
    x = opening(x, makeBrush(5, shape='disc'))
    x = fillHull(x) # Fills holes in objects
    x = bwlabel(x) # labels images
    
    #Watershed image segmentation - this improves segmentation since it looks at topology of image
    x = watershed(distmap(x), 2)
    
    return(x)
    
  })
  
  
  # Reactive Object: Datatable that summarizes all features from objects in the window
  calc.features <- reactive({
    
    # Import necessary information
    x = EBI.mask.im()
    y = EBI.im()
    
    #compute features 
    ft = computeFeatures(x, y, xname="nucleus")
    
    return(ft)
    
  })
  
  
  
  # Action Button: When clicked it will segment the current image
  observeEvent(input$segc1, {
    
    # Output: create a map that shows the cell segmentation
    output$seg.im <- renderPlot({
      
      if (is.null(img()) == TRUE){
        return(NULL)
      }else if (input$RenderRaster != 'No choices here yet'){
        
        x = EBI.mask.im()
        y = EBI.im()
        
        segmented = paintObjects(x, y, col='#ff00ff')
        plot(segmented, all = TRUE)
      }
      
      
    })
    
    # Output: Render datatable with all the features <- can we create graphs
    output$features <- renderDataTable({
      
      calc.features()
      
    })
    
    # Output: csv file with all the features for all images
    output$downloadData <- downloadHandler(
      filename = function() {
        paste('data-', Sys.Date(), '.csv', sep='')
      },
      content = function(con) {
        
        if (is.null(img()) == TRUE){
          data = NULL
        }else if (input$RenderRaster != 'No choices here yet'){
          
          im.all <- img()
          data = data.frame()
          
          withProgress(message = 'Calculating Features', value = 0, {
            # Number of times we'll go through the loop
            n <- length(im.all)
            
            for (i in 1:n) {
              # Each time through the loop, add another row of data. This is
              # a stand-in for a long-running computation.
              
              im.c <- im.all[i]
              n.im.c <- names(im.all)[i]
              #print(n.im.c)
              
              # Normalize by the maximum intensity 
              im <- im.c[[1]]
              im <- raster::as.matrix(im)
              im <- im/max(im)
              
              # working with EBIimage
              y <- EBImage::Image(data = im, dim = dim(im)) # make it an EBI object missing color greadient
              
              x = thresh(y, input$x_mask, input$y_mask, 0.05)
              x = opening(x, makeBrush(5, shape='disc'))
              x = fillHull(x) # Fills holes in objects
              x = bwlabel(x) # labels images
              
              #Watershed image segmentation - this improves segmentation since it looks at topology of image
              x = watershed(distmap(x), 2)
              
              #compute features 
              ft = as.data.frame(computeFeatures(x, y, xname="nucleus"))
              ft$file = rep(n.im.c, nrow(ft))
              
              #add to final
              data = rbind(data, ft)
              
              
              # Increment the progress bar, and update the detail text.
              incProgress(1/n, detail = paste("Doing part", i))
              
              # Pause for 0.1 seconds to simulate a long computation.
              Sys.sleep(0.1)
            }
          })
          
        }
        
        write.csv(data, con)
      }
    )
    
    
  })
  
  
  

  
  
  
  
  
  
}

shinyApp(ui, server)