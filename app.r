source("global.R")
options(shiny.maxRequestSize=30*1024^2) 
ui <- shinyUI(navbarPage("Application to Explore Spatial Regression Techniques", theme = shinytheme("lumen"),
                         windowTitle = "Cap App",
                         tabPanel("Load Data",
                                  sidebarLayout(
                                    sidebarPanel(
                                      conditionalPanel(condition="input.tabs1 == null || input.tabs1 =='Data'",
                                      fileInput("file",'Input file 1',accept=c('.csv', '.shp','.dbf','.sbn','.sbx','.shx',".prj",'.geojson'), multiple=TRUE),
                                      tags$hr(),
                                      fileInput("file2",'Input file 2',accept=c('.csv', '.shp','.dbf','.sbn','.sbx','.shx',".prj",'.geojson'), multiple=TRUE)),
                                      conditionalPanel(condition="input.tabs1=='GWRmap'",
                                                         selectInput("method", "Method",
                                                                   c("cv", "aic")
                                                       ),
                                                       selectInput("inVar", label = h5("Select Variable for GWR Analysis"), ""),
                                                       actionButton("rungwr", "Run GWR"),
                                                       downloadButton('downloadGWR', 'Download GWR Map', class="dlButton"),
                                                         imageOutput("plot")),
                                      conditionalPanel(condition = "input.tabs1=='Nearest Neighbor Map'",
                                                       radioButtons('radio', label = 'Select a neighborhood type:',  
                                                                    choices = list("Queen's case contiguity" = 1, "Rook's case contiguity" = 2, 
                                                                                   'K-nearest neighbors' = 3, 'Distance' = 4), 
                                                                    selected = 1), 
                                                       conditionalPanel(
                                                         condition = "input.radio == 3", 
                                                         sliderInput("knn_slider", 'Select number of neighbors', 
                                                                     min = 1, max = 217, value = 8)
                                                       ), 
                                                       conditionalPanel(
                                                         condition = "input.radio == 4", 
                                                         sliderInput("dist_slider", "Select a distance threshold in km", 
                                                                     min = 0.1, max = 30, step = 0.1, value = 1)
                                                       )),
                                      conditionalPanel(condition = "input.tabs1=='OLS Residuals Map'",
                                      selectInput("colors", "Color Scheme",
                                                  rownames(subset(brewer.pal.info, category %in% c("seq", "div"))), selected = "YlOrRd"),
                                      actionButton("runols", "Run OLS"),
                                      downloadButton('downloadOLS', 'Download OLS Map', class="dlButton"),
                                      verbatimTextOutput("olsText"),
                                      actionButton("addmoran", "Add Moran's I Summary"),
                                      actionButton("removemoran", "Remove Moran's I Summary")
                                      )
                                      
                                      
                                      ),
                                    mainPanel(
                                      tags$style(type = "text/css", "#map {height: calc(100vh - 80px) !important;}"),
                                      uiOutput("tb")
                                      
                                      # use below code if you want the tabset programming in the main panel. If so, then tabset will appear when the app loads for the first time.
                                      #       tabsetPanel(tabPanel("Summary", verbatimTextOutput("sum")),
                                      #                   tabPanel("Data", tableOutput("table")))
                                    )
                                    
                                  )
                         )
)
)

server <- shinyServer(function(input,output, session){
  dsnames <- c()
  # This reactive function will take the inputs from UI.R and use them for read.table() to read the data from the file. It returns the dataset in the form of a dataframe.
  # file$datapath -> gives the path of the file
  data <- reactive({
    
    myshape<- input$file
    if (is.null(myshape)) 
      return(NULL)       
    
    dir<-dirname(myshape[1,4])
    
    for ( i in 1:nrow(myshape)) {
      file.rename(myshape[i,4], paste0(dir,"/",myshape[i,1]))}
    
    getshp <- list.files(dir, pattern="*.shp", full.names=TRUE)
    shape<-st_transform(st_read(getshp), "+proj=longlat +datum=WGS84 +no_defs")
    shape$id <- 1:nrow(shape) ## This is for spatial neighbors
    shape[is.na(shape)] <- 0
    shape
  })
  data2 <- reactive({
    myshape<- input$file2
    if (is.null(myshape)) 
      return(NULL)       
    
    dir<-dirname(myshape[1,4])
    
    for ( i in 1:nrow(myshape)) {
      file.rename(myshape[i,4], paste0(dir,"/",myshape[i,1]))}
    
    getshp <- list.files(dir, pattern="*.shp", full.names=TRUE)
    shape<-st_transform(st_read(getshp), "+proj=longlat +datum=WGS84 +no_defs")
    shape
  })
  
  aq <- reactive({
    data <- as(data(), "Spatial")
    aq <- poly2nb(data, queen = TRUE, row.names = data$id)
    aq
  })
  
  ar <- reactive({
    data <- as(data(), "Spatial")
    ar <- poly2nb(data, queen = FALSE, row.names = data$id)
    ar
  })
  
  knn <- reactive({
    if (input$radio == 3) {
      data <- as(data(), "Spatial")
      k <- knearneigh(coordinates(data), k = input$knn_slider, longlat = TRUE)
      return(k$nn)
    } else {
      return(NULL)
    }
  })
  
  dist <- reactive({
    if (input$radio == 4) {
      data <- as(data(), "Spatial")
      return(dnearneigh(coordinates(data), 0, input$dist_slider, longlat = TRUE))
    } else {
      return(NULL)
    }
  })
  
  spweights <- reactive({
    if (input$radio == 1) {
      return(nb2listw(aq()))
    } else if (input$radio == 2) {
      return(nb2listw(ar()))
    }
  })
  
  morantest <- reactive({
    req(runRegression())
    req(spweights())
    x <- lm.morantest(runRegression(), spweights())
    x
  })
  
  output$moran <- renderPrint({
    if(!is.null(input$indep)){
      morantest()
    } else {
      print(data.frame(Warning="Please select Model Parameters."))
    }
  })
  
  observeEvent(input$addmoran, {
    insertUI(
      selector = "#addmoran",
      where = "afterEnd",
      ui = tabPanel("OLS Residuals Map", verbatimTextOutput("moran"))
    )
  })
  
  observeEvent(input$removemoran, {
    removeUI(
      selector = "#moran"
    )
  })
  
  click_tract <- eventReactive(input$map_shape_click, {
    return(input$map_shape_click$id)
  })
  
  focal_tract <- reactive({
    req(click_tract())
    return(data()[data()$id == click_tract(), ])
  })
  
  neighbors <- reactive({
    req(click_tract())
    if (input$radio == 1) {
      return(data()[data()$id %in% aq()[[click_tract()]], ])
    } else if (input$radio == 2) {
      return(data()[data()$id %in% ar()[[click_tract()]], ])
    } else if (input$radio == 3) {
      v <- knn()[click_tract(), ]
      return(data()[data()$id %in% v, ])
    } else if (input$radio == 4) {
      v <- dist()[[click_tract()]]
      if (v == 0) {
        return(NULL)
      } else {
        return(data()[data()$id %in% v, ])
      }
    }
  })

  
  observe({
    req(click_tract())
    proxy <- leafletProxy('map')
    if (!is.null(neighbors())) {
      proxy %>%
        removeShape('focal') %>%
        clearGroup('neighbors') %>%
        addPolygons(data = neighbors(), fill = FALSE, color = '#FFFF00',
                    group = 'neighbors', opacity = 1) %>%
        addPolygons(data = focal_tract(), color = '#00FFFF', 
                    opacity = 1, layerId = 'focal', fillColor = 'transparent')
    } else {
      proxy %>%
        removeShape('focal') %>%
        clearGroup('neighbors') %>%
        addPolygons(data = focal_tract(), color = '#00FFFF', 
                    opacity = 1, layerId = 'focal', fillColor = 'transparent')
    }
  })
  
  # this reactive output contains the summary of the dataset and display the summary in table format
  output$filedf <- renderPlot({
    if(is.null(data())){return ()}
    plot(data())
  })
  
  # this reactive output contains the summary of the dataset and display the summary in table format
  output$sum <- renderPrint({
    if(is.null(data())){return ()}
    summary(data())
  })
  
  ## Selecting Variable Section
  observe({
    updateSelectInput(
      session,
      "indep",
      choices=names(data()))
  })
  
  observe({
    updateSelectInput(
      session,
      "depend",
      choices=names(data()))
  })
  
  observe({
    req(input$depend)
    updateSelectInput(session, "gwrdepend",
                      choices = input$depend
    )
  })
  
  formula <- reactive({
    paste(input$indep," ~ ",paste(input$depend,collapse="+"))

  })
  
  runRegression <- reactive({
    lm(as.formula(formula()),data=data())
    })
  
  output$reg2 <- renderUI({
    if(!is.null(input$indep)){
      data.ols2 <- sjt.lm(runRegression(), emph.p = TRUE)
      HTML(data.ols2$knitr)
    } else {
      print(data.frame(Warning="Please select Model Parameters."))
    }
  })
  
  observeEvent(input$add, {
      insertUI(
        selector = "#add",
        where = "afterEnd",
        ui = tabPanel("Regression", htmlOutput("reg2"))
      )
  })
  
  observeEvent(input$remove, {
    removeUI(
      selector = "#reg2"
    )
  })
  # This reactive output contains the dataset and display the dataset in table format
  output$table <- renderDataTable({
    if(is.null(data())){return ()}
    as.data.frame(data())
  })
  # This reactive output contains the dataset and display the dataset in table format
  output$map <- renderLeaflet({
    icons <- awesomeIcons(
      icon = 'ion-android-bicycle',
      library = 'ion'
    )
    icons2 <- awesomeIcons(
      icon = 'ion-android-car',
      markerColor = 'red',
      library = 'ion'
    )
    html_legend <- "<img src='http://leafletjs.com/docs/images/leaf-green.png'>green<br/>
    <img src='http://leafletjs.com/docs/images/leaf-red.png'>red"
    if(!is.null(data())){return (
      leaflet(data()) %>% addTiles() %>% 
        addPolygons(layerId = ~id,fillColor = 'transparent', 
                    color = 'blue', weight = 0.5, smoothFactor = 0.1)
    )}
    if(!is.null(data())&!is.null(data2())){return (
      if (class(st_geometry(data()))[1] == "sfc_POLYGON"){
        leaflet() %>% addTiles() %>% 
          addAwesomeMarkers(data = data2(),
                            icon = icons,
                            popup = data()$addressStr,
                            group ="data 1") %>%
          addPolygons(data = data()) %>%
          addControl(html = html_legend, position = "bottomleft")}
      else if (class(st_geometry(data()))[1] == "sfc_POINT"){
        leaflet() %>% addTiles() %>% 
          addAwesomeMarkers(data = data(),
                            icon = icons,
                            popup = data()$addressStr,
                            group ="data 1") %>%
          addMarkers(data = data()) %>% 
          addControl(html = html_legend, position = "bottomleft")}
    )
    }
    if(!is.null(data3())){return (
      leaflet() %>% addPolygons(data = data3())
      %>%  addTiles()  
    )}
  })
  
  atx2 <- reactive({
    as(data(), "Spatial")
  })
  tract_weights <- reactive({
    return(nb2listw(include.self(poly2nb(atx2()))))
  })
  
  gi_data <- reactive({
    g <- localG(atx2()$MEDHINC_CY, tract_weights())
    atx2()$g <- g
    return(atx2())
  })
  
  data3 <- eventReactive(input$runols, {
    data3 <- as(data(), "Spatial")
    data3.ols<-runRegression()
    data3@data$data3.ols.res<-resid(data3.ols)
    data3
  })
  
  colorpal <- reactive({
    colorBin(input$colors, data3()$data3.ols.res, bins=5)
  })
  olsmap <- eventReactive(input$runols, {
    req(input$indep)
    req(input$depend)
    pal <- colorpal()
    
    if(!is.null(data3())){return (
      leaflet(data3()) %>% addProviderTiles("CartoDB.Positron") %>% 
        addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
                    opacity = 1.0, fillOpacity = 0.5,
                    fillColor = ~pal(data3.ols.res),
                    highlightOptions = highlightOptions(color = "white", weight = 2,
                                                        bringToFront = TRUE)) %>% 
        addLegend("bottomright", pal = pal, values = ~data3.ols.res,
                  title = "OLS Residuals",
                  opacity = 1)
    )}
  })
  output$gimap <- renderLeaflet({
    olsmap()
  })
  
  observe({
    vars <- all.vars(parse(text=as.character(input$depend)))
    vars <- as.list(vars)
    updateSelectInput(session=session, inputId = "inVar", choices = vars)
  })
  
  etext <- eventReactive(input$rungwr,{
    input$inVar
  })
  
  print <- eventReactive(input$savemap,{
    saveWidget(map, "map.html")
  })
  
  output$eText <- renderText({
    paste0("You have just run a GWR analysis with the dependent variable ", etext())
  })
  
  output$olsText <- renderText({
    paste0("The map model is : ", formula())
  })
  
  wght <- reactive({
    input$gweight
  })
  
  gwrdata <- eventReactive(input$rungwr, {
    req(input$indep)
    req(input$depend)
    x <- as.character(input$inVar)
    y <- input$indep
    sp_shape <- as(data(), "Spatial")
    bwG <- gwr.sel(formula(), data=sp_shape, gweight = gwr.Gauss, method = input$method, verbose = FALSE)
    gwrG <- gwr(formula(),data=sp_shape, bandwidth = bwG, gweight = gwr.Gauss, hatmatrix = TRUE)
    xse <- paste0(x, "_se")
    sigTest95 <-  abs(gwrG$SDF[[x]]) -2 * gwrG$SDF[[xse]]
    sigTest99 <-  abs(gwrG$SDF[[x]]) -3 * gwrG$SDF[[xse]]
    sf_gwr <- st_as_sf(gwrG$SDF)
    sigTest95[sigTest95>0] <- 1
    sigTest95[sigTest95<=0] <- 0
    sigTest99[sigTest99>0] <- 1
    sigTest99[sigTest99<=0] <- 0
    sf_gwr$sigTest95 <- sigTest95
    sf_gwr$sigTest99 <- sigTest99   
    bins<-3
    sf_gwr<-mutate(sf_gwr,parBin=cut2(sf_gwr[[x]],g=bins,levels.mean = TRUE))
    sf_gwr<-mutate(sf_gwr,sigBin=cut2(sigTest95,g=bins,levels.mean = TRUE))
    bvColors=c("#e8e8e8","#dfb0d6","#be64ac","#ace4e4","#a5add3","#8c62aa","#5ac8c8","#5698b9","#3b4994")
    levels(sf_gwr$parBin)<-1:bins
    sf_gwr$sigBin <- sigTest95 + sigTest99 + 1
    sf_gwr <- mutate(sf_gwr,value=paste(parBin,'-',sigBin,sep=''))
    index <- c(1, 2, 3)
    values <- c("Under 95%", "95%", "99%")
    sf_gwr$signper <- values[match(sf_gwr$sigBin, index)]
    sf_gwr
  })

  
  colorpal2 <- reactive({
    binpal <- colorFactor(c("#e8e8e8","#dfb0d6","#be64ac","#ace4e4","#a5add3","#8c62aa","#5ac8c8","#5698b9","#3b4994"), gwrdata()$value)
  })
  
  ### Bivariate Legend
  output$plot <- renderPlot({
    bvColors=c("#e8e8e8","#dfb0d6","#be64ac","#ace4e4","#a5add3","#8c62aa","#5ac8c8","#5698b9","#3b4994")
    legendGoal=melt(matrix(1:9,nrow=3))
    lg<-ggplot(legendGoal, aes(Var2,Var1,fill = as.factor(value)))+ geom_tile()
    lg<- lg + scale_fill_manual(name="",values=bvColors)
    lg<-lg+theme(legend.position="none")
    lg<-lg + theme(axis.title.x=element_text(size=rel(1),color=bvColors[3])) + xlab(" Greater Parameter -->")
    lg<-lg + theme(axis.title.y=element_text(size=rel(1),color=bvColors[3])) + ylab("   Higher Significance -->")
    lg<-lg+theme(axis.text=element_blank())
    lg<-lg+theme(line=element_blank())
    lg
    
  })
  
  gwrmap1 <- eventReactive(input$rungwr, {
  req(input$indep)
  req(input$depend)
  x <- as.character(input$inVar)
  y <- input$indep
  mapdata <- gwrdata()
  pup <- mapdata[[x]]
  pal <- colorpal2()
  if(!is.null(gwrdata())){return (
    leaflet(gwrdata()) %>% addProviderTiles("CartoDB.Positron") %>% 
      addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
                  opacity = 0.1, fillOpacity = 0.75,
                  fillColor = ~pal(value),
                  popup = paste0("Parameter Estimate: ", pup,"<br>",
                                 "Significance: ", as.character(gwrdata()$signper),"<br>",
                                 "Bin: ", as.character(gwrdata()$value)),
                  highlightOptions = highlightOptions(color = "white", weight = 2,
                                                      bringToFront = TRUE)) %>% 
      addLegend("bottomright", pal = pal, values = ~value,
                title = "GWR Residuals Residuals",
                opacity = 1)
  )}
  })
  ### Leaflet Map of GWR Analysis
  output$gwrmap <- renderLeaflet({
    gwrmap1()
  })

  
  output$selection <- renderPrint(
    input$mychooser)
  
  observeEvent(input$button, {
    toggle("slider1")
    
  })
  
  ### Knitr - Print results
  ##cpal <- reactive({
    ##colorNumeric(input$colors1, data()[[vari()]])
  ##})
  
  vari <- reactive({
    vari <- as.character(input$variable1)
  })
  
  mymap <- eventReactive(input$loadMap, {
    x <- data()
    ##pal <- cpal()
    facts_popup <- paste0("<br><strong>Selected Variable: </strong>", 
                   x[[vari()]])

    mapview(x, zcol = vari(), Legend = TRUE)@map
    
  })
  
  output$varimap <- renderLeaflet({
    mymap()
  })
  
  
  ### Download Handlers for each Map
  
  output$downloadModel <- downloadHandler(
    filename = function(){
      paste0(vari(), "-variable-map.html")
      },
    content = function(file) {
      saveWidget(mymap(), file)
    }
  )
  
  output$downloadOLS <- downloadHandler(
    filename = function(){
      paste0(formula(), "-ols-map.html")
    },
    content = function(file) {
      saveWidget(olsmap(), file)
    }
  )
  
  output$downloadGWR <- downloadHandler(
    filename = function(){
      paste0(input$inVar, "-gwr-map.html")
    },
    content = function(file) {
      saveWidget(gwrmap1(), file)
    }
  )
  
  ### Render UI - Only display after a file is uploaded 
  output$tb <- renderUI({
    if(is.null(data()))
    {return()}
    else
      tabsetPanel(id ="tabs1", 
                  ###tabPanel("About file", plotOutput("filedf")),
                  ###tabPanel("Summary", verbatimTextOutput("sum")),
                  tabPanel("Data", dataTableOutput("table")),
                  ###tabPanel("Attribute", verbatimTextOutput("selection"),  
                  ###chooserInput("mychooser", "Available frobs", "Selected frobs",
                  ###names(data()), names(data2()), size = 10, multiple = TRUE),
                  ###h5("Move Variable to the right for use in Spatial Regression Model"),
                  ###actionButton("button", "Select Variable")),
                  tabPanel("Variable Map", 
                           leafletOutput("varimap", width = "1000px", height = "600px"),
                           absolutePanel(top=60, right=20, draggable = TRUE,
                           ##selectInput("colors1", "Color Scheme",
                                       ##rownames(subset(brewer.pal.info, category %in% c("seq", "div")))
                           ##),
                           selectInput("variable1", "Variable",
                                       names(data())
                           ),
                           actionButton("loadMap", "Load the map"),
                           downloadButton('downloadModel', 'Download Variable Map', class="dlButton")
                  )),
                  tabPanel("Nearest Neighbor Map",
                           leafletOutput("map", width = "1000px", height = "600px")),
                  tabPanel("Regression Selection", 
                           selectInput("indep", label = h5("Select Dependent"), "", multiple = TRUE),
                           selectInput("depend", label = h5("Select Explanatory"), "", multiple = TRUE),
                           actionButton("add", "Add Regression Summary"),
                           actionButton("remove", "Remove Regression Summary")
                  ),
                  tabPanel("OLS Residuals Map", 
                           leafletOutput("gimap", width = "1000px", height = "600px")
                           ),
                  tabPanel("GWRmap", 
                           leafletOutput("gwrmap", width = "1000px", height = "600px"),
                           verbatimTextOutput("eText")


                           ##selectInput("colors", "Color Scheme",
                                       ##rownames(subset(brewer.pal.info, category %in% c("seq", "div"))), selected = "YlOrRd"
                           ##),
                           )
      )
  })
  })

shinyApp(server = server, ui = ui)