library(shiny)
library(SPIAT)
library(tidyverse)
library(ggsci)

options(repos = BiocManager::repositories())

#################################################################
##                         Sample data                         ##
#################################################################

images <- tibble(n = c(1:4),
                 filename = c("HGS_B_4_B+T_Scan1_Core[2,3,J]_[19600,46803]_component_data.tif",
                              "HGS_B_4_B+T_Scan1_Core[4,2,E]_[14543,54060]_component_data.tif",
                              "HGS_D_4_B+T_Scan1_Core[2,2,B]_[13421,44528]_component_data.tif",
                              "HGS_D_4_B+T_Scan1_Core[4,4,C]_[14067,54839]_component_data.tif"))
ihc_data <- read.csv("ihc_data.csv.gz")

#################################################################
##                        UI function                          ##
#################################################################

ui <- fluidPage(titlePanel("SPIAT package demo"), #Title 
                sidebarLayout(
                  #Sidebar Panel
  sidebarPanel(selectInput("image", "Select image", choices = setNames(images$filename, images$n)), #Image selection
               sliderInput("n_size", "Minimum neighborhood size", value = 12, min = 0, max = 30), #Neighborhood size input
               sliderInput("radius", "Maximum cell-cell distance", value = 15, min = 0, max = 30), #Radius input
               submitButton("Update View", icon("refresh"))), #Update view button
               #Main Panel
  mainPanel(imageOutput("selectedImage"), #Selected image output
            plotOutput("neighborhood_map", width = 400, height=400) # Image output
  )
                ))


#################################################################
##                       Server function                       ##
#################################################################

server <- function(input, output) {

  # Filter sample data based on which image is chosen
  single_image_data <- reactive({
    ihc_data %>% filter(Image == input$image) %>%
    mutate(ID = row_number())
    })


  # Create a spatial object based on the sample data from chosen image
  spe_object <- reactive({
    coord_x = single_image_data()$Centroid_x
    coord_y = single_image_data()$Centroid_y
    phenotypes = single_image_data()$Class
    dummy_intensity = rep(0, nrow(single_image_data()))
    intensity_matrix = matrix(dummy_intensity, nrow=1, ncol=nrow(single_image_data()))
    colnames(intensity_matrix) = single_image_data()$ID
    spe_object = format_image_to_spe(intensity_matrix=intensity_matrix,
                                   phenotypes = phenotypes, coord_x = coord_x,coord_y = coord_y)
  })

# Render image based on which image is chosen
output$selectedImage <- renderImage({
  list(
    src = file.path("images", paste0(input$image, ".png")),
    contentType = "image/jpeg",
    width = 400,
    height = 400
    )
  }, deleteFile = FALSE)

# Output neighborhood map based on image chosen, neighborhood radius and minimum neighborhood size

# Step 1 - use the SPIAT::identify_neighborhoods function to create a dataframe that tells me which cells are in which neighborhood
output$neighborhood_map <- renderPlot({
  plot1 <- identify_neighborhoods(spe_object(),
                                  cell_types_of_interest = c("CD20p", "CD3pCD8p", "CD3pCD8n", "CD79ApCD20n", "CD3pCD8nFoxP3p",
                                                             "CD3pCD8pFoxP3p"),
                                  method="hierarchical",
                                  radius=input$radius,
                                  min_neighborhood_size = input$n_size,
                                  feature_colname = "Phenotype")
  
  # Step 2 - extract the cell coordinates from the neighborhood output
  coordinates <- data.frame(spatialCoords(plot1)) %>%
    mutate(Cell.ID = row_number())
  
  # Step 3 - invert y-axis so that it matches the image (just a QuPath coordinates quirk)
  coordinates$Cell.Y.Position <- max(coordinates$Cell.X.Position) - coordinates$Cell.Y.Position

  # Step 4 - extract the phenotypes from the neighborhood output (this also has the neighborhood ID)
  phenotypes <- colData(plot1) %>%
    as.data.frame() %>%
    mutate(Cell.ID = as.numeric(Cell.ID)) %>%
    arrange(Cell.ID)

  # Join coordinates and neighborhood id
  coordinates <- left_join(coordinates, phenotypes)
  
  # I want to colour only the cells that are in a neighborhood, so I subset the previous df
  in_neighborhoods <- subset(coordinates, Neighborhood != "Free_cell" & !is.na(Neighborhood))

  # Calculate where to plot the label - this is part of the source code for identify_neighborhood function
  label_location <- vector()
  for (Clusternumber in unique(in_neighborhoods$Neighborhood)) {
    cells_in_Cluster <- in_neighborhoods[in_neighborhoods$Neighborhood ==
                                           Clusternumber, ]
    minX <- min(cells_in_Cluster$Cell.X.Position)
    maxX <- max(cells_in_Cluster$Cell.X.Position)
    minY <- min(cells_in_Cluster$Cell.Y.Position)
    maxY <- max(cells_in_Cluster$Cell.Y.Position)
    averageX <- (minX + maxX)/2
    averageY <- (minY + maxY)/2
    label_location <- rbind(label_location, c(Clusternumber,
                                              averageX, averageY))
  }

  label_location <- as.data.frame(label_location)
  
  # Clean up the label locations (convert to numeric, add column names)
  label_location[,2:3] <- sapply(label_location[,2:3], as.numeric)
  colnames(label_location) <- c("Cluster", "Cell.X.Position", "Cell.Y.Position")
  label_location$Cluster <- gsub("Cluster_", "", label_location$Cluster)

  #Make a colour vector (some colours are repeated but that is A-okay!)
  colors = rep(ggsci::pal_npg()(10), 100)[1:length(label_location$Cluster)]
  
  # Plot the cells, colored by neighborhood with text label with neighborhood ID
  ggplot(coordinates, aes(Cell.X.Position, Cell.Y.Position))+
    geom_point(alpha = 0.2)+
    geom_point(data=in_neighborhoods, aes(color=Neighborhood))+
    geom_text(data=label_location, aes(x=Cell.X.Position, y=Cell.Y.Position,
                                       label=Cluster), size=14)+
    scale_color_manual(values=colors)+
    theme_void()+
    theme(legend.position="none")
})

output$neighborhood_composition <- renderPlot({plot(1,1)})

}


shinyApp(ui = ui, server = server)
