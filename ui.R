library(shiny)

shinyUI(fluidPage(
  titlePanel('Discovering Biodiversity in Antarctica and the Southern Ocean'),

  sidebarLayout(
    sidebarPanel(
      helpText("This app visualizes the accumulative curve of discovered species belong to same higher taxa.
               Two user-input parameters, 'taxa' and 'rank', are required.
               'taxa' refers to a specific taxa in a kingdom, a phylum, a class, a family or a genus to which all inquired species belong to.
               'rank' refers to the taxonomic rank below kingdom (i.e. phylum, class, family, genus or species) that the user would like to visualize in the accumulative curve.
               Please note that the rank of user-input 'rank' must be lower than that of user-input 'taxa'.
               The default setting exemplifies the accumulative numbers of discovered phyla belong to taxa 'Animalia' over time."),
      textInput("taxa", "enter taxa:","Animalia"),
      textInput("rank", "enter rank:","Phylum"),
      helpText("By selecting fitting curve, you chose to visualise the fitting curve and errors of a logistic regression model."),
      checkboxInput("model", label = "fitting curve (logistic regression)", value = FALSE)#,

      #actionButton("newplot", "Update View")
    ),
    mainPanel(
      tableOutput("dataview"),
      plotOutput("taxacurve"),
      plotOutput("fittingcurve")
    )
)
))
