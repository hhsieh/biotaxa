library(shiny)

shinyUI(fluidPage(
  titlePanel(''),

  sidebarLayout(
    sidebarPanel(
      helpText("This app visualizes the accumulative curve of discovered species belong to same higher taxa.
               Two user-input parameters, 'taxa' and 'rank', are required.
               'taxa' refers to a specific taxa to which all inquired species belong to.
               'rank' refers to the rank that the user would like to visualize in the the accumulative curve.
               The default setting visualises the accumulative numbers of discovered phyla belong to taxa 'Animalia' over time.
               Please note that the rank of entered 'rank' needs to be lower than that of entered 'taxa'."),
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
