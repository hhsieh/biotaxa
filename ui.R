library(shiny)

shinyUI(fluidPage(
  titlePanel(''),

  sidebarLayout(
    sidebarPanel(
      textInput("taxa", "enter taxa:","Animalia"),
      textInput("rank", "enter rank:","Phylum"),
      checkboxInput("model", label = "Fitting curve (logistic regression)", value = FALSE)#,

      #actionButton("newplot", "Update View")
    ),
    mainPanel(
      tableOutput("dataview"),
      plotOutput("taxacurve"),
      plotOutput("fittingcurve")
    )
)
))
