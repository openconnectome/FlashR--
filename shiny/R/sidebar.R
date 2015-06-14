makeSidebar <- function()
{
  sidebarPanel(
  p("Select a graph to load. Currently only graphml formats are supported."),
  radioButtons(inputId="Source",label="Input Source",
      choices=c("Local Disk","Open Connectome"),
      selected="Local Disk"),
  conditionalPanel(
           condition = "input.Source == 'Local Disk'",
              fileInput(inputId="graphFile",label="File",accept="graphml")),
  conditionalPanel(
           condition = "input.Source == 'Open Connectome'",
  selectInput(inputId="openconnectome",
                                 label="Open Connectome Graph",
                     choices=openconnectome.graphs,
                     selected=openconnectome.graphs$worm[[1]],
                     selectize=TRUE)),
  numericInput(inputId='seed',label="Random Number Seed",value=seed,
      min=1,max=1000)
  )
}
