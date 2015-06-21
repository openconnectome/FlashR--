makeSidebar <- function()
{
  sidebarPanel(
  p(paste("Select a graph to load.",
          "Currently only graphml formats have been tested.",
          "However, it should be able to process:",
          paste(paste(eval(formals(read.graph)$format),collapse=", "),
                '.',sep=''),
       "It can take a while to load larger graphs,",
       "particularly ones with lots of attributes.",
       "Zipped graphs (humans) do not seem to load correctly at this time.",
          collapse=" ")),
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
      min=1,max=1000),
            p("Clicking on `Plotting Parametes' opens/closes a menu of options for the plot."),
            plottingOpts()
  )
}
