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
        fileInput(inputId="graphFile",label="File",accept="graphml"),
  numericInput(inputId='seed',label="Random Number Seed",value=seed,
      min=1,max=100000,step=1),
            p("Clicking on `Plotting Parameters' opens/closes a menu of options for the plot."),
            plottingOpts(),
	 bsCollapse(
		 bsCollapsePanel("Save State",
           p(paste("Type in a base filename for the saved file.",
                   "This should be a pure root filename, no directory or",
                   "extensions.",
                   "If no filename is given, the default",
                   "is a file with a unique date/time in the name.")),
           textInput(inputId='outputfilename',label='Output Base Filename',
			           value=""),
               p(paste("The state file will be saved in the directory",
                       "set in your browser as the downloads directory",
                       "The file will be an RData file.")),
               downloadButton('downloadState','Download'),
		 id='downloadParameters'),id='downloadCollapse'),
	 bsCollapse(
		 bsCollapsePanel("Restore State",
               p(paste("Restoring the state restores all the variables",
                       "selected in the GUI with the exception of the",
                       "graph file.")),
              fileInput(inputId="saveFile",label="Restore File",
                        accept="RData"),
		 id='uploadParameters'),id='uploadCollapse')

  )
}
