graphplotTab <- function()
{
   tabPanel("Graph Plotting",
      tabsetPanel("Plot Vertices/Edges",
         tabPanel("2D",
      p("Plotting the graph may take a long time unless the fast plotting option is checked."),
            bsCollapse(
                  bsCollapsePanel(title="Plotting Parameters",
               selectInput(inputId="plotMethod",label="Plot Method",
                  choices=plot.methods,selected=plot.methods[1],selectize=FALSE),
               bsTooltip(id='plotMethod',
                         title="The type of layout for the graph.",
                         placement='top'),
               conditionalPanel(
                    condition = "input.plotMethod == 'Fruchterman Reingold' || input.plotMethod == 'Fruchterman Reingold Grid'",
                       selectInput(inputId="FRniter",label="Number of Iterations",
                            choices=c(10,100,500,1000),selected=500),
                       selectInput(inputId="FRcoolexp",label="Cooling exponent",
                            choices=c(0.1,0.5,0.99,1,3,5,10),selected=3)
               ),
               conditionalPanel(
                    condition = "input.plotMethod == 'Reingold Tilford'",
                         checkboxInput('circular',"Circular",FALSE)
               ),
               conditionalPanel(
                  condition = "input.plotMethod == 't-SNE'",
                     sliderInput(inputId='theta',
                         label="Theta",value=0.5,min=0,max=1,step=0.1)
                     
               ),
               conditionalPanel(
                    condition = "input.plotMethod == 'RDPG'",
                         radioButtons(inputId='u',label="Vectors",choices=c("U","V","UV"),
                                      selected="U")
               ),
               conditionalPanel(
                    condition = "input.plotMethod == 'star'",
                         numericInput(inputId='star.center',
                                      label="Center",min=1,max=1000,
                                      value=1)
               ),
               conditionalPanel(
                    condition = "input.plotMethod == 'Laplacian'",
                    checkboxInput(inputId="scaleLaplacian","Scaled Laplacian",
                                 TRUE)
               ),
               conditionalPanel(
                    condition = "input.plotMethod == 'Kamada Kawai'",
                       selectInput(inputId="KKniter",label="Number of Iterations",
                            choices=c(10,100,500,1000),selected=500),
                       selectInput(inputId="KKinittemp",label="Initial Temperature",
                            choices=c(0.1,0.5,0.99,1,3,5,10,20,50,100),selected=10),
                       selectInput(inputId="KKcoolexp",label="Cooling exponent",
                            choices=c(0.1,0.5,0.99,1,3,5,10),selected=0.99)
               ),
               checkboxInput('fast',"Fast Plotting",TRUE),
                checkboxInput('UseAlpha',"Use Alpha Blending",TRUE),
               conditionalPanel(
                    condition = "input.UseAlpha == true",
                    sliderInput(inputId="alphaLevel",label="Alpha Level",
                            min=0.0005,max=.25,value=0.1,step=0.0005)
               ),
               bsTooltip(id='fast',
                   title="Fast plotting only plots vertices and edges, no attributes or colors.",
                         placement='top'),
         sliderInput(inputId='vertexSize',
                       label="Vertex Size",min=1,max=15,
                       value=1,step=1),
      bsTooltip(id='vertexSize',
                title="The size of the dot representing the vertex.",
                placement='top')
            )),
      conditionalPanel(
           condition = "input.fast == false",
      selectInput(inputId="vertexLabel",label="Vertex Label",
         choices="None",selected="None",selectize=FALSE),
      selectInput(inputId="edgeLabel",label="Edge Label",
         choices="None",selected="None",selectize=FALSE),
      selectInput(inputId="edgeColor",label="Color Edges by",
         choices="None",selected="None",selectize=FALSE),
       checkboxInput('useWeights',"Weight Edges",FALSE)),
                  plotOutput("plotgraph", height="800px")),
         tabPanel("3D",
                checkboxInput('UseAlpha3D',"Use Alpha Blending",TRUE),
               conditionalPanel(
                    condition = "input.UseAlpha3D == true",
                    sliderInput(inputId="alphaLevel3D",label="Alpha Level",
                            min=0.0005,max=.25,value=0.1,step=0.0005)
               ),
                checkboxInput('randomZ',
                   label="Random Z (if layout doesn't provide z)",TRUE),
                  webGLOutput("plotgraph3d", height="800px",width="800px"),
                  br()),
         tabPanel("Adjacency Matrix",
            checkboxInput('contour',"Contour Plot",FALSE),
            bsTooltip(id='contour',
                title="A contour plot of a subsampled version of the adjacency matrix. Unlikely to be of much interest for small graphs.",
                placement='top'),
            conditionalPanel(
                 condition = "input.contour == true",
                     sliderInput(inputId='subsample',
                          label="Subsampling",min=2,max=3,
                          value=2,step=1)
            ),
            plotOutput("adjacencyPlot", height="600px",width="600px"),
            br())
      )
   )
}
