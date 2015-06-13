
shinyUI(pageWithSidebar(
  
  ###  Application title
  headerPanel("Graph Explorer"),
  
  ### Sidebar with text box for input search term
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
  ),
  ### Main Panel
  mainPanel(
    tabsetPanel(
      tabPanel("At/ributes",
       tabsetPanel(
            tabPanel("Graphical",
            wellPanel(
               fluidRow(
                  column(6,
                     h3("Vertex Attributes"),
                     selectInput(inputId="vertexAtts",
                                 label="Attribute to plot",
                     choices='None',selected='None',selectize=FALSE),
                     plotOutput("plotVA", height="250px")),
                  column(6,
                     h3("Edge Attributes"),
                     selectInput(inputId="edgeAtts",label="Attribute to plot",
                     choices='None',selected='None',selectize=FALSE),
                     plotOutput("plotEA", height="250px"))))
            ),
         tabPanel("Tabular",
      p("Tables of attributes. Note that for large graphs an error may be thrown due to the limitations of the datatable."),
     p("It can take a while to load larger graphs, particularly ones with lots of attributes."),
            wellPanel(
               fluidRow(
         h3("Graph Attributes"),
         DT::dataTableOutput("graphAtt")
         ),
         fluidRow(
         h3("Vertex Attributes"),
         DT::dataTableOutput("vertexAtt")
         ),
         fluidRow(
         h3("Edge Attributes"),
         DT::dataTableOutput("edgeAtt")
         )))
      )),
      tabPanel("Graph",
         p("Plotting the graph may take a long time unless the fast plotting option is checked."),
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
                      min=0.0005,max=.25,value=0.005,step=0.0005)
         ),
         bsTooltip(id='fast',
             title="Fast plotting only plots vertices and edges, no attributes or colors.",
                   placement='top'),
         tabsetPanel(
            tabPanel("2D",
            sliderInput(inputId='vertexSize',
                          label="Vertex Size",min=1,max=15,
                          value=1,step=1),
         bsTooltip(id='vertexSize',
                   title="The size of the dot representing the vertex.",
                   placement='top'),
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
                     webGLOutput("plotgraph3d", height="800px",width="800px"),
                     br())
         )
      ),
      tabPanel("Invariants",
          tabsetPanel(
            tabPanel("Tabular",
               bsCollapse(
                  bsCollapsePanel(title="Select Invariants",
                     checkboxGroupInput(inputId='invariantsList',
                         label="Invariants to Compute",
                         choices=invariants$Invariant,
                         selected=c("Order","Size","#Components",
                                    "Max Component Size")))
                  ),
                  DT::dataTableOutput("mytable",width='50%'),
                  br()
            ),
            tabPanel("Graphical",
               selectInput(inputId="invariants",label="Invariant to Plot",
                      choices=plot.invariants,selected="Degree Distribution"),
               plotOutput("plotinvariants", height="500px"),
         conditionalPanel(
              condition = "input.invariants == 'Alpha Centrality'",
                 numericInput(inputId="alpha",label="Alpha",
                      min=0,max=1,value=0.5,step=0.1)
         ),
         conditionalPanel(
              condition = "input.invariants == 'Page Rank'",
                 sliderInput(inputId="damping",label="Damping",
                      min=0,max=1,value=0.85,step=0.05)
         )
         )
         )
      ),
      tabPanel("Communities",
         checkboxInput(inputId='dendPlot',label='Plot Dendrogram',FALSE),
         bsTooltip(id='dendPlot',
                   title="Plot a dendrogram of the communities. Not available for Infomap, Spinglass, Multilevel, Label Propation, Leading Eigenvalue t-SNE or Laplacian.",
                   placement='top'),
         conditionalPanel(
              condition = "input.dendPlot == false",
         selectInput(inputId="CplotMethod",label="Plot Method",
               choices=plot.methods,selected="Auto",selectize=FALSE),
         sliderInput(inputId='CvertexSize',
                          label="Vertex Size",min=1,max=15,
                          value=5,step=1)
         ),
         selectInput(inputId="CvertexLabel",label="Vertex Label",
               choices="None",selected="None",selectize=FALSE),
         selectInput(inputId="communities",label="Communities to Compute",
                      choices=communities.list,selected="Fast Greedy"),
         conditionalPanel(
              condition = "input.communities == 'Laplacian' || input.communities == 'RDPG'",
         sliderInput(inputId='Cd',
                          label="Embedding Dimension",min=2,max=15,
                          value=5,step=1),
         sliderInput(inputId='CG',
                          label="Communities Range",min=2,max=50,
                          value=c(2,10),step=1)
         ),
         conditionalPanel(
            condition = "input.communities == 't-SNE' || input.CplotMethod == 't-SNE'",
               sliderInput(inputId='Ctheta',
                   label="Theta",value=0.5,min=0,max=1,step=0.1)
               
         ),
         plotOutput("communityPlot", height="1800px",width="1800px"),
         br(),br()
      )
   )
  )
))

