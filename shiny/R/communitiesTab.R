communitiesTab <- function()
{
   tabPanel("Communities",
      tabsetPanel(
      tabPanel("Single Algorithm",
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
   ),
   tabPanel("Compare Two Algorithms",
      checkboxInput(inputId='tanglegram',label="Tanglegram",value=TRUE),
      conditionalPanel(
         condition = "input.tanglegram==true",
      checkboxInput(inputId='sort',label="Sort",value=TRUE)),
      h2("Compare Fastgreedy to Edge Betweenness Communities"),
      plotOutput("communityComp1", height="800px",width="500px"),
      h2("Compare Fastgreedy to Walktrap Communities"),
      plotOutput("communityComp2", height="800px",width="500px"),
      h2("Compare Edge Betweenness to Walktrap Communities"),
      plotOutput("communityComp3", height="800px",width="500px")
   ),
   tabPanel("Compare All Algorithms",
      h2("Compare All Communities"),
      plotOutput("communityCompM", height="800px",width="500px")
   )
   )
   )
}
