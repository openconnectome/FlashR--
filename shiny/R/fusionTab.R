fusionTab <- function()
{
   tabPanel("Fusion",
      h2("Fusion of graph and vertex attributes and/or covariates."),
      selectInput(inputId="fusion",label="Embedding Algorithm",
                   choices=c("RDPG","Laplacian"),selected="RDPG"),
      sliderInput(inputId='Cd',
                       label="Embedding Dimension",min=2,max=15,
                       value=5,step=1),
      sliderInput(inputId='CG',
                       label="#Mclust models",min=2,max=50,
                       value=c(2,10),step=1),
      plotOutput("fusionPlot", height="800px",width="800px"),
      br(),br()
   )
}
