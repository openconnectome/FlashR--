communitiesTab <- function()
{
   tabPanel("Communities",
      h2("Community structure in the graph."),
      tabsetPanel(
      tabPanel("Single Algorithm",
      checkboxInput(inputId='dendPlot',label='Plot Dendrogram',FALSE),
      bsTooltip(id='dendPlot',
                title="Plot a dendrogram of the communities. Not available for Infomap, Spinglass, Multilevel, Label Propation, Leading Eigenvalue t-SNE or Laplacian.",
                placement='top'),
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
#   tabPanel("Compare Two Algorithms",
#      checkboxInput(inputId='tanglegram',label="Tanglegram",value=TRUE),
#      conditionalPanel(
#         condition = "input.tanglegram==true",
#      checkboxInput(inputId='sort',label="Sort",value=TRUE)),
#      h2("Compare Fastgreedy to Edge Betweenness Communities"),
#      plotOutput("communityComp1", height="800px",width="500px"),
#      h2("Compare Fastgreedy to Walktrap Communities"),
#      plotOutput("communityComp2", height="800px",width="500px"),
#      h2("Compare Edge Betweenness to Walktrap Communities"),
#      plotOutput("communityComp3", height="800px",width="500px")
#   ),
   tabPanel("Compare All Algorithms",
      h2("Compare All Communities"),
      p(paste(
           "This computes all the community structures currently implemented.",
           "This will take a while to compute, and will show a heatmap",
           "of a comparison between the communities.")),
      p(paste(
           "This uses the variation of information index of Meila,",
           "Journal of Multivariate Analysis, 98, 2007, 873-895",
           "Columns are sorted by an equality distance -- it counts",
           "the number of times two columns contain the same community",
           "index for each row.")),
      p(paste(
         "To the right of a designator of the community algorithm",
         "is a number indicating the number of communities",
         "found by the given algorithm.",
         "Note that gray level between rows is basically meaningless.")),
      wellPanel(
         fluidRow(
            column(6,
               h3("Heatmap comparison"),
               plotOutput("communityCompM", height="800px",width="500px")
            ),
            column(6,
               h3("Community Names"),
               p("Laplac: Laplacian embedding + Mclust."),
               p("RDPG: RDPG embedding + Mclust."),
               p("t-SNE: t-SNE embedding + Mclust."),
               p("Fast: fast greedy."),
               p("Edge: edge betweenness."),
               p("Walk: walktrap."),
               p("LEigen: leading eigenvector."),
               p("LabelP: label progagation."),
               p("SpinG: spinglass."),
               p("MulitL: multilevel."),
               p("InfoM: infomap.")
            )
         )
      )
      )
   )
)
}
