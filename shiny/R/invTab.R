invTab <- function()
{
   tabPanel("Graph Invariants and Statistics",
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
   )
}
