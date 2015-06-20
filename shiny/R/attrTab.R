attrTab <- function(){
   tabPanel("Attributes of the Graph",
    h2("The graph, edge and vertex attributes."),
    tabsetPanel(
         tabPanel("Graphical",
         h3("Bar plots of the attributes."),
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
         h3("Tables of the attributes."),
   p("Note that for large graphs an error may be thrown due to the limitations of the datatable function."),
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
))
}
