attrTab <- function(){
   tabPanel("Attributes of the Graph",
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
))
}
