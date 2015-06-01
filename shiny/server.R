

graph.spectral.embedding <- function(graph=g,no=2,
   cvec = degree(graph)/(vcount(graph) - 1))
{
	A <- get.adjacency(graph)
	B <- A+Diagonal(x=cvec)
	svds(B,k=no)
}

getLayout <- function(g,input)
{
  if(is.null(g)) return(NULL)
  layout <- paste('layout',gsub(" ",".",tolower(input$plotMethod)),
                  sep=".")
  if(layout == 'layout.auto'){
     layout <- layout.auto(g,dim=2)
  } else if(layout == 'layout.laplacian'){
     A <- graph.laplacian(as.undirected(g,mode='collapse'))
     z <- eigs(A,k=3,which="SM")
     d <- rev(z$values)
     layout <- z$vectors[,(length(d)-1):1]
  } else if(layout == 'layout.rdpg'){
     z <- graph.spectral.embedding(g)
     if(input$u=="U") {
        x <- z$u
     } else if(input$u=="V"){
        x <- z$v
     } else {
        x <- cbind(z$u[,1],z$v[,1])
     }
     layout <- x
  } else if(layout == 'layout.fruchterman.reingold'){
     layout <- layout.fruchterman.reingold(g,niter=input$FRniter,
                                          coolexp=input$FRcoolexp)
  } else if(layout == 'layout.fruchterman.reingold.grid'){
     layout <- layout.fruchterman.reingold.grid(g,niter=input$FRniter,
                                          coolexp=input$FRcoolexp)
  } else if(layout == 'layout.reingold.tilford'){
     layout <- layout.reingold.tilford(g,circular=input$circular)
  } else if(layout == 'layout.star'){
     layout <- layout.star(g,center=min(input$star.center,input$n))
  } else if(layout == 'layout.kamada.kawai'){
     layout <- layout.kamada.kawai(g,niter=input$KKniter,
                                          inittemp=input$KKinittemp,
                                          coolexp=input$KKcoolexp)
  } else if(input$plotMethod=='Coordinates'){
     if(input$linegraph == FALSE){
        layout <- getData(input)
     } else {
        layout <- layout.auto
     }
  } else {
     layout <- get(layout)
  }
  layout

}

### Main shinyServer function
shinyServer(function(input, output, session) {
   
  gGraph <- reactive({
     g <- NULL
     if(!is.null(input$graphFile)){
        g <- read.graph(input$graphFile$datapath,format="graphml")
        vertexLabels <- c("None",list.vertex.attributes(g))
        updateSelectInput(session,inputId="vertexLabel",
            choices=vertexLabels,selected="None")
        edgeLabels <- c("None",list.edge.attributes(g))
        updateSelectInput(session,inputId="edgeLabel",
            choices=edgeLabels,selected="None")
        updateSelectInput(session,inputId="edgeColor",
            choices=edgeLabels,selected="None")
     }
     g
  })

  layout <- reactive({
     g <- gGraph()
     getLayout(g,input)
  })

  ### Generate plot output
  output$plotgraph <- renderPlot({  
     g <- gGraph()
     if(!is.null(g)){
        x <- layout()
        vl <- input$vertexLabel
        if(vl=='None') vl <- NA
        else vl <- get.vertex.attribute(g,vl)
        el <- input$edgeLabel
        if(el=='None') el <- NA
        else el <- get.edge.attribute(g,el)
        weight <- 1
        if(input$useWeights) {
           if('weight' %in% list.edge.attributes(g)){
              weight <- get.edge.attribute(g,'weight')
           }
        }
        color <- 1
        ec <- input$edgeColor
        if(ec=='None') {
           ec <- 1
        } else {
           att <- get.edge.attribute(g,ec)
           uatt <- unique(att)
           ec <- ((match(att,uatt)-1) %% 8) + 1
        }
        plot(g,vertex.size=input$vertexSize,vertex.label=vl,
             edge.width=weight,edge.label=el,
             edge.color=ec,
             layout=x)
     }
  })

  output$plotinvariants <- renderPlot({  
     g <- gGraph()
     if(!is.null(g)){
        if(input$invariants=="Degree Distribution"){
           if(is.directed(g)){
              a <- degree.distribution(g,mode='out')
              b <- degree.distribution(g,mode='in')
              ylim <- range(c(a,b))
              plot(a,xlab="Degree",ylab="Proportion",log='x',
                   ylim=ylim,pch=20)
              points(b,pch=20,col=2)
              legend(0.9*length(a),max(ylim),legend=c("Out","In"),pch=20,col=1:2)
              
           } else {
              plot(degree.distribution(g),xlab="Degree",ylab="Proportion",
                   log='x',pch=20)
           }
        } else if(input$invariants=='Alpha Centrality'){
           plot(alpha.centrality(g,alpha=input$alpha),xlab="Vertex",
                ylab=expression(alpha*~centrality),pch=20)
        } else if(input$invariants=='Authority Score'){
           plot(authority.score(g)$vector,xlab="Vertex",
                ylab="Authority Score",pch=20)
        } else if(input$invariants=='Hub Score'){
           plot(hub.score(g)$vector,xlab="Vertex",
                ylab="Hub Score",pch=20)
        } else if(input$invariants=='Articulation Points'){
           a <- articulation.points(g)
           cols <- rep(1,vcount(g))
           cols[a] <- 2
           x <- layout()
           plot(g,layout=x,vertex.color=cols,
                     vertex.frame.color=cols,vertex.size=3)
           title("Articulation Points")
        } else if(input$invariants=='Boncich Power Centrality'){
           plot(bonpow(g),xlab="Vertex",
                ylab="Boncich Centrality",pch=20)
        } else if(input$invariants=='Betweenness'){
           plot(betweenness(g),xlab="Vertex",
                ylab="Betweenness",pch=20)
        } else if(input$invariants=='Eigenvalue Centrality'){
           plot(evcent(g)$vector,xlab="Vertex",
                ylab="Eigenvalue Centrality",pch=20)
        } else if(input$invariants=='Page Rank'){
           plot(page.rank(g,damping=input$damping)$vector,xlab="Vertex",
                ylab="Page Rank",pch=20)
        } else if(input$invariants=='Average KNN Degree'){
           plot(graph.knn(simplify(g))$knn,xlab="Vertex",
                ylab="Average KNN Degree",pch=20)
        } else if(input$invariants=='Closeness'){
           plot(closeness(g),xlab="Vertex",
                ylab="Closeness",pch=20)
        } else if(input$invariants=='Diversity'){
           plot(graph.diversity(g),xlab="Vertex",
                ylab="Diversity",pch=20)
        }
     }
  })


  output$mytable <- DT::renderDataTable({
      g <- gGraph()
      cl <- clusters(g)
      nc <- cl$no
      nc1 <- max(cl$csize)
      s <- ecount(g)
      inv <- data.frame(Invariant=c("Order","Size",
                                    "#Components","Max Component Size",
                                    "Diameter","Girth",
                                    "Density",
                                    "Average Path Length",
                                    "Centralization Degree",
                                    "Reciprocity"),
                      Value=c(vcount(g),ecount(g),
                              nc,nc1,
                              diameter(g),girth(g)$girth,
                              round(graph.density(g),3),
                              round(average.path.length(g),3),
                              round(centralization.degree(g)$centralization,3),
                              ifelse(is.directed(g),
                                     round(reciprocity(g),3),"NA")),
                      stringsAsFactors=FALSE)
      datatable(inv,rownames=FALSE)
  })

  output$graphAtt <- DT::renderDataTable({
      g <- gGraph()
      if(is.null(g)) return(NULL)
      att <- list.graph.attributes(g)
      if(is.null(att)) return(NULL)
      inv <- data.frame(Attribute=att,Value=rep("",length(att)),
                stringsAsFactors=FALSE)
      for(i in 1:length(att)){
         inv[i,2] <- get.graph.attribute(g,att[i])
      }
      datatable(inv,rownames=FALSE)
  })

  output$vertexAtt <- DT::renderDataTable({
      g <- gGraph()
      if(is.null(g)) return(NULL)
      att <- list.vertex.attributes(g)
      if(is.null(att)) return(NULL)
      inv <- as.data.frame(matrix("",nrow=vcount(g),ncol=length(att)),
                stringsAsFactors=FALSE)
      colnames(inv) <- att
      for(i in 1:length(att)){
         inv[,i] <- get.vertex.attribute(g,att[i])
      }
      datatable(inv,rownames=FALSE)
  })

  output$edgeAtt <- DT::renderDataTable({
      g <- gGraph()
      if(is.null(g)) return(NULL)
      att <- list.edge.attributes(g)
      if(is.null(att)) return(NULL)
      inv <- as.data.frame(matrix("",nrow=ecount(g),ncol=length(att)),
                stringsAsFactors=FALSE)
      colnames(inv) <- att
      for(i in 1:length(att)){
         inv[,i] <- get.edge.attribute(g,att[i])
      }
      datatable(inv,rownames=FALSE)
  })


})
