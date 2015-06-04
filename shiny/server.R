
computeInvariants <- function(g)
{
   cat("Computing Invariants:\n")
   inv <- data.frame(Invariants=invariants$Invariant,
                     Values=rep(0.0,nrow(invariants)),
                     Timings=rep(0.0,nrow(invariants)),
                     stringsAsFactors=FALSE)
   for(i in 1:nrow(invariants)){
      cat("Computing",invariants$Function[i],"\n")
      t1 <- system.time(inv[i,2] <- 
          round(do.call(invariants$Function[i],args=list(g=g))),3)
      inv[i,3] <- round(t1['elapsed'],3)
      cat("\t",inv[i,2],inv[i,3],"\n")
   }
   datatable(inv,rownames=FALSE)
}

fastPlot <- function(g,layout,alpha)
{
   plot(layout,pch=20,axes=FALSE,xlab="",ylab="")
   edges <- get.edgelist(g,names=FALSE)
   col <- alpha('black',alpha)
   if(is.directed(g)){
      arrows(layout[edges[,1],1],layout[edges[,1],2],
             layout[edges[,2],1],layout[edges[,2],2],length=0.1,col=col)
   } else {
      segments(layout[edges[,1],1],layout[edges[,1],2],
             layout[edges[,2],1],layout[edges[,2],2],col=col)
   }
   points(layout,pch=20,col=2)
}

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
  } else if(layout =='layout.coordinates'){
     x <- get.vertex.attribute(g,'x')
     if(!is.null(x)){
        y <- get.vertex.attribute(g,'y')
        layout <- cbind(x,y)
     } else {
        layout <- layout.auto(g,dim=2)
     }
  } else {
     layout <- get(layout)
  }
  layout

}

### Main shinyServer function
shinyServer(function(input, output, session) {

observe({
   cat("Open Connectome graph:",input$openconnectome,"\n")
})
   
  gGraph <- reactive({
     g <- NULL
     cat("Graph:",input$openconnectome,"\n")
     file <- openconnectome.dir
     for(i in 1:length(openconnectome.graphs)){
        if(any(openconnectome.graphs[[i]] == input$openconnectome)){
           file <- paste(file,names(openconnectome.graphs)[i],"/",
                         input$openconnectome,sep="")
           break
        }
     }
     cat("Getting",file,"\n")
     format <- "graphml"
     ex <- rev(strsplit(gsub(".zip","",input$openconnectome),
                        split="\\.")[[1]])[[1]]
     cat("File Extension:",ex,"\n")
     if(ex %in% c("edgelist", "pajek", "ncol", "lgl",
                  "graphml", "dimacs", "graphdb", "gml", "dl")){
        format <- ex
      }
     t1 <- system.time( 
        g <- read.graph(file,format=format))
        cat("File read\n")
        print(t1)
     vertexLabels <- union("None",list.vertex.attributes(g))
     updateSelectInput(session,inputId="vertexLabel",
         choices=vertexLabels,selected="None")
     updateSelectInput(session,inputId="vertexAtts",
         choices=vertexLabels,selected="None")
     edgeLabels <- union("None",list.edge.attributes(g))
     updateSelectInput(session,inputId="edgeLabel",
         choices=edgeLabels,selected="None")
     updateSelectInput(session,inputId="edgeColor",
         choices=edgeLabels,selected="None")
     updateSelectInput(session,inputId="edgeAtts",
         choices=edgeLabels,selected="None")
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
        if(input$fast){
           fastPlot(g,x,input$alpha)
        } else {
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
     }
  })

  output$plotVA <- renderPlot({  
     g <- gGraph()
     if(!is.null(g)){
        if(input$vertexAtts != 'None'){
           a <- get.vertex.attribute(g,input$vertexAtts)
           ta <- table(a)
           ta <- sort(ta,decreasing=TRUE)
           m <- min(30,length(ta))
           ta <- rev(ta[1:m])
           mar <- par('mar')
           par(mar=c(2,7,2,2))
           barplot(ta,horiz=TRUE,names=names(ta),xlab="",las=2,cex.axis=.75)
           par(mar=mar)
        }
     }
  })

  output$plotEA <- renderPlot({  
     g <- gGraph()
     if(!is.null(g)){
        if(input$edgeAtts != 'None'){
           a <- get.edge.attribute(g,input$edgeAtts)
           ta <- table(a)
           ta <- sort(ta,decreasing=TRUE)
           m <- min(30,length(ta))
           ta <- rev(ta[1:m])
           mar <- par('mar')
           par(mar=c(2,7,2,2))
           barplot(ta,horiz=TRUE,names=names(ta),xlab="",las=2,cex.axis=.75)
           par(mar=mar)
        }
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
        } else if(input$invariants=="Burt's Constraint"){
           plot(constraint(g),xlab="Vertex",
                ylab="Burt's Constraint",pch=20)
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
      d <- NULL
      if(!is.null(g)){
         d <- computeInvariants(g)
      }
      d
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
