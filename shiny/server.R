### Main shinyServer function
shinyServer(function(input, output, session) {

set.seed(seed)

observe({
   cat("Open Connectome graph:",input$openconnectome,"\n")
})

observe({
   set.seed(input$seed)
})
   
  gGraph <- reactive({
     g <- NULL
     if(input$Source=='Local Disk'){
        if(!is.null(input$graphFile)){
           format <- "graphml"
           ex <- rev(strsplit(basename(input$graphFile$name),
                         split="\\.")[[1]])[[1]]
           cat("File Extension:",ex,"\n")
           if(ex %in% c("edgelist", "pajek", "ncol", "lgl",
                        "graphml", "dimacs", "graphdb", "gml", "dl")){
              format <- ex
            }
            cat("name:",input$graphFile$name,"\n")
            cat("datapath:",input$graphFile$datapath,"\n")
            if(grepl('\\.zip$',input$graphFile$name)){
              tf <- input$graphFile$datapath
              gr <- sub("\\.zip$","",input$graphFile$name)
              t1 <- system.time(g <- read.graph(unz(tf,gr),format=format))
            } else {
              t1 <- system.time(g <- read.graph(input$graphFile$datapath,
                    format=format))
            }
           cat("File read\n")
           print(t1)
        }
     } else {
        if(!is.null(input$openconnectome)){
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
              if(grepl("\\.zip$",file)){
                 progress <- Progress$new(session,min=1,max=10)
                 on.exit(progress$close())
                 progress$set(message = 'Download in progress',
                           detail='This may take a while...')
                 progress$set(value=1)
                 tf <- paste(tempfile(),"zip",sep='.')
                 t2 <- system.time(download.file(file,tf))
                 print(t2)
                 wd <- getwd()
                 setwd(tempdir())
                 gr <- sub("\\.zip$","",input$openconnectome)
                 progress$set(message = 'Unzipping graph')
                 cat("unzipping",tf,"extracting",gr,"\n")
                 unzip(tf,gr)
                 unlink(tf)
                 progress$set(value=3)
                 progress$set(message = 'Reading in graph',
                           detail='This too may take a while...')
                 cat("reading",gr,"\n")
                 t2 <- system.time(g <- read.graph(gr,format=format))
                 print(t2)
                 unlink(gr)
                 setwd(wd)
              } else {
                 g <- read.graph(file,format=format)
              }
           )
              cat("File read\n")
              print(t1)
        }
     }
     if(!is.null(g)){
        vertexLabels <- union("None",list.vertex.attributes(g))
        updateSelectInput(session,inputId="vertexLabel",
            choices=vertexLabels,selected="None")
        updateSelectInput(session,inputId="CvertexLabel",
            choices=c("None","Community",vertexLabels),selected="Community")
        updateSelectInput(session,inputId="vertexAtts",
            choices=vertexLabels,selected="None")
        edgeLabels <- union("None",list.edge.attributes(g))
        updateSelectInput(session,inputId="edgeLabel",
            choices=edgeLabels,selected="None")
        updateSelectInput(session,inputId="edgeColor",
            choices=edgeLabels,selected="None")
        updateSelectInput(session,inputId="edgeAtts",
            choices=edgeLabels,selected="None")
     }
     g
  })

  layout <- reactive({
     set.seed(input$seed)
     g <- gGraph()
     getLayout(g, 
               plotMethod=input$plotMethod, 
               u=input$u, 
               FRniter=input$FRniter,
               FRcoolexp=input$FRcoolexp, 
               circular=input$circular, 
               star.center=input$star.center,
               n=input$n, 
               KKniter=input$KKniter, 
               KKinittemp=input$KKinittemp, 
               KKcoolexp=input$KKcoolexp,
               scaleLaplacian=input$scaleLaplacian,
               dim=3,plotOnly=TRUE,theta=input$theta)
  })

  ## to allow for different layouts on the graph tab than on the
  ## communities tab. I've chosen to only do normalized Laplacian,
  ## and for the RDPG it is cbind(u,v) unless the graph is undirected,
  ## in which case it is u. Note that it doesn't scale by d, but that
  ## only changes the component variances, and this is being passed to
  ## mclust.
  layoutC <- reactive({
     set.seed(input$seed)
     g <- gGraph()
     getLayout(g, 
               plotMethod=input$CplotMethod, 
               u=ifelse(is.directed(g),"UV", "U"),
               FRniter=input$FRniter,
               FRcoolexp=input$FRcoolexp, 
               circular=input$circular, 
               star.center=input$star.center,
               n=input$n, 
               KKniter=input$KKniter, 
               KKinittemp=input$KKinittemp, 
               KKcoolexp=input$KKcoolexp,
               scaleLaplacian=TRUE,dim=input$Cd,plotOnly=FALSE,
               theta=input$Ctheta)
  })

  ### Generate plot output
  output$plotgraph <- renderPlot({  
     g <- gGraph()
     if(!is.null(g)){
        x <- layout()
        cat("Layout (plot):",dim(x),vcount(g),"\n")
        if(input$fast){
           fastPlot(g,x,input$UseAlpha,input$alphaLevel,input$vertexSize,
           color=2)
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
           if(input$UseAlpha) ec <- alpha(ec,input$alphaLevel)
           plot(g,vertex.size=input$vertexSize,vertex.label=vl,
                edge.width=weight,edge.label=el,
                edge.color=ec,
                layout=x)
        }
     }
  })

  output$plotgraph3d <- renderWebGL({  
     g <- gGraph()
     if(!is.null(g)){
        progress <- Progress$new(session,min=1,max=10)
        on.exit(progress$close())
        progress$set(message = 'Computing the 3d plot',
                     detail='Please be patient...')
        x <- layout()
        progress$set(value=1)
        fastPlot3D(g,x,input$UseAlpha,input$alphaLevel)
        progress$set(value=10)
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
         d <- computeInvariants(g,input$invariantsList)
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

  getCommunities <- reactive({
      set.seed(input$seed)
      g <- gGraph()
      if(is.null(g)) return(NULL)
      if(input$communities =="Fast Greedy"){
         z <- fastgreedy.community(as.undirected(simplify(g)))
      } else if(input$communities=="Leading Eigenvector"){
         z <- leading.eigenvector.community(as.undirected(simplify(g)))
      } else if(input$communities=="Label Propagation"){
         z <- label.propagation.community(g)
      } else if(input$communities=="Multilevel"){
         z <- multilevel.community(as.undirected(simplify(g)))
      } else if(input$communities=="Edge Betweenness"){
         z <- edge.betweenness.community(g)
      } else if(input$communities=="Infomap"){
         z <- infomap.community(g)
      } else if(input$communities=="Spinglass"){
         z <- spinglass.community(g)
      } else if(input$communities=="Walktrap"){
         z <- walktrap.community(g)
      } else if(input$communities=="Laplacian"){
         x <- computeLaplacian(g,d=input$Cd,normalize=TRUE)
         if(vcount(g)<=1000){
            z <- Mclust(x,G=input$CG[1]:input$CG[2])
         } else {
            init <- sample(vcount(g),1000)
            z <- Mclust(x,G=input$CG[1]:input$CG[2],
                    initialization=list(subset=init))
         }
      } else if(input$communities=="RDPG"){
         x <- graph.spectral.embedding(g,no=input$Cd)
         if(is.directed(g)) {
            y <- cbind(x$u,x$v)
         } else {
            y <- x$u
         }
         if(vcount(g)<=1000){
            z <- Mclust(y,G=input$CG[1]:input$CG[2])
         } else {
            init <- sample(vcount(g),1000)
            z <- Mclust(y,G=input$CG[1]:input$CG[2],
                    initialization=list(subset=init))
         }
      } else if(input$communities=="t-SNE"){
         x <- graph.spectral.embedding(g,no=25)
         if(is.directed(g)) {
            y <- cbind(x$u,x$v)
         } else {
            y <- x$u
         }
         x <- Rtsne(y,pca=FALSE,theta=input$Ctheta)$Y
         if(vcount(g)<=1000){
            z <- Mclust(x,G=input$CG[1]:input$CG[2])
         } else {
            init <- sample(vcount(g),1000)
            z <- Mclust(x,G=input$CG[1]:input$CG[2],
                    initialization=list(subset=init))
         }
      }
      z
  })

  ### Generate plot output
  output$communityPlot <- renderPlot({  
      dp <- input$dendPlot
      if(input$communities == "Infomap" ||
         input$communities == "Spinglass" ||
         input$communities == "RDPG" ||
         input$communities == "t-SNE" ||
         input$communities == "Multilevel" ||
         input$communities == "Label Propagation" ||
         input$communities == "Leading Eigenvalue" ||
         input$communities == "Laplacian" 
      ){
         updateCheckboxInput(session,inputId="dendPlot",value=FALSE)
         dp <- FALSE
         warning("Cannot dendPlot this community type\n")
      }
      g <- gGraph()
      if(is.null(g)) return(NULL)
      x <- layoutC()
      cat("Layout (community):",dim(x),"|V|:",vcount(g),"\n")
      z <- getCommunities()
      if(input$communities == "RDPG" ||
         input$communities == "t-SNE" ||
         input$communities == "Laplacian" 
      ){
         m <- z$classification
      } else {
         m <- membership(z)
      }
      vl <- input$CvertexLabel
      if(vl == 'None') {
         labels <- NULL
      } else if(vl == 'Community'){
         labels <- m
      } else {
         labels <- get.vertex.attribute(g,vl)
      }
      if(dp){
         dendPlot(z,labels=labels,main=paste(max(m),"Communities"))
      } else {
         ec <- input$edgeColor
         if(ec=='None') {
            ec <- 1
         } else {
            att <- get.edge.attribute(g,ec)
            uatt <- unique(att)
            ec <- ((match(att,uatt)-1) %% 8) + 1
         }
         col <- colors.list[((m-1) %% length(colors.list))+1]
         plot(g,layout=x[,1:2],edge.color=ec,vertex.color=col,
              vertex.size=input$CvertexSize,vertex.label=labels,
              main=paste(max(m),"Communities"))
      }
  })

})
