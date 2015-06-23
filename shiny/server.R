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
        updateSelectInput(session,inputId="vertexAttsSize",
            choices=vertexLabels,selected="None")
        updateSelectInput(session,inputId="vertexAttsColor",
            choices=vertexLabels,selected="None")
        updateSelectInput(session,inputId="CvertexLabel",
            choices=c("None","Community",vertexLabels),selected="Community")
        updateSelectInput(session,inputId="vertexAtts",
            choices=vertexLabels,selected="None")
        updateSelectInput(session,inputId="coordinates",
            choices=vertexLabels[-grep("None",vertexLabels)])
        edgeLabels <- union("None",list.edge.attributes(g))
        updateSelectInput(session,inputId="edgeLabel",
            choices=edgeLabels,selected="None")
        updateSelectInput(session,inputId="edgeColor",
            choices=edgeLabels,selected="None")
        updateSelectInput(session,inputId="edgeAtts",
            choices=edgeLabels,selected="None")
        m <- max(3,floor(vcount(g)/40))
        updateSliderInput(session,inputId="subsample",
            min=2,max=m,value=min(4,m-1),step=1)
     }
     g
  })

  ## For the RDPG method, the input$u field controls which singular
  ## vectors are plotted:
  ## if we are plotting in 2D, it concatenates u[,1] and v[,1].
  ## if in 3D, it concatenates u[,1:2] and v[,1].
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
               dim=as.numeric(input$layoutD),plotOnly=TRUE,theta=input$theta,
               coords=input$coordinates)
  })

  getSubsampled <- reactive({
     g <- gGraph()
     if(!is.null(g)){
        A <- get.adjacency(g)
        cat("Subsample:",input$subsample,"\n")
        if(input$subsample>1){
           if(input$contour || nrow(A)>2000){
              w <- input$subsample
              cat("Window:",w,"\n")
              B <- apply(A,1,slideFunct,window=w,step=round(w/2))
              A <- apply(B,2,slideFunct,window=w,step=round(w/2))
           }
        }
     }
     A
  })

  ### Generate plot output
  output$adjacencyPlot <- renderPlot({  
     g <- gGraph()
     if(!is.null(g)){
        A <- getSubsampled()
        if(input$contour){
           contour(1:nrow(A),1:ncol(A),A,xlab="",ylab="",axes=FALSE)
        } else {
           image(t(A)[nrow(A):1,])
        }
     }
  })

  ### Generate plot output
  output$plotgraph <- renderPlot({  
     g <- gGraph()
     if(!is.null(g)){
        x <- layout()
        if(is.null(x)) return(NULL)
        cat("Layout (plot):",dim(x),vcount(g),"\n")
        if(input$sizeByVar && input$vertexAttsSize != 'None'){
           size <- get.vertex.attribute(g,input$vertexAttsSize)
           if(all(is.numeric(size))){
              size <- 3*size/max(size)
           } else {
              warning("sizing vertices using a non-numeric label")
              a <- sort(unique(size))
              size <- match(size,a)
              size <- 3*size/max(size)
           }
        } else {
           size <- input$vertexSize
        }
        color <- "SkyBlue"
        if(input$colorByVar && input$vertexAttsColor != 'None'){
           color <- get.vertex.attribute(g,input$vertexAttsColor)
           vars <- color
           if(all(is.numeric(color))){
              if(diff(range(color))!=0) {
                 color <- gray((max(color)-color)/(max(color)-min(color)))
              }
           } else {
              a <- sort(unique(color))
              color <- match(color,a)
           }
           legnd <- unique(data.frame(name=vars,color=color,
                           stringsAsFactors=FALSE))
           legnd <- legnd[order(legnd$name),]
        } 
        if(input$fast){
           fastPlot(g,x,input$UseAlpha,input$alphaLevel,size,
           color=color)
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
           ec <- input$edgeColor
           if(ec=='None') {
              ec <- 1
           } else {
              att <- get.edge.attribute(g,ec)
              uatt <- unique(att)
              ec <- ((match(att,uatt)-1) %% 8) + 1
           }
           if(input$UseAlpha) ec <- alpha(ec,input$alphaLevel)
           plot(g,vertex.size=size,vertex.label=vl,
                vertex.color=color,
                edge.width=weight,edge.label=el,
                edge.color=ec,
                layout=x)
        }
        if(input$showLegend){
           n <- nrow(legnd)
           if(n>40){
              nc <- 5
           } else if(n>30) {
              nc <- 4
           } else if(n>20) {
              nc <- 3
           } else if(n>10) {
              nc <- 2
           } else {
              nc <- 1
           }
           legend(x="topright",legend=legnd$name,col=legnd$color,
                  ncol=nc,
                  pch=20)
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
        if(is.null(x)) return(NULL)
        progress$set(value=1)
        fastPlot3D(g,x,input$UseAlpha3D,input$alphaLevel3D,input$randomZ)
        progress$set(value=10)
     }
  })

  output$plotVA <- renderPlot({  
     g <- gGraph()
     if(!is.null(g)){
        if(input$vertexAtts != 'None'){
           a <- get.vertex.attribute(g,input$vertexAtts)
           if(class(a)=='numeric'){
              hist(a,xlab=input$vertexAtts,main="")
           } else {
              ta <- table(a)
              ta <- sort(ta,decreasing=TRUE)
              m <- min(30,length(ta))
              ta <- rev(ta[1:m])
              mar <- par('mar')
              par(mar=c(2,7,2,2))
              barplot(ta,horiz=TRUE,names=names(ta),xlab="",
                      las=2,cex.axis=.75)
              par(mar=mar)
           }
        } 
     } 
  })

  output$plotEA <- renderPlot({  
     g <- gGraph()
     if(!is.null(g)){
        if(input$edgeAtts != 'None'){
           a <- get.edge.attribute(g,input$edgeAtts)
           if(class(a)=='numeric'){
              hist(a,xlab=input$edgeAtts,main="")
           } else {
              ta <- table(a)
              ta <- sort(ta,decreasing=TRUE)
              m <- min(30,length(ta))
              ta <- rev(ta[1:m])
              mar <- par('mar')
              par(mar=c(2,7,2,2))
              barplot(ta,horiz=TRUE,names=names(ta),xlab="",
                      las=2,cex.axis=.75)
              par(mar=mar)
           }
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
           if(is.null(x)) return(NULL)
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
      if(is.null(att) || length(att)==0) return(NULL)
      inv <- data.frame(Attribute=c(att,"|V|","|E|"),
                Value=rep("",length(att)+2),
                stringsAsFactors=FALSE)
      inv[length(att)+1,'Value'] <- vcount(g)
      inv[length(att)+2,'Value'] <- ecount(g)
      for(i in 1:length(att)){
         inv[i,2] <- get.graph.attribute(g,att[i])
      }
      datatable(inv,rownames=FALSE)
  })

  output$vertexAtt <- DT::renderDataTable({
      g <- gGraph()
      if(is.null(g)) return(NULL)
      att <- list.vertex.attributes(g)
      if(is.null(att) || length(att)==0) return(NULL)
      inv <- data.frame(Attribute=att,
                Type=rep("numeric",length(att)),
                Values=rep("",length(att)),
                stringsAsFactors=FALSE)
      for(i in 1:length(att)){
         a <- get.vertex.attribute(g,att[i])
         inv[i,"Type"] <- class(a)
         if(class(a)=='numeric'){
            inv[i,"Values"] <- paste("range:",
                                 paste(round(range(a,na.rm=TRUE),3),
                                       collapse=" -- "))
         } else {
            inv[i,"Values"] <- paste("#unique:",length(unique(a)),collapse=" ")
         }
      }
      datatable(inv,rownames=FALSE)
  })

  output$edgeAtt <- DT::renderDataTable({
      g <- gGraph()
      if(is.null(g)) return(NULL)
      att <- list.edge.attributes(g)
      if(is.null(att) || length(att)==0) return(NULL)
      inv <- data.frame(Attribute=att,
                Type=rep("numeric",length(att)),
                Values=rep("",length(att)),
                stringsAsFactors=FALSE)
      for(i in 1:length(att)){
         a <- get.edge.attribute(g,att[i])
         inv[i,"Type"] <- class(a)
         if(class(a)=='numeric'){
            inv[i,"Values"] <- paste("range:",
                                 paste(round(range(a,na.rm=TRUE),3),
                                       collapse=" -- "))
         } else {
            inv[i,"Values"] <- paste("#unique:",length(unique(a)),collapse=" ")
         }
      }
      datatable(inv,rownames=FALSE)
  })

  getCommunities <- reactive({
                 progress <- Progress$new(session,min=1,max=10)
                 on.exit(progress$close())
                 progress$set(message = 'Computing communities',
                           detail='This may take a while...')
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
         x <- Rtsne(y,pca=FALSE,theta=input$Ctheta,dims=input$TSNEdimension)$Y
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
      x <- layout()
      if(is.null(x)) return(NULL)
      cat("Layout (community):",dim(x),"|V|:",vcount(g),"\n")
      z <- getCommunities()
      if(input$communities == "RDPG" ||
         input$communities == "t-SNE" ||
         input$communities == "Laplacian" 
      ){
         m <- z$classification
      } else {
         m <- membership(z)
         a <- order(as.numeric(names(m)))
         m <- m[a]
      }
      vl <- input$vertexLabel
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
              vertex.size=input$vertexSize,vertex.label=labels,
              main=paste(max(m),"Communities"))
      }
  })

  getCommunitiesMatrix <- reactive({
        g <- gGraph()
        progress <- Progress$new(session,min=1,max=9)
        on.exit(progress$close())
        progress$set(message = 'Computing all communities',
                  detail='This will take a while...')
        progress$set(value=1)
         x <- computeLaplacian(g,d=input$Cd,normalize=TRUE)
         if(vcount(g)<=1000){
            z <- Mclust(x,G=input$CG[1]:input$CG[2])
         } else {
            init <- sample(vcount(g),1000)
            z <- Mclust(x,G=input$CG[1]:input$CG[2],
                    initialization=list(subset=init))
         }
         lap <- z$classification
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
         rd <- z$classification
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
         ts <- z$classification
        M <- rbind(lap,rd,ts,
               membership(fastgreedy.community(as.undirected(simplify(g)))),
                   membership(edge.betweenness.community(g)),
                   membership(walktrap.community(g)),
                   membership(leading.eigenvector.community(
                        as.undirected(simplify(g)))),
                   membership(label.propagation.community(g)),
                   membership(spinglass.community(g)),
                   membership(multilevel.community(  
                      as.undirected(simplify(g)))),
                   membership(infomap.community(g)))
     a <- apply(M,1,max)
     rownames(M) <- paste(c("Laplac","RDPG","t-SNE",
                      "Fast","Edge","Walk",
                      "LEigen","LabelP","SpinG","MultiL","InfoM"),a)
     a <- hclust(eqDist(t(M)),method='ward.D2')
     M[,a$order]
  })

  output$communityCompM <- renderPlot({  
     heatmap(getCommunitiesMatrix(),labCol=NA,col=gray((255:0)/255),
             distfun=meila,Colv=NA)
  })

})
