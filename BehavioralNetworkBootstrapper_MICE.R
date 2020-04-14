


#library(grid)
library(gridBezier)
#Incorporated into ProcessBOP_Sept4_2017.R on Sept 4, 2017


library(igraph)
library(asnipe)
library(assortnet)
library(vegan)
library(markovchain)
library(network)

#Plot matrix for quick vizualization
library(colorRamps)
colfunc<-colorRampPalette(rev(c("red","orange","yellow","springgreen","aquamarine","mediumblue","darkorchid4","gray40","black")))
vizframe<-function(matrix,camera){
 # if(camera=="A13"){mousecolors<-colorRampPalette(c("goldenrod3", "darkorange", "darkorchid3","blue3"))[]}
 # if(camera=="A24"){mousecolors<-c("goldenrod3", "darkorange", "darkgreen","darkred")}
  
  matrix[is.na(matrix)]<-0
  t3.d<-matrix
  t3.d.r<-t3.d[rev(1:nrow(t3.d)),]
  image(1:ncol(t3.d.r), 1:nrow(t3.d.r), t(t3.d.r), 
        breaks=c(0,1,20,50,90,300,500),
        col = c("black","yellow", "orange","blue3","gray40","darkorchid4"),
        zlim=c(0,200),
        xlim=c(1,640),ylim=c(1,480),
        axes = FALSE,asp=1, xlab='',ylab='',main="") 
}
 




rgb2hex <- function(r,g,b) rgb(r, g, b, maxColorValue = 255)
#################################
#
#datastandardizer takes quad-specific dataframes, and then creates standarized (1000 equally distributed measurements),
# and normalized (max = 1) dataframes for cumulative and rolling pee measures 
#output = list with 2 elements (cumulative pee, rolling pee{both standardized and normalized})
datastandardizer<-function(A1,A3,A2,A4,tnm,strain,p.day,thistrial,totalplot=TRUE){
      #smoothing
      #  http://r-statistics.co/Loess-Regression-With-R.html
    a1sex<-as.character(thistrial[which(thistrial$quad=="a1"),"sex"])
    a3sex<-as.character(thistrial[which(thistrial$quad=="a3"),"sex"])
    a2sex<-as.character(thistrial[which(thistrial$quad=="a2"),"sex"])
    a4sex<-as.character(thistrial[which(thistrial$quad=="a4"),"sex"])
  
  
  
  
    tlist<-list(A1,A3,A2,A4)
    anames<-c("a1","a3","a2","a4")
    
    ploflag<-1
    plotlist<-list()
    
    for(AAA in 1:length(tlist)){
      
      Aset<-tlist[[AAA]]
      wA<-anames[[AAA]]
      mousesex<-as.character(thistrial[which(thistrial$quad==wA),"sex"])
      
      thelabel<-paste(tnm,strain,mousesex,p.day,paste0(wA,"Normalized"),sep=".")
      
      if(nrow(Aset)>0){
        
          cumname<-paste0(wA,"cumpee")
          pee.vol.name<-paste0(wA,".pee.vol")
          
          ####CUMULATIVE PEE VOLUME
          num1<-approx(x=c(0,100),n=nrow(Aset))$y
          cumpeescaled.Aset<-cbind(num1,(Aset[,cumname]/sum(Aset[,pee.vol.name])))#scale to max 1
          # cumpeescaled.A1<-myFun(cumpeescaled.A1)
          cs5<-loess(cumpeescaled.Aset[,2]~cumpeescaled.Aset[,1], span=0.05) # 5% smoothing span
          sm5<-predict(cs5)
          sm5[sm5<0]<-0
          cumpeescaled.Aset[,2]<-sm5
          cumpeescaled<-cbind(approx(cumpeescaled.Aset,n=101)$x,approx(cumpeescaled.Aset,n=101)$y)
          
             
          ###########ROLL PEE
          roll.pee.name<-paste0(wA,".roll.pee")
          
          Aset[,roll.pee.name]<-rollapply(Aset[,pee.vol.name],width=1200,function(x) sum(x,na.rm=TRUE),partial=TRUE)
          rollingpeescaled.Aset<-Aset[,roll.pee.name]/sum(Aset[,pee.vol.name])#scale to area = 1
          rollingpeescaled.Aset<-cbind(num1,rollingpeescaled.Aset)
          loessMod5 <- loess(rollingpeescaled.Aset[,2]~rollingpeescaled.Aset[,1], span=0.05) # 2.5% smoothing span
          smoothed5 <- predict(loessMod5)
          smoothed5[smoothed5<0]<-0
          rollingpeescaled.Aset[,2]<-smoothed5        
          rollingpeescaled<-cbind(approx(rollingpeescaled.Aset,n=101)$x,approx(rollingpeescaled.Aset,n=101)$y)#reshape to 1000 points
          

          #################ROLL PEE COUNTS
          roll.peecount.name<-paste0(wA,".roll.peecount")
          
          didpee<-ifelse(Aset[,pee.vol.name]>0,1,0)
          
          Aset[,roll.peecount.name]<-rollapply(didpee,width=1200,function(x) sum(x,na.rm=TRUE),partial=TRUE)
          rollingpeecountscaled.Aset<-Aset[,roll.peecount.name]/sum(didpee)#scale to area = 1
          rollingpeecountscaled.Aset<-cbind(num1,rollingpeecountscaled.Aset)
          loessMod5 <- loess(rollingpeecountscaled.Aset[,2]~rollingpeecountscaled.Aset[,1], span=0.05) # 2.5% smoothing span
          smoothed5 <- predict(loessMod5)
          smoothed5[smoothed5<0]<-0
          rollingpeecountscaled.Aset[,2]<-smoothed5        
          rollingpeecountscaled<-cbind(approx(rollingpeecountscaled.Aset,n=101)$x,approx(rollingpeecountscaled.Aset,n=101)$y)#reshape to 1000 points


          if(totalplot==TRUE){
            

            plotlist[[ploflag]]<-list(cumpeescaled,rollingpeescaled,rollingpeecountscaled)
            names(plotlist)[ploflag]<-thelabel
            ploflag<-ploflag+1
          }
          
        ######################################################
      }else {
        cumpeescaled<-cbind(rep(NA,1000),rep(NA,1000))
        rollingpeescaled<-cbind(rep(NA,1000),rep(NA,1000))
        rollingpeecountscaled<-cbind(rep(NA,1000),rep(NA,1000))
      }
      
      if(AAA==1){
        standardized.cumulative<-cumpeescaled[,2]
        standardized.rolling<-rollingpeecountscaled[,2]
        standardized.rollcount<-rollingpeecountscaled[,2]
      } else {
        standardized.cumulative<-rbind(standardized.cumulative,cumpeescaled[,2])
        standardized.rolling<-rbind(standardized.rolling,rollingpeescaled[,2])
        standardized.rollcount<-rbind(standardized.rollcount,rollingpeecountscaled[,2])
      }

    }
  

    
    rownames(standardized.cumulative)<-c(paste(tnm,strain,a1sex,p.day,"a1NormalizedCumPee",sep="."),
                                          paste(tnm,strain,a3sex,p.day,"a3NormalizedCumPee",sep="."),
                                          paste(tnm,strain,a2sex,p.day,"a2NormalizedCumPee",sep="."),
                                          paste(tnm,strain,a4sex,p.day,"a4NormalizedCumPee",sep="."))
    
    rownames(standardized.rolling)<-c(paste(tnm,strain,a1sex,p.day,"a1NormalizedRolling",sep="."),
                                         paste(tnm,strain,a3sex,p.day,"a3NormalizedRolling",sep="."),
                                         paste(tnm,strain,a2sex,p.day,"a2NormalizedRolling",sep="."),
                                         paste(tnm,strain,a4sex,p.day,"a4NormalizedRolling",sep="."))
    
    rownames(standardized.rollcount)<-c(paste(tnm,strain,a1sex,p.day,"a1NormalizedRollCount",sep="."),
                                      paste(tnm,strain,a3sex,p.day,"a3NormalizedRollCount",sep="."),
                                      paste(tnm,strain,a2sex,p.day,"a2NormalizedRollCount",sep="."),
                                      paste(tnm,strain,a4sex,p.day,"a4NormalizedRollCount",sep="."))
    
    
    
    
    return(list(standardized.cumulative,standardized.rolling,standardized.rollcount,plotlist))
 
}


#peeOVERtimeplot
peeOvertimeplot<-function(listthingy,nameoflistthingy,boxt="l",linewidth=5,labelsize=10){
  
  cumpeescaled<-listthingy[[1]]
  rollingpeescaled<-listthingy[[2]]
  rollingpeecountscaled<-listthingy[[3]]
  
  plot(cumpeescaled,xlab="Time step",ylab="Normalized value",
       main=nameoflistthingy,ylim=c(0,1),bty=boxt,cex.main=labelsize)
  lines(cumpeescaled, col="blue",lwd=linewidth)

  #plot(rollingpeescaled)
  lines(rollingpeescaled,col='green',lwd=linewidth)

  #plot(rollingpeecountscaled)
  lines(rollingpeecountscaled,col='purple',lwd=linewidth,lty=2)

  legend(x="right", legend=c("Cumulative Vol", "Rolling Vol","Rolling Count"),
         col=c("blue", "green","purple"), lty=c(4,4,3), cex=0.8,
         box.lty=0,bg=rgb(1,1,1,alpha=0.1))
}


#########################################################################################
#MOUSE ARENA PLOTTER
PLOTARENA<-function(regionnames=c("none","full","short"),AL,lineweighting=10,
                    showme.looped,totalpees,plotsquarearena=FALSE,LABEL=individual.label){
  
  require(grDevices)
  require(colorRamps)
  require(geometry)
  require(plyr)
  require(reshape)
  ####################################################################################################################################
  library(MASS)
  library(plotrix)
  
  
  Arena<-read.table('C:/Users/Rusty/Amazon Drive/MICE/ArenaCoords.csv', sep = ",", header=TRUE)
  colnames(Arena)[1]<-"Region"
  Arena$Region<-factor(Arena$Region)
  
  labelslocations<-read.table('C:/Users/Rusty/Amazon Drive/MICE/centers.csv', sep = ",", header=FALSE);colnames(labelslocations)<-c("x","y")
  labelslocations$region<-c("SScorner","SSbarrier","centralcorner","OSbarrier","OScorner","OSwall",  "OSwater", "outercorner","SSwall",
                            "SSwater","center")
  labelslocations<-rbind(labelslocations,c(300,250,"generalmid"))
  labelslocations$shortnames<-c("ssc","ssb","cc","osb","osc","osw","osh","fc","ssw","ssh","c","gm")
  labelslocations<-data.frame(labelslocations)
  labelslocations[,c(1,2)]<-apply(labelslocations[,c(1,2)],2,function(x){as.numeric(as.character(x))})
  
  #xlim=c(0,600),ylim=c(100,700)
  
  if(plotsquarearena==TRUE){
    
    plot(Arena[,c("x","y")],xlim=c(0,700),ylim=c(750,50),xlab="",ylab="", axes=F,col="white")
    polygon(Arena[which(Arena$Region=='SSC'),c("x","y")],col=rgb(153,217,234,maxColorValue = 255), lty=1, lwd=1, border="black")
    par(new=T)
    polygon(Arena[which(Arena$Region=='SSB'),c("x","y")], col=rgb(112,146,190,maxColorValue = 255), lty=1, lwd=1, border="black")
    par(new=T)
    polygon(Arena[which(Arena$Region=='CC'),c("x","y")], col=rgb(181,230,29,maxColorValue = 255), lty=1, lwd=1, border="black")
    par(new=T)
    polygon(Arena[which(Arena$Region=='OSB'),c("x","y")], col=rgb(255,201,14,maxColorValue = 255), lty=1, lwd=1, border="black")
    par(new=T)
    polygon(Arena[which(Arena$Region=='OSC'),c("x","y")], col=rgb(34,177,76,maxColorValue = 255), lty=1, lwd=1, border="black")
    par(new=T)
    polygon(Arena[which(Arena$Region=='OSW'),c("x","y")], col=rgb(185,122,87,maxColorValue = 255), lty=1, lwd=1, border="black")
    par(new=T)
    polygon(Arena[which(Arena$Region=='OSH'),c("x","y")], col=rgb(237,28,36,maxColorValue = 255), lty=1, lwd=1, border="black")
    par(new=T)
    polygon(Arena[which(Arena$Region=='FC'),c("x","y")], col=rgb(127,127,127,maxColorValue = 255), lty=1, lwd=1, border="black")
    par(new=T)
    polygon(Arena[which(Arena$Region=='SSW'),c("x","y")], col=rgb(0,162,232,maxColorValue = 255), lty=1, lwd=1, border="black")
    par(new=T)
    polygon(Arena[which(Arena$Region=='SSH'),c("x","y")], col=rgb(239,228,176,maxColorValue = 255), lty=1, lwd=1, border="black")
    par(new=T)
    polygon(Arena[which(Arena$Region=='C'),c("x","y")], col=rgb(255,174,201,maxColorValue = 255), lty=1, lwd=1, border="black")
    par(new=T)
   
    
    polygon(rbind(c(14,121),c(14,697),c(585,697),c(585,121)), col=rgb(255,10,10,alpha=0.1,maxColorValue = 255), lty=1, lwd=1, border="black")
    par(new=T)
  } else {
    #plot(Arena[,c("x","y")],xlim=c(0,750),ylim=c(775,25),xlab="",ylab="", axes=F,col="white")
    #par(new=T)
  }
  #VERTEX LOCATIONS (originally 'captured' with:
  #
  #tkplot(showme.looped) #position verts manually
  #katchem<-tk_coords(tkp.id, norm = FALSE) #save coords
  #vertlocations<-dput(katchem)
  #
  
  
  vertlocations<-structure(c(387, 389, 137, 193, 249, 20, 33.4794520547945, 54, 
                             212, 385, 292, 20, 212, 5, 18, 75, 411, 389.669421487603, 91.4876033057851, 
                             200, 228, 396, 117, 0), .Dim = c(12L, 2L))
  
  rownames(vertlocations)<-c("OSbarrier","OScorner","OSwall","OSwater",
                             "SSbarrier","SScorner","SSwall","SSwater",
                             "center","centralcorner","generalmid","outercorner")
  
  flagme<-0
  for(jb in 1:nrow(vertlocations)){
    matcher<-grep(rownames(vertlocations)[jb],V(showme.looped)$name)
    if(length(matcher)==0){
      
    } else {
      
      keep<-vertlocations[jb,,drop=FALSE]
      if(flagme==0){
        flagme<-1
        newvert<-keep
      } else {
        newvert<-rbind(newvert,keep)
      }
    }
  }
  
  vertlocations<-newvert
  
    
  l <- layout_in_circle(showme.looped)
  
  colormat<-c(rgb2hex(255,201,14),rgb2hex(34,177,76),rgb2hex(185,122,87),rgb2hex(237,28,36),
    rgb2hex(112,146,190),rgb2hex(153,217,234),rgb2hex(0,162,232),rgb2hex(239,228,176),
    rgb2hex(255,174,201),rgb2hex(181,230,29),rgb2hex(150,150,150),rgb2hex(127,127,127))
  
  names(colormat)<-c("OSbarrier","OScorner","OSwall","OSwater",
                       "SSbarrier","SScorner","SSwall","SSwater",
                       "center","centralcorner","generalmid","outercorner")
  
  totalpees2<-totalpees[c("OSbarrier","OScorner","OSwall","OSwater",
                          "SSbarrier","SScorner","SSwall","SSwater",
                          "center","centralcorner","generalmid","outercorner")]
  
  tp2<-unlist(lapply(totalpees2,function(x){if(is.na(x)){bb<-0}else{bb<-x}}))
  names(tp2)<-c("OSbarrier","OScorner","OSwall","OSwater",
                "SSbarrier","SScorner","SSwall","SSwater",
                "center","centralcorner","generalmid","outercorner")
  totalpees2<-tp2
  

  totalpees2<-totalpees2[totalpees2>0]
  
  #set color
  V(showme.looped)$color<-colormat[V(showme.looped)$name]
  #Rather Let’s color the edges of the graph based on their source node color. 
  edge.start <- ends(showme.looped, es=E(showme.looped), names=F)[,1]
  edge.col <- V(showme.looped)$color[edge.start]
  
  
  
  plot.igraph(showme.looped,
              edge.width = E(showme.looped)$weight*100,
              # vertex.label=c("osb","osc","osw","osh",
              #                "ssb","ssc","ssw","ssh",
              #                "c","cc","gm","fc"),
              
              vertex.color = V(showme.looped)$color,
              
              edge.color=V(showme.looped)$color[edge.start],
              #edge.color='black',
              edge.curved = 0.2, 
              edge.arrow.size=.4,
              vertex.size=c(.31*totalpees2),
              edge.loop.angle=rep(5.7,length(E(showme.looped))),
              vertex.shape = "circle",
              label.col="black",
              vertex.frame.color="purple",
              layout=vertlocations,
              main=LABEL)
  
  par(new=T)
  polygon(rbind(c(-5,-5),c(-5,5),c(5,5),c(5,-5)), col=rgb(0,0,0,alpha=.1), lty=0, lwd=1, border="black")
  par(new=T)
  plot.igraph(showme.looped,
              edge.width = E(showme.looped)$weight*100,
              # vertex.label=c("osb","osc","osw","osh",
              #                "ssb","ssc","ssw","ssh",
              #                "c","cc","gm","fc"),
              
              vertex.color = rgb(0,0,0,alpha=0),
              
              edge.color=V(showme.looped)$color[edge.start],
              #edge.color='black',
              edge.curved = 0.2, 
              edge.arrow.size=.4,
              vertex.size=c(.31*totalpees2),
              edge.loop.angle=rep(5.7,length(E(showme.looped))),
              vertex.shape = "circle",
              label.col="black",
              vertex.frame.color="purple",
              layout=vertlocations,
              main=LABEL)
  
  # dropnamesslightly<-labelslocations[,c(1:2)]
  # dropnamesslightly$y<-dropnamesslightly$y-20
  # 
  # if(regionnames=="full"){
  #   text(dropnamesslightly,labels=labelslocations$region)
  # }
  # if(regionnames=="short"){
  #   text(dropnamesslightly,labels=labelslocations$shortnames)
  # }
  
}

#Return Center locations of ROIs
ROIcenters<-function(){
  
  Arena<-read.table('C:/Users/Rusty/Amazon Drive/MICE/ArenaCoords.csv', sep = ",", header=TRUE)
  colnames(Arena)[1]<-"Region"
  Arena$Region<-factor(Arena$Region)
  
  labelslocations<-read.table('C:/Users/Rusty/Amazon Drive/MICE/centers.csv', sep = ",", header=FALSE);colnames(labelslocations)<-c("x","y")
  labelslocations$region<-c("SScorner","SSbarrier","centralcorner","OSbarrier","OScorner","OSwall",  "OSwater", "outercorner","SSwall",
                            "SSwater","center")
  labelslocations<-rbind(labelslocations,c(300,250,"generalmid"))
  labelslocations$shortnames<-c("ssc","ssb","cc","osb","osc","osw","osh","fc","ssw","ssh","c","gm")
  labelslocations<-data.frame(labelslocations)
  labelslocations[,c(1,2)]<-apply(labelslocations[,c(1,2)],2,function(x){as.numeric(as.character(x))})
  

  return(labelslocations)
 
}


###############################
#runs chi-square tests relative to 'area.expected', 
#though a compatible df describing alt. expected values (e.g. time spent) can work too
chisquare.arena<-function(z,area.expected){
  summarytable<-(as.data.frame(table(z)))
  lookatboth<-merge(area.expected,summarytable,by.x="ROI",by.y="z",all.x=TRUE)
  
  lookatboth<-lookatboth[(match(area.expected$ROI,lookatboth$ROI)),]
  lookatboth[is.na(lookatboth)]<-0

  return(chisq.test(lookatboth$Freq,p=lookatboth$area.expected))
  
}

###############################
# creates a vector with the same length as nrow(compdataframe)
# containing numeric values indicating the duration, at every time point,
# SINCE the last of a given event (the timing of which is described with 'indices')
#and the first set uses 'avgtimesince' to populate before the first obs of that type
time.since.maker<-function(indices,compdataframe){
  #Could DO this if want to interpolate: function(indices,avgtimesince,compdataframe)
  
  for(ez in 1:length(indices)){
    #for first, just repeat NAs
    if(ez==1){
      timesincevalues<-c(rep(NA,(indices[ez]-1)))
      
      #This part of the code is only potentially useful if we're confident to extrapolate
      # missing 'time-since' data
      # if(length(timesincevalues)>round(avgtimesince)){
      #   
      #   
      # } else {
      #   #if the original number of NAs is less than the average (rounded) time between events
      #   #this will create a new vector of the appropriate length, but counting up an extrapolated
      #   #time from last event
      #   timesincevalues<-rev(seq(round(avgtimesince),(round(avgtimesince)-length(timesincevalues)+1),-1))
      # }
    } 
    #if this is the last one
    if(indices[ez]==last(indices)){
      thiswater<-indices[ez]
      lastwater<-nrow(compdataframe)
      
      #if this the last observation AND it occurs at the very end
      if(thiswater==lastwater){
        timesincevalues<-c(timesincevalues,0)
      } else {
        timesincevalues<-c(timesincevalues,0,c(seq(1,(lastwater-thiswater),1)))
      }
      
    } else {
      thiswater<-indices[ez]
      nextwater<-indices[(ez+1)]
      #if consecutive
      if((nextwater-thiswater)==1){
        timesincevalues<-c(timesincevalues,0)
      } else {
        timesincevalues<-c(timesincevalues,0,c(seq(1,(nextwater-thiswater-1),1)))
      }
      
    }
  }
  return(timesincevalues)
}

###############################
# creates a vector with the same length as nrow(compdataframe)
# containing numeric values indicating the distance from the last pee event, at every time point,
# FROM the last of a given event (the location and timing of which is described with 'peelocationframe')
#and the first set uses NA
area.dist.pee.funx<-function(mousequad,rawpeeloc,peelocationframe,compdataframe,Dcon){
  #Could DO this if want to interpolate: function(indices,avgtimesince,compdataframe)


  
  if(nrow(peelocationframe)!=nrow(compdataframe)){
    #cat("Nrows of peelocationframe and compdataframe do not match (and they need to)")
    
    justrightpee<-merge(peelocationframe,compdataframe[,c("dailysecond",paste0(mousequad,"x0"),paste0(mousequad,"y0"))],
                          by="dailysecond",all.y=TRUE)
      
    if(length(which(is.na(justrightpee[,paste0(mousequad,"x0")])))>0){
      justrightpee[which(is.na(justrightpee[,paste0(mousequad,"x0")])),c(paste0(mousequad,"x0"),paste0(mousequad,"y0"))]<-justrightpee[(which(is.na(justrightpee[,paste0(mousequad,"x0")]))-1),c(paste0(mousequad,"x0"),paste0(mousequad,"y0"))]
    }
    
  } else {
    
      justrightpee<-cbind(peelocationframe,compdataframe[,c(paste0(mousequad,"x0"),paste0(mousequad,"y0"))])
  }#close else
      
      #hold first obs before first pee
      holdfirstNAs<-justrightpee[which(justrightpee$dailysecond<rawpeeloc$dailysecond[1]),];
      holdfirstNAs$last.pee.true.x<-NA;holdfirstNAs$last.pee.true.y<-NA;holdfirstNAs$distfromlastpee<-NA
      
      #just use obs AFTER first pee event
      justrightpee<-justrightpee[which(justrightpee$dailysecond>=rawpeeloc$dailysecond[1]),]
      
      #Take care of Xs
      justrightpee$last.pee.true.x<-na.locf(justrightpee$true.x) #na.locf uses values from previous !is.na values to fill NA vals
      
      #this subs ensures that when we calculate the distance from the last pee, at the moment of a new pee, we use the LAST pee, rather than the current one
      justrightpee$last.pee.true.x[which(!is.na(justrightpee$true.x))]<-c(NA,justrightpee$last.pee.true.x[(which(!is.na(justrightpee$true.x))-1)])
      
      #Take care of Ys
      justrightpee$last.pee.true.y<-na.locf(justrightpee$true.y)
      #this subs ensures that when we calculate the distance from the last pee, at the moment of a new pee, we use the LAST pee, rather than the current one
      justrightpee$last.pee.true.y[which(!is.na(justrightpee$true.y))]<-c(NA,justrightpee$last.pee.true.y[(which(!is.na(justrightpee$true.y))-1)])
      
      
      
      justrightpee$distfromlastpee<-apply(justrightpee,1,function(x){dist(rbind(x[c(4,5)],x[c(6,7)]))})
      
      reallyjustright<-rbind(holdfirstNAs,justrightpee)
      
      distancefromlastpee<-reallyjustright$distfromlastpee/Dcon#cam-specific scale funxn
      
      compdataframe$dist.from.last.pee<-distancefromlastpee
      ###############

      #install.packages("adehabitatHR")
      library(adehabitatHR)
      library(sp)
      
      # 
      # 
      # RunningConvexPolygonArea<-function(xys){
      #  
      #   xys<-as.data.frame(xys)
      #   xys<-xys[-nrow(xys),]
      #   
      #   if(nrow(xys)<5){
      #     thisarea<-geometry::polyarea(x=xys[,1],y=xys[,2])
      #   } else {
      #     g<-SpatialPoints(xys, proj4string=CRS(as.character(NA)), bbox = NULL)
      #     
      #     heck<-mcp(g, percent = 100,unin = "m",unout = "m2")
      #     thisarea<-heck$area
      #   }
      #   return(thisarea)
      # }  
      
      
        x1<-na.omit(justrightpee[,c("dailysecond","true.x","true.y")])

        for(o in 1:nrow(x1)){
          
          if(o<11){
            sr=1
          } else {
            sr<-o-10
          }
          
          xys<-x1[c(sr:o),c(2:3)]
          
          if(nrow(xys)<5){
            thisarea<-geometry::polyarea(x=xys[,1],y=xys[,2])
          } else {
            g<-SpatialPoints(xys, proj4string=CRS(as.character(NA)), bbox = NULL)
            
            heck<-mcp(g, percent = 100,unin = "m",unout = "m2")
            thisarea<-heck$area
          }
          
          if(o<3){
            x1roll.A<-c(NA,NA)#can't calc. area for a point or line
          } else {
            x1roll.A<-c(x1roll.A,thisarea)
          }
          
          
        }
        
        x1rolls<-x1roll.A
        lookrollpolygon<-cbind(x1,x1rolls)
         
        # x1rolls<-rollapply(x1[,c(2:3)],width=11,FUN=function(x) {
        #   x<-as.data.frame(x);
        #   RunningConvexPolygonArea(x)},
        #   by.column=FALSE,
        #   fill=NA,
        #   align='right'
        #   )
        # 
        # 
        # lookrollpolygon<-cbind(x1,x1rolls)
        # 
        # for(z in 3:10){
        #   xys<-lookrollpolygon[c(1:z),c(1:2)]
        #   if(nrow(xys)<5){
        #     thisarea<-geometry::polyarea(x=xys[,1],y=xys[,2])
        #   } else {
        #     g<-SpatialPoints(xys, proj4string=CRS(as.character(NA)), bbox = NULL)
        #     
        #     heck<-mcp(g, percent = 100,unin = "m",unout = "m")
        #     thisarea<-heck$area
        #   }
        #   if(z==3){
        #     zz<-c(NA,NA,thisarea)
        #   } else {
        #     zz<-c(zz,thisarea)
        #   }
        # }
        # 
        # lookrollpolygon$x1rolls[c(1:10)]<-zz
        
        
       compdataframe2<-merge(reallyjustright,lookrollpolygon,by=c("dailysecond","true.y","true.x"),all=TRUE) 
       
       
       ##############
       #hold first obs before first pee
       holdfirstMCPNAS<-compdataframe2[which(compdataframe2$dailysecond<lookrollpolygon$dailysecond[3]),];
       holdfirstMCPNAS$x1rolls<-NA
       holdfirstMCPNAS$arealast10<-NA
       
       #just use obs AFTER first pee event
       compdataframe3<-compdataframe2[which(compdataframe2$dailysecond>=lookrollpolygon$dailysecond[3]),]
       
       
       #Take care of MIN CONVEX POLYGON PEE AREAS
       compdataframe3$arealast10<-na.locf(compdataframe3$x1rolls)
       #this subs ensures that when we calculate the distance from the last pee, at the moment of a new pee, we use the LAST pee, rather than the current one
       compdataframe3$arealast10[which(!is.na(compdataframe3$true.y))]<-c(NA,compdataframe3$arealast10[(which(!is.na(compdataframe3$true.y))-1)])
       
       
       
       
       megapee.polygons<-rbind(holdfirstMCPNAS,compdataframe3)

       
       #pix/cm
       megapee.polygons$arealast10<-((megapee.polygons$arealast10)/(Dcon^2))
        
       compdataframe$AreaLast10pees<-megapee.polygons$arealast10
      
      
      
      return(compdataframe)
    
  
}


#####################
# takes input kl dataframe, and corrected coordinate dataframe, and plugs in the corrected
# coordinates into the kl dataframe IF dimensions (nrow) match
coord.corrector<-function(in_kl,cor_frames){
  if(nrow(in_kl)!=nrow(cor_frames)){
    cat("Rows/obs for keeplong and correctedframes input DO NOT MATCH")
  } else {
      for(j in 1:ncol(cor_frames)){
        thiscol<-colnames(cor_frames)[j]
        in_kl[,thiscol]<-cor_frames[,thiscol]
      }
     return(in_kl)
  }
}




#takes dataframe (use peeroller), and trial/day specific roi information to add the roi where the mouse is
#RETURNS dataframe with same dimensions as input dataframe, but with ROI-identified locations for instances that were previously NA
peeroi.fixer<-function(dfwp,inforoi){
  library(sp)
  
  quaddos<-c("a1","a3","a2","a4")
  for(vero in 1:4){
    quad<-quaddos[vero]
    fixvec<-which(is.na(dfwp[,paste0(quad,"roi")]))
    quadroi<-dfwp[fixvec,]
    
    #pulls just the roi info for the correct camera
    if(quad=="a1"){inforoi2<-inforoi[which(inforoi$arena=="a1"),]}
    if(quad=="a3"){inforoi2<-inforoi[which(inforoi$arena=="a3"),]}
    if(quad=="a2"){inforoi2<-inforoi[which(inforoi$arena=="a2"),]}
    if(quad=="a4"){inforoi2<-inforoi[which(inforoi$arena=="a4"),]}
    
    #if no NAs, nothing to worry about!
    if(length(fixvec)==0){
      
    } else {
        errorflag<-1
        for(r in 1:length(fixvec)){
          ne<-0
          for(s in 1:nrow(inforoi2)){
            
            e<-point.in.polygon(quadroi[r,paste0(quad,"x0")],quadroi[r,paste0(quad,"y0")],pol.x=c(inforoi2[s,c(6,8,10,12)]),pol.y=c(inforoi2[s,c(7,9,11,13)]))
            #e2<-ifelse(e==1 & dfwp2[r,"quad"]==inforoi[s,"arena"],1,0)
            ne<-c(ne,e)
          }
          ne<-ne[-1] #drop first entry, which was created as 0 before s loop
          roispossible<-inforoi2$ROI_name[ne==1]
          roispossible<-unique(roispossible)
          fullpossiblerois<-paste(as.character(roispossible),collapse=".")
          if(length(roispossible)>1){#possible for a) waters, which are also on walls, and water is given precedence
            # and b) corners, and corners given precedence
            #print(as.character(roispossible))
            if(length(roispossible[grep("water",roispossible)])>0){
              roispossible<-roispossible[grep("water",roispossible)]
            }
            if(length(roispossible[grep("barrier",roispossible)])==2){
              roinam<-as.character(roispossible[1])
              roispossible<-paste0(strsplit(roinam,"_")[[1]][1],"_central_corner")
            }
            if(length(roispossible[grep("corner",roispossible)])>0){
              roispossible<-roispossible[grep("corner",roispossible)]
            }
            if(((length(roispossible[grep("barrier",roispossible)])>0)+(length(roispossible[grep("wall",roispossible)])))==2){
              roispossible<-roispossible[grep("barrier",roispossible)]
            }
            
          }
          
          #if multiple 'corner' rois possible, use the one that matches the quad identified by me, earlier in the code
          if(length(roispossible[grep("corner",roispossible)])>1){
            roispossible<-roispossible[grep(quad,roispossible)]
            #if the above trim doesn't work, length will still be >1
            #IN which case we take the lead quad
            if(length(roispossible)>1){
              dropindex<-which(roispossible==roispossible[grep(paste0("_",dfwp2[r,"quad"]),roispossible)])
              roispossible<-roispossible[-dropindex]
            }
          }
          
          if(length(roispossible)==0){
            roispossible<-"generalmid"
          }
          if(length(fullpossiblerois)==0){
            fullpossiblerois<-"NOMATCH"
          }
          
          roi.location<-as.character(roispossible)
          
          actual.location<-fixvec[r]
          
          dfwp[actual.location,paste0(quad,"roi")]<-roi.location
          
 
        }# close 'r' loopping through fixvec, for specific quad

    } #close 'else' statement activated if there ARE any NAs
    
  }

  return(dfwp)
}

#############################################
StandardizePeeTransitionMatrices<-function(testone){
  
  roi.names<-c("OSbarrier","OScorner","OSwall","OSwater",
    "SSbarrier","SScorner","SSwall","SSwater",
    "center","centralcorner","generalmid","outercorner")
  
  #find any rois not represented in roi matrix
  not.in.testone<-roi.names[-which(roi.names %in% rownames(testone))]
  
  #if there is at least 1 not represented, add it/them, and match var order
  if(length(not.in.testone)>0){
    for(nit in 1:length(not.in.testone)){
      testone<-cbind(testone,0)
      testone<-rbind(testone,0)
      last(colnames(testone))<-last(rownames(testone))<-not.in.testone[nit]
    }
    testone<-testone[roi.names,roi.names] #standardize order of rows/cols
    
  } else {
    testone<-testone[roi.names,roi.names]#standardize order of rows/cols
    
  }
  
  return(testone)
 
}


######################
# https://www.h2o.ai/blog/finally-you-can-plot-h2o-decision-trees-in-r/
# createDataTree(H2OTree) created that traverses a tree and translates it 
# from H2OTree into data.tree accumulating information about decision tree 
# splits and predictions into the node and edge attributes of a tree:
library(data.tree)

createDataTree <- function(h2oTree) {
  h2oTreeRoot = h2oTree@root_node
  dataTree = Node$new(h2oTreeRoot@split_feature)
  dataTree$type = 'split'
  addChildren(dataTree, h2oTreeRoot)
  return(dataTree)
}

addChildren <- function(dtree, node) {
  
  if(class(node)[1] != 'H2OSplitNode') return(TRUE)
  
  feature = node@split_feature
  id = node@id
  na_direction = node@na_direction
  
  if(is.na(node@threshold)) {
    leftEdgeLabel = printValues(node@left_levels, 
                                na_direction=='LEFT', 4)
    rightEdgeLabel = printValues(node@right_levels, 
                                 na_direction=='RIGHT', 4)
  }else {
    leftEdgeLabel = paste("<", node@threshold, 
                          ifelse(na_direction=='LEFT',',NA',''))
    rightEdgeLabel = paste(">=", node@threshold, 
                           ifelse(na_direction=='RIGHT',',NA',''))
  }
  
  left_node = node@left_child
  right_node = node@right_child
  
  if(class(left_node)[[1]] == 'H2OLeafNode')
    leftLabel = paste("prediction:", left_node@prediction)
  else
    leftLabel = left_node@split_feature
  
  if(class(right_node)[[1]] == 'H2OLeafNode')
    rightLabel = paste("prediction:", right_node@prediction)
  else
    rightLabel = right_node@split_feature
  
  if(leftLabel == rightLabel) {
    leftLabel = paste(leftLabel, "(L)")
    rightLabel = paste(rightLabel, "(R)")
  }
  
  dtreeLeft = dtree$AddChild(leftLabel)
  dtreeLeft$edgeLabel = leftEdgeLabel
  dtreeLeft$type = ifelse(class(left_node)[1] == 'H2OSplitNode', 'split', 'leaf')
  
  dtreeRight = dtree$AddChild(rightLabel)
  dtreeRight$edgeLabel = rightEdgeLabel
  dtreeRight$type = ifelse(class(right_node)[1] == 'H2OSplitNode', 'split', 'leaf')
  
  addChildren(dtreeLeft, left_node)
  addChildren(dtreeRight, right_node)
  
  return(FALSE)
}

printValues <- function(values, is_na_direction, n=4) {
  l = length(values)
  if(l == 0)
    value_string = ifelse(is_na_direction, "NA", "")
  else
    value_string = paste0(paste0(values[1:min(n,l)], collapse = ', '),
                          ifelse(l > n, ",...", ""),
                          ifelse(is_na_direction, ", NA", ""))
  return(value_string)
}



################################
prop.table.excludeNAs<- function(x) {x/sum(x, na.rm=TRUE)}

Network.Quick<-function(z,plotpretty=c('yes','arena','no'),
                        individual.label,regionnames=c("none","full","short")){
  # regionnames refers to options for ArenaPlotting
  summarytable<-table(z)
  totalpees<-c(summarytable)
  proportiontable<-prop.table(summarytable)
  uniquebehaviors<-length(summarytable)
  
  #TRUE SHANNON DIVERSITY
  dftable<-t(as.data.frame(summarytable))
  colnames(dftable)<-dftable[1,]
  dftable<-dftable[-1,,drop=FALSE]
  rownames(dftable) <- c()
  dftable<-data.frame(dftable)
  indx<-c(1:uniquebehaviors)
  dftable[indx]<- lapply(dftable[indx], function(x) as.numeric(as.character(x)))
  Shannon<-vegan::diversity(dftable,index = "shannon")
  trueShannon<-exp(Shannon)
  
  
  #TRANSITIONSSSSSSSS  
  ##########################################################################################################
  c.no.OFFs<-as.character(z)
  
  #uses "createSequenceMatrix" function (from *markovchain*) to calculate transition matrix
  TransitionMatrix.cum<-createSequenceMatrix(c.no.OFFs,sanitize=FALSE)
  
  
  noselfs<-TransitionMatrix.cum
  diag(noselfs)<-NA
  
  prop.table.excludeNAs<- function(x) {x/sum(x, na.rm=TRUE)}
  
  proportiontable.trans<-prop.table.excludeNAs(TransitionMatrix.cum) 
  
  ###########TRUE SHANNON DIVERSITY FOR TRANSITIONS
  if(uniquebehaviors>1){
    val3.s<-matrix(nrow=(((uniquebehaviors^2)-uniquebehaviors)/2),ncol=2)
    rownum<-1
    for (u in 1:uniquebehaviors){
      for(v in u:uniquebehaviors){#iteratively reduces columns analyzed in next loop to avoid double counting transitions 
        if (u!=v){
          val3.s[rownum,2]<-as.numeric(TransitionMatrix.cum[u,v]+TransitionMatrix.cum[v,u])
          val3.s[rownum,1]<-paste(row.names(TransitionMatrix.cum)[u],colnames(TransitionMatrix.cum)[v],sep="-")
          rownum<-rownum+1
        }
      }
    }
    val4.s<-as.data.frame(val3.s)
    val4.s[,2]<-as.numeric(as.character(val4.s[,2])) 
    
    dftable.T<-t(val4.s)
    colnames(dftable.T)<-dftable.T[1,]
    dftable.T<-dftable.T[-1,,drop=FALSE]
    rownames(dftable.T) <- c()
    dftable.T<-data.frame(dftable.T)
    indx<-c(1:length(dftable.T))
    dftable.T[indx]<- lapply(dftable.T[indx], function(x) as.numeric(as.character(x)))
    Shannon.T<-diversity(dftable.T,index = "shannon")
    trueShannon.T<-exp(Shannon.T)
    propShannon.T<-Shannon.T/log(specnumber(dftable.T))
    redundancy.T<-1-propShannon.T
    uniquetransitions<-length(val3.s)
    
  } else {
    trueShannon.T<-NA
    propShannon.T<-NA
    redundancy.T<-NA
    uniquetransitions<-1
  }
  
  
  #########################################################
  PropSelfTransitioning<-sum(diag(TransitionMatrix.cum))/sum(TransitionMatrix.cum)
  
  
  network.TransitionMatrix.cum<-network(proportiontable.trans,directed=TRUE,loops=TRUE,hyper=FALSE,matrix.type="adjacency")
  
  
  #Calculate behavioral group membership network
  community.observed <- fastgreedy.community(graph.adjacency(proportiontable.trans,mode="undirected",weighted=TRUE))
  ncommunities<-length(unique(community.observed$membership))
  
  # NETWORK SUMMARY VARIABLES
  showme.looped<- graph_from_adjacency_matrix(proportiontable.trans,mode="directed",weighted=TRUE,diag=TRUE)
  
  
  selfprops<-diag(proportiontable.trans)
  centersfornetworks<-ROIcenters()
  adjacencyList <- reshape2::melt(TransitionMatrix.cum)  # Convert to list of ties only
  adjacencyList <- adjacencyList[adjacencyList$value > 0, ] # remove zeros along diagonal
  adjacencyList[,c(1,2)] <- apply(adjacencyList[,c(1,2)],2,function(x){as.character(x)})
  
  maxweight<-max(adjacencyList$value)
  E(showme.looped)$width <- E(showme.looped)$weight
  
  
  
  showme.looped2<- graph_from_adjacency_matrix(proportiontable.trans,mode="undirected",weighted=TRUE,diag=TRUE)
  
  
  
  
  if(plotpretty=='yes'){
    l <- layout_with_fr(showme.looped)
    
    
  #set color
  V(showme.looped)$color<-roicolors[V(showme.looped)$name]
  
  #set node size
  sizes<-colSums(TransitionMatrix.cum)
  V(showme.looped)$size <- sizes*2
  
  #set edge width
  E(showme.looped)$width <- E(showme.looped)$weight *20
  
  #change arrow size and edge color:
  E(showme.looped)$arrow.size <- .6
  E(showme.looped)$edge.color <- "gray80"
  
  
  #Rather Let’s color the edges of the graph based on their source node color. 
  edge.start <- ends(showme.looped, es=E(showme.looped), names=F)[,1]
  edge.col <- V(showme.looped)$color[edge.start]
  
  plot(showme.looped,vertex.color=roicolors[V(showme.looped)$name],edge.curved=.2 ,edge.color=edge.col,layout=l)
  legend("bottomleft", names(sizes), pch=21,
         col="#777777", pt.bg=roicolors[V(showme.looped)$name], pt.cex=2, cex=.8, bty="n", ncol=1)
  title(individual.label)
  
  
  }
  
  if(plotpretty=='arena'){
    PLOTARENA(regionnames="short",AL=adjacencyList,lineweighting=10,showme.looped,totalpees,plotsquarearena=FALSE,LABEL=individual.label)
  }
  ########################################################
  
  # Klein DJ, Randić M. Resistance distance. Journal of Mathematical Chemistry 1993; 12: 81–95.
  # Ellens W, Kooij RE. 2013. Graph measures and network robustness. arXiv:1311.5064v1
  
  # Klein and Randić [30] found that the effective graph resistance of a connected network can be written as a function of all 
  # non-zero Laplacian eigenvalues of the network.
  
  # Yang et al. 2016. The Rationality of Four Metrics of Network Robustness: A Viewpoint of Robust Growth of 
  # Generalized Meshes. PLOS ONE. https://doi.org/10.1371/journal.pone.0161077
  
  lap.mat<-laplacian_matrix(showme.looped2, norm=FALSE, sparse=FALSE)
  ev.lm<-eigen(lap.mat)
  ev.lm<-unlist(ev.lm[1])
  sorted.Laplacian.eigenvalues<-sort(ev.lm)
  value.ER<-0
  for(lam in 2:length(sorted.Laplacian.eigenvalues)){
    value.ER<-value.ER+1/sorted.Laplacian.eigenvalues[lam]
  }
  
  EffectiveGraphResistance<-value.ER*length(sorted.Laplacian.eigenvalues)
  
  
  
  
  
  ############################################################################
  #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
  #6.1 Density
  
  #The proportion of present edges from all possible edges in the network.
  #Edge density: % of edges compared to maximum
  #  densityscore<-ecount(showme.looped)/(vcount(showme.looped)*(vcount(showme.looped)-1)) #for a directed network
  densityscore<-ecount(showme.looped)/(vcount(showme.looped)^2)#for a directed network WITH self-loops
  
  
  #
  deg <- igraph::degree(showme.looped,loops=TRUE)
  #Average degree: Avg. number of links
  mean.degree<-mean(deg)
  
  #MAverage path length: avg of shortest pathes between reachable nodes
  mean_path_length<-mean_distance(showme.looped, directed=T)
  
  #network diameter is the longest geodesic distance (length of the shortest path between two nodes) in the network.
  network.diameter<-diameter(showme.looped, directed=T,weights=NA)
  
  #average clustering coefficient
  avg.cluster.coeff<-transitivity(showme.looped)
  
  #Modularity
  Modularity.Value<-modularity(community.observed)
  
  
  
  ############################################
  #SMALL WORLDNESS --- video-wise, no filter
  #number of nodes/vertices in graph
  vertices<- uniquebehaviors
  #number of edges in G(n,m) graph
  edges<- sum(TransitionMatrix.cum!=0)
  
  rando.network<-sample_gnm(vertices, edges, directed = TRUE, loops = TRUE)
  Trobserved<-avg.cluster.coeff
  mean.Trrandom<-transitivity(rando.network)
  SPobserved<-mean_path_length
  mean.SPrandom<-mean_distance(rando.network, directed=T)  
  
  Smallworldness<- (Trobserved/mean.Trrandom)/(SPobserved/mean.SPrandom)
  ############################################
  
  orig.printme<-data.frame(ncommunities,uniquebehaviors,uniquetransitions,
                           Shannon,Shannon.T,
                           trueShannon,trueShannon.T,
                           EffectiveGraphResistance,
                           densityscore,mean.degree,mean_path_length,network.diameter,avg.cluster.coeff,
                           PropSelfTransitioning,Smallworldness,Modularity.Value)
  
  
  return(list(orig.printme,community.observed))
}

#BOOTSTRAP.networksummarizer performs network summary, with bootstrap iterations
# to test robustness.
# Returns 'RemovalSummary' DF which gives a measure of 'robustness'
# (EffectiveGraphResistance)
############################################################
# Effective graph resistance measures and concept based on:
# Klein DJ, Randić M. Resistance distance. Journal of Mathematical Chemistry 1993; 12: 81–95.
# Ellens W, Kooij RE. 2013. Graph measures and network robustness. arXiv:1311.5064v1

# "Klein and Randić [30] found that the effective graph resistance of a connected network can be written as a function of all 
# non-zero Laplacian eigenvalues of the network." --- from Yang et al. 2016

# Yang et al. 2016. The Rationality of Four Metrics of Network Robustness: A Viewpoint of Robust Growth of 
# Generalized Meshes. PLOS ONE. https://doi.org/10.1371/journal.pone.0161077
##########################################################################################
# z = a nx1 matrix (1 column), containing behavioral observations
# n.bootstraps = n of repetitions per removal proportion
# pr.thresh = threshold for reduced data summary variable/unmanipulated summary variable 
# show. prog = logical (default TRUE) to print progress for each
BOOTSTRAP.networksummarizer<-function(z,n.bootstraps=100,pr.thresh=0.25,
                                      which.var.ref="mean_path_length",show.prog=TRUE){
  
  
  OriginalFullNetwork<-Network.Quick(z)
  OriginalFullNetwork<-data.frame(1, OriginalFullNetwork);colnames(OriginalFullNetwork)[1]<-"HowMuch"
  
  remove.dis<-seq(0.01,0.99,.01)
  
  #Bootstrap ---subset and recalculate
  md<-1
  
  while(md<100){

    HowMuch<-1-remove.dis[md]
    
    if(show.prog==TRUE)
      print(paste("Remove ",remove.dis[md]," of data",sep=''))
    
    
    for(xar in 1:n.bootstraps){
      new.z<-sample(z,round(nrow(z)*HowMuch),replace=FALSE)
      SubNetwork<-Network.Quick(new.z)
      SubNetwork<-data.frame(HowMuch,SubNetwork)
      
      if(xar==1)
        SubFrame<-SubNetwork
      
      if(xar>1)
        SubFrame<-rbind(SubFrame,SubNetwork)
      
      
      if(xar==n.bootstraps){
        
        relative<-t(matrix((unlist(apply(SubFrame, 1,function(x) x/OriginalFullNetwork))),ncol=100))
        colnames(relative)<-colnames(SubFrame)
        
        Line.Item<-apply(SubFrame,2,function(x) median(x))
        SD.Lines<-((apply(SubFrame,2,sd)));names(SD.Lines)<-paste(names(SD.Lines),".SD",sep='')
        Line.Item2<-c(Line.Item,SD.Lines[c(2:16)])
      }
    }
    
    #Line.Item2
    
    RatioForRemoved<-data.frame(Line.Item/OriginalFullNetwork)
    
    if(md==1)
      RemovalSummary<-Line.Item2
    
    if(md>1)
      RemovalSummary<-rbind(RemovalSummary,Line.Item2)
    
    md<-md+1
    #IF ratio for defined variable falls below threshold, this will stop the loop
    if(RatioForRemoved[,which.var.ref]<pr.thresh){
        md<-100
    }
    
  }
    
  return(RemovalSummary)
  
}





###################
RAL.hack<-function (wp.P, piclist,woodpeckerscalar=0.02,optimizeangle=FALSE,rightontip=TRUE,ladder=TRUE) 
{
  #plot(tree.p,type="fan")
  tree_view<-as.phylo((hclust(wp.P)))
  df <- tree_view$tip.label
  idx <- match(names(piclist), df)
  
  if(anyNA(idx))
    print("Names of PNGs and names used in distance calculation do not match")
  
  if(ladder==TRUE){
    plot(ladderize(as.phylo((hclust(wp.P))), cex = 0.05, label.offset = 1),type = "fan")
  } else {
    plot(as.phylo((hclust(wp.P)), cex = 0.05, label.offset = 1),type = "fan")
    
  }
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tip <- 1:lastPP$Ntip
  XX <- lastPP$xx[tip]
  YY <- lastPP$yy[tip]
  
  xrange<-max(XX)-min(XX)
  yrange<-max(YY)-min(YY)
  
  maxrange<-max(c(xrange,yrange))
  
  for(h in 1:length(piclist)){
    img<-piclist[[h]]
    myX<-XX[h]
    myY<-YY[h]
    
    
    
    if(rightontip==TRUE){
      
      xR<-(myX)+(woodpeckerscalar*maxrange/2)
      xL<-(myX)-(woodpeckerscalar*maxrange/2)
      yT<-(myY)+(woodpeckerscalar*maxrange/2)
      yB<-(myY)-(woodpeckerscalar*maxrange/2)
    }
    if(rightontip==FALSE){
      
      if(myY==0 && myX>0){ #middle right
        xL<-myX
        xR<-myX+woodpeckerscalar*maxrange
        yT<-woodpeckerscalar*maxrange/2
        yB<-yT*-1
      } 
      
      
      
      if(myY>0 && myX>0){ #top right
        xL<-myX
        xR<-myX+woodpeckerscalar*maxrange
        yB<-myY
        yT<-yB+woodpeckerscalar*maxrange
      } 
      
      if(myY>0 && myX<0){#top left
        xL<-myX-woodpeckerscalar*maxrange
        xR<-myX     
        yB<-myY
        yT<-yB+woodpeckerscalar*maxrange
      }
      if(myY<0 && myX>0){#bottom right
        xL<-myX
        xR<-myX+woodpeckerscalar*maxrange
        yT<-myY
        yB<-yT-woodpeckerscalar*maxrange
      }    
      if(myY<0 && myX<0){#bottom left
        xL<-myX-woodpeckerscalar*maxrange
        xR<-myX
        yT<-myY
        yB<-yT-woodpeckerscalar*maxrange
      }
      
    }
    
    ############################ Calculate angle between right horizontal and bird/tip location, to provide angle param for plot
    xa<-c(1,0,myX)
    ya<-c(0,0,myY)
    d <- diff(complex(real = xa, imaginary = ya))
    
    if(optimizeangle==TRUE){
      angle1<-(diff(Arg(d)) %% (2*pi)) * 360/(2*pi)
    } else {
      angle1<-0
    }
    ###########
    
    rasterImage(img, xleft = xL, ybottom=yB, xright=xR,ytop=yT,angle=angle1)
    
    
  }
}  

# WoodpeckerSpace ---------------------------------------------------------
#Function for plotting Woodpeckers in space
WoodpeckerSpace<-function(SpaceCoords,piclist,
                          woodpeckerscalar.X=0.02,
                          woodpeckerscalar.Y=0.02,
                          woodpeckerscalar=NA,
                          textcolor="black",textsize=1,
                          label.x=NA,
                          label.y=NA,
                          ylabelplacement=2,xlabelplacement=2,
                          ynumplacement=2,xnumplacement=2,
                          axisnumsize=1,axislabsize=1,xlabelcolor="green",ylabelcolor="blue",bxtype="l"){
  
  SpaceCoords<-data.frame(SpaceCoords)
  
  if(!is.na(label.x)){
    plot(SpaceCoords,cex.lab=axislabsize,cex.axis=axisnumsize,col.lab = xlabelcolor,bty=bxtype,
         xlab = '',ylab='',
         col.axis=xlabelcolor,yaxt='n',xaxt='n')
    # plot(SpaceCoords,cex.lab=axislabsize,cex.axis=axisnumsize,col.lab = xlabelcolor,bty=bxtype,
    #      xlab = label.x,ylab=label.y,
    #      col.axis=xlabelcolor,yaxt='n')
    par(mgp=c(3,ynumplacement,0))
    axis(2, col.lab = ylabelcolor, col.axis = ylabelcolor,cex.axis=axisnumsize)
    title(ylab=label.y,col.lab=ylabelcolor,cex.lab=axislabsize,cex.axis=axisnumsize,line=ylabelplacement)
    
    par(mgp=c(3,xnumplacement,0))
    axis(1, col.lab = xlabelcolor, col.axis = xlabelcolor,cex.axis=axisnumsize)
    title(xlab=label.x,col.lab=xlabelcolor,cex.lab=axislabsize,cex.axis=axisnumsize,line=xlabelplacement)
    
    
    
    
  } else {
    
    plot(SpaceCoords,cex.lab=axislabsize,cex.axis=axisnumsize,col.lab = xlabelcolor,bty=bxtype,
         col.axis=xlabelcolor,yaxt='n',ylab="")
    axis(2, col.lab = ylabelcolor, col.axis = ylabelcolor,cex.axis=axisnumsize)
    title(ylab=colnames(SpaceCoords)[2],col.lab=ylabelcolor,cex.lab=axislabsize,cex.axis=axisnumsize)
  }
  
  
  if(textcolor=="OFF"){
    
  } else {
    text(SpaceCoords,labels=row.names(SpaceCoords),col=textcolor,cex=textsize)
  }
  
  maxrange<-max(c(range(SpaceCoords[,1])[1],range(SpaceCoords[,1])[2]))-  min(c(range(SpaceCoords[,1])[1],range(SpaceCoords[,1])[2]))
  maxYrange<-max(c(range(SpaceCoords[,2])[1],range(SpaceCoords[,2])[2]))-  min(c(range(SpaceCoords[,2])[1],range(SpaceCoords[,2])[2]))
  
  if(maxYrange==0)
    only1Y<-max(c(range(SpaceCoords[,2])[1],range(SpaceCoords[,2])[2]))
  
  #calculate max ranges
  
  if(all(names(piclist)==rownames(SpaceCoords))){
    
    for(h in 1:length(piclist)){
      img<-piclist[[h]]
      myX<-SpaceCoords[h,1]
      myY<-SpaceCoords[h,2]
      
      
      if(is.na(woodpeckerscalar)){
        xR<-(myX)+(woodpeckerscalar.X*maxrange/2)
        xL<-(myX)-(woodpeckerscalar.X*maxrange/2)
        yT<-(myY)+(woodpeckerscalar.Y*maxYrange/2)
        yB<-(myY)-(woodpeckerscalar.Y*maxYrange/2)
        
        if(maxYrange==0){
          ybottom.a=only1Y-(woodpeckerscalar.Y/2)
          ytop.a=only1Y+(woodpeckerscalar.Y/2)
          rasterImage(img, xleft = xL, ybottom=ybottom.a, xright=xR,ytop=ytop.a)
        } else {
          rasterImage(img, xleft = xL, ybottom=yB, xright=xR,ytop=yT)
        }
      } else {
        xR<-(myX)+(woodpeckerscalar*maxrange/2)
        xL<-(myX)-(woodpeckerscalar*maxrange/2)
        yT<-(myY)+(woodpeckerscalar*maxYrange/2)
        yB<-(myY)-(woodpeckerscalar*maxYrange/2)
        
        if(maxYrange==0){
          ybottom.a=only1Y-(woodpeckerscalar/2)
          ytop.a=only1Y+(woodpeckerscalar/2)
          rasterImage(img, xleft = xL, ybottom=ybottom.a, xright=xR,ytop=ytop.a)
        } else {
          rasterImage(img, xleft = xL, ybottom=yB, xright=xR,ytop=yT)
        }
      }
      
      
      
    }
    
    
  } else {
    print("Names of images do not match rownames")
  }
  reset_par()
}



#RESET PAR FUNCTION----
reset_par <- function(){
  op <- structure(list(xlog = FALSE, ylog = FALSE, adj = 0.5, ann = TRUE,
                       ask = FALSE, bg = "transparent", bty = "o", cex = 1, cex.axis = 1,
                       cex.lab = 1, cex.main = 1.2, cex.sub = 1, col = "black",
                       col.axis = "black", col.lab = "black", col.main = "black",
                       col.sub = "black", crt = 0, err = 0L, family = "", fg = "black",
                       fig = c(0, 1, 0, 1), fin = c(6.99999895833333, 6.99999895833333
                       ), font = 1L, font.axis = 1L, font.lab = 1L, font.main = 2L,
                       font.sub = 1L, lab = c(5L, 5L, 7L), las = 0L, lend = "round",
                       lheight = 1, ljoin = "round", lmitre = 10, lty = "solid",
                       lwd = 1, mai = c(1.02, 0.82, 0.82, 0.42), mar = c(5.1, 4.1,
                                                                         4.1, 2.1), mex = 1, mfcol = c(1L, 1L), mfg = c(1L, 1L, 1L,
                                                                                                                        1L), mfrow = c(1L, 1L), mgp = c(3, 1, 0), mkh = 0.001, new = FALSE,
                       oma = c(0, 0, 0, 0), omd = c(0, 1, 0, 1), omi = c(0, 0, 0,
                                                                         0), pch = 1L, pin = c(5.75999895833333, 5.15999895833333),
                       plt = c(0.117142874574832, 0.939999991071427, 0.145714307397962,
                               0.882857125425167), ps = 12L, pty = "m", smo = 1, srt = 0,
                       tck = NA_real_, tcl = -0.5, usr = c(0.568, 1.432, 0.568,
                                                           1.432), xaxp = c(0.6, 1.4, 4), xaxs = "r", xaxt = "s", xpd = FALSE,
                       yaxp = c(0.6, 1.4, 4), yaxs = "r", yaxt = "s", ylbias = 0.2), .Names = c("xlog",
                                                                                                "ylog", "adj", "ann", "ask", "bg", "bty", "cex", "cex.axis",
                                                                                                "cex.lab", "cex.main", "cex.sub", "col", "col.axis", "col.lab",
                                                                                                "col.main", "col.sub", "crt", "err", "family", "fg", "fig", "fin",
                                                                                                "font", "font.axis", "font.lab", "font.main", "font.sub", "lab",
                                                                                                "las", "lend", "lheight", "ljoin", "lmitre", "lty", "lwd", "mai",
                                                                                                "mar", "mex", "mfcol", "mfg", "mfrow", "mgp", "mkh", "new", "oma",
                                                                                                "omd", "omi", "pch", "pin", "plt", "ps", "pty", "smo", "srt",
                                                                                                "tck", "tcl", "usr", "xaxp", "xaxs", "xaxt", "xpd", "yaxp", "yaxs",
                                                                                                "yaxt", "ylbias"))
  par(op)
}
###############
OnetoHundred.1000<-c(0, 0.1001001001001, 0.2002002002002, 0.3003003003003, 0.4004004004004, 
                     0.500500500500501, 0.600600600600601, 0.700700700700701, 0.800800800800801, 
                     0.900900900900901, 1.001001001001, 1.1011011011011, 1.2012012012012, 
                     1.3013013013013, 1.4014014014014, 1.5015015015015, 1.6016016016016, 
                     1.7017017017017, 1.8018018018018, 1.9019019019019, 2.002002002002, 
                     2.1021021021021, 2.2022022022022, 2.3023023023023, 2.4024024024024, 
                     2.5025025025025, 2.6026026026026, 2.7027027027027, 2.8028028028028, 
                     2.9029029029029, 3.003003003003, 3.1031031031031, 3.2032032032032, 
                     3.3033033033033, 3.4034034034034, 3.5035035035035, 3.6036036036036, 
                     3.7037037037037, 3.8038038038038, 3.9039039039039, 4.004004004004, 
                     4.1041041041041, 4.2042042042042, 4.3043043043043, 4.4044044044044, 
                     4.5045045045045, 4.6046046046046, 4.7047047047047, 4.8048048048048, 
                     4.9049049049049, 5.00500500500501, 5.10510510510511, 5.20520520520521, 
                     5.30530530530531, 5.40540540540541, 5.50550550550551, 5.60560560560561, 
                     5.70570570570571, 5.80580580580581, 5.90590590590591, 6.00600600600601, 
                     6.10610610610611, 6.20620620620621, 6.30630630630631, 6.40640640640641, 
                     6.50650650650651, 6.60660660660661, 6.70670670670671, 6.80680680680681, 
                     6.90690690690691, 7.00700700700701, 7.10710710710711, 7.20720720720721, 
                     7.30730730730731, 7.40740740740741, 7.50750750750751, 7.60760760760761, 
                     7.70770770770771, 7.80780780780781, 7.90790790790791, 8.00800800800801, 
                     8.10810810810811, 8.20820820820821, 8.30830830830831, 8.40840840840841, 
                     8.50850850850851, 8.60860860860861, 8.70870870870871, 8.80880880880881, 
                     8.90890890890891, 9.00900900900901, 9.10910910910911, 9.20920920920921, 
                     9.30930930930931, 9.40940940940941, 9.50950950950951, 9.60960960960961, 
                     9.70970970970971, 9.80980980980981, 9.90990990990991, 10.01001001001, 
                     10.1101101101101, 10.2102102102102, 10.3103103103103, 10.4104104104104, 
                     10.5105105105105, 10.6106106106106, 10.7107107107107, 10.8108108108108, 
                     10.9109109109109, 11.011011011011, 11.1111111111111, 11.2112112112112, 
                     11.3113113113113, 11.4114114114114, 11.5115115115115, 11.6116116116116, 
                     11.7117117117117, 11.8118118118118, 11.9119119119119, 12.012012012012, 
                     12.1121121121121, 12.2122122122122, 12.3123123123123, 12.4124124124124, 
                     12.5125125125125, 12.6126126126126, 12.7127127127127, 12.8128128128128, 
                     12.9129129129129, 13.013013013013, 13.1131131131131, 13.2132132132132, 
                     13.3133133133133, 13.4134134134134, 13.5135135135135, 13.6136136136136, 
                     13.7137137137137, 13.8138138138138, 13.9139139139139, 14.014014014014, 
                     14.1141141141141, 14.2142142142142, 14.3143143143143, 14.4144144144144, 
                     14.5145145145145, 14.6146146146146, 14.7147147147147, 14.8148148148148, 
                     14.9149149149149, 15.015015015015, 15.1151151151151, 15.2152152152152, 
                     15.3153153153153, 15.4154154154154, 15.5155155155155, 15.6156156156156, 
                     15.7157157157157, 15.8158158158158, 15.9159159159159, 16.016016016016, 
                     16.1161161161161, 16.2162162162162, 16.3163163163163, 16.4164164164164, 
                     16.5165165165165, 16.6166166166166, 16.7167167167167, 16.8168168168168, 
                     16.9169169169169, 17.017017017017, 17.1171171171171, 17.2172172172172, 
                     17.3173173173173, 17.4174174174174, 17.5175175175175, 17.6176176176176, 
                     17.7177177177177, 17.8178178178178, 17.9179179179179, 18.018018018018, 
                     18.1181181181181, 18.2182182182182, 18.3183183183183, 18.4184184184184, 
                     18.5185185185185, 18.6186186186186, 18.7187187187187, 18.8188188188188, 
                     18.9189189189189, 19.019019019019, 19.1191191191191, 19.2192192192192, 
                     19.3193193193193, 19.4194194194194, 19.5195195195195, 19.6196196196196, 
                     19.7197197197197, 19.8198198198198, 19.9199199199199, 20.02002002002, 
                     20.1201201201201, 20.2202202202202, 20.3203203203203, 20.4204204204204, 
                     20.5205205205205, 20.6206206206206, 20.7207207207207, 20.8208208208208, 
                     20.9209209209209, 21.021021021021, 21.1211211211211, 21.2212212212212, 
                     21.3213213213213, 21.4214214214214, 21.5215215215215, 21.6216216216216, 
                     21.7217217217217, 21.8218218218218, 21.9219219219219, 22.022022022022, 
                     22.1221221221221, 22.2222222222222, 22.3223223223223, 22.4224224224224, 
                     22.5225225225225, 22.6226226226226, 22.7227227227227, 22.8228228228228, 
                     22.9229229229229, 23.023023023023, 23.1231231231231, 23.2232232232232, 
                     23.3233233233233, 23.4234234234234, 23.5235235235235, 23.6236236236236, 
                     23.7237237237237, 23.8238238238238, 23.9239239239239, 24.024024024024, 
                     24.1241241241241, 24.2242242242242, 24.3243243243243, 24.4244244244244, 
                     24.5245245245245, 24.6246246246246, 24.7247247247247, 24.8248248248248, 
                     24.9249249249249, 25.025025025025, 25.1251251251251, 25.2252252252252, 
                     25.3253253253253, 25.4254254254254, 25.5255255255255, 25.6256256256256, 
                     25.7257257257257, 25.8258258258258, 25.9259259259259, 26.026026026026, 
                     26.1261261261261, 26.2262262262262, 26.3263263263263, 26.4264264264264, 
                     26.5265265265265, 26.6266266266266, 26.7267267267267, 26.8268268268268, 
                     26.9269269269269, 27.027027027027, 27.1271271271271, 27.2272272272272, 
                     27.3273273273273, 27.4274274274274, 27.5275275275275, 27.6276276276276, 
                     27.7277277277277, 27.8278278278278, 27.9279279279279, 28.028028028028, 
                     28.1281281281281, 28.2282282282282, 28.3283283283283, 28.4284284284284, 
                     28.5285285285285, 28.6286286286286, 28.7287287287287, 28.8288288288288, 
                     28.9289289289289, 29.029029029029, 29.1291291291291, 29.2292292292292, 
                     29.3293293293293, 29.4294294294294, 29.5295295295295, 29.6296296296296, 
                     29.7297297297297, 29.8298298298298, 29.9299299299299, 30.03003003003, 
                     30.1301301301301, 30.2302302302302, 30.3303303303303, 30.4304304304304, 
                     30.5305305305305, 30.6306306306306, 30.7307307307307, 30.8308308308308, 
                     30.9309309309309, 31.031031031031, 31.1311311311311, 31.2312312312312, 
                     31.3313313313313, 31.4314314314314, 31.5315315315315, 31.6316316316316, 
                     31.7317317317317, 31.8318318318318, 31.9319319319319, 32.032032032032, 
                     32.1321321321321, 32.2322322322322, 32.3323323323323, 32.4324324324324, 
                     32.5325325325325, 32.6326326326326, 32.7327327327327, 32.8328328328328, 
                     32.9329329329329, 33.033033033033, 33.1331331331331, 33.2332332332332, 
                     33.3333333333333, 33.4334334334334, 33.5335335335335, 33.6336336336336, 
                     33.7337337337337, 33.8338338338338, 33.9339339339339, 34.034034034034, 
                     34.1341341341341, 34.2342342342342, 34.3343343343343, 34.4344344344344, 
                     34.5345345345345, 34.6346346346346, 34.7347347347347, 34.8348348348348, 
                     34.9349349349349, 35.035035035035, 35.1351351351351, 35.2352352352352, 
                     35.3353353353353, 35.4354354354354, 35.5355355355355, 35.6356356356356, 
                     35.7357357357357, 35.8358358358358, 35.9359359359359, 36.036036036036, 
                     36.1361361361361, 36.2362362362362, 36.3363363363363, 36.4364364364364, 
                     36.5365365365365, 36.6366366366366, 36.7367367367367, 36.8368368368368, 
                     36.9369369369369, 37.037037037037, 37.1371371371371, 37.2372372372372, 
                     37.3373373373373, 37.4374374374374, 37.5375375375375, 37.6376376376376, 
                     37.7377377377377, 37.8378378378378, 37.9379379379379, 38.038038038038, 
                     38.1381381381381, 38.2382382382382, 38.3383383383383, 38.4384384384384, 
                     38.5385385385385, 38.6386386386386, 38.7387387387387, 38.8388388388388, 
                     38.9389389389389, 39.039039039039, 39.1391391391391, 39.2392392392392, 
                     39.3393393393393, 39.4394394394394, 39.5395395395395, 39.6396396396396, 
                     39.7397397397397, 39.8398398398398, 39.9399399399399, 40.04004004004, 
                     40.1401401401401, 40.2402402402402, 40.3403403403403, 40.4404404404404, 
                     40.5405405405405, 40.6406406406406, 40.7407407407407, 40.8408408408408, 
                     40.9409409409409, 41.041041041041, 41.1411411411411, 41.2412412412412, 
                     41.3413413413413, 41.4414414414414, 41.5415415415415, 41.6416416416416, 
                     41.7417417417417, 41.8418418418418, 41.9419419419419, 42.042042042042, 
                     42.1421421421421, 42.2422422422422, 42.3423423423423, 42.4424424424424, 
                     42.5425425425425, 42.6426426426426, 42.7427427427427, 42.8428428428428, 
                     42.9429429429429, 43.043043043043, 43.1431431431431, 43.2432432432432, 
                     43.3433433433433, 43.4434434434434, 43.5435435435435, 43.6436436436436, 
                     43.7437437437437, 43.8438438438438, 43.9439439439439, 44.044044044044, 
                     44.1441441441441, 44.2442442442442, 44.3443443443443, 44.4444444444444, 
                     44.5445445445446, 44.6446446446447, 44.7447447447448, 44.8448448448449, 
                     44.944944944945, 45.0450450450451, 45.1451451451452, 45.2452452452453, 
                     45.3453453453454, 45.4454454454455, 45.5455455455456, 45.6456456456457, 
                     45.7457457457458, 45.8458458458459, 45.945945945946, 46.0460460460461, 
                     46.1461461461462, 46.2462462462463, 46.3463463463464, 46.4464464464465, 
                     46.5465465465466, 46.6466466466467, 46.7467467467468, 46.8468468468469, 
                     46.946946946947, 47.0470470470471, 47.1471471471472, 47.2472472472473, 
                     47.3473473473474, 47.4474474474475, 47.5475475475476, 47.6476476476477, 
                     47.7477477477478, 47.8478478478479, 47.947947947948, 48.0480480480481, 
                     48.1481481481482, 48.2482482482483, 48.3483483483484, 48.4484484484485, 
                     48.5485485485486, 48.6486486486487, 48.7487487487488, 48.8488488488489, 
                     48.948948948949, 49.0490490490491, 49.1491491491492, 49.2492492492493, 
                     49.3493493493494, 49.4494494494495, 49.5495495495496, 49.6496496496497, 
                     49.7497497497498, 49.8498498498499, 49.94994994995, 50.0500500500501, 
                     50.1501501501502, 50.2502502502503, 50.3503503503504, 50.4504504504505, 
                     50.5505505505506, 50.6506506506507, 50.7507507507508, 50.8508508508509, 
                     50.950950950951, 51.0510510510511, 51.1511511511512, 51.2512512512513, 
                     51.3513513513514, 51.4514514514515, 51.5515515515516, 51.6516516516517, 
                     51.7517517517518, 51.8518518518519, 51.951951951952, 52.0520520520521, 
                     52.1521521521522, 52.2522522522523, 52.3523523523524, 52.4524524524525, 
                     52.5525525525526, 52.6526526526527, 52.7527527527528, 52.8528528528529, 
                     52.952952952953, 53.0530530530531, 53.1531531531532, 53.2532532532533, 
                     53.3533533533534, 53.4534534534535, 53.5535535535536, 53.6536536536537, 
                     53.7537537537538, 53.8538538538539, 53.953953953954, 54.0540540540541, 
                     54.1541541541542, 54.2542542542543, 54.3543543543544, 54.4544544544545, 
                     54.5545545545546, 54.6546546546547, 54.7547547547548, 54.8548548548549, 
                     54.954954954955, 55.0550550550551, 55.1551551551552, 55.2552552552553, 
                     55.3553553553554, 55.4554554554555, 55.5555555555556, 55.6556556556557, 
                     55.7557557557558, 55.8558558558559, 55.955955955956, 56.0560560560561, 
                     56.1561561561562, 56.2562562562563, 56.3563563563564, 56.4564564564565, 
                     56.5565565565566, 56.6566566566567, 56.7567567567568, 56.8568568568569, 
                     56.956956956957, 57.0570570570571, 57.1571571571572, 57.2572572572573, 
                     57.3573573573574, 57.4574574574575, 57.5575575575576, 57.6576576576577, 
                     57.7577577577578, 57.8578578578579, 57.957957957958, 58.0580580580581, 
                     58.1581581581582, 58.2582582582583, 58.3583583583584, 58.4584584584585, 
                     58.5585585585586, 58.6586586586587, 58.7587587587588, 58.8588588588589, 
                     58.958958958959, 59.0590590590591, 59.1591591591592, 59.2592592592593, 
                     59.3593593593594, 59.4594594594595, 59.5595595595596, 59.6596596596597, 
                     59.7597597597598, 59.8598598598599, 59.95995995996, 60.0600600600601, 
                     60.1601601601602, 60.2602602602603, 60.3603603603604, 60.4604604604605, 
                     60.5605605605606, 60.6606606606607, 60.7607607607608, 60.8608608608609, 
                     60.960960960961, 61.0610610610611, 61.1611611611612, 61.2612612612613, 
                     61.3613613613614, 61.4614614614615, 61.5615615615616, 61.6616616616617, 
                     61.7617617617618, 61.8618618618619, 61.961961961962, 62.0620620620621, 
                     62.1621621621622, 62.2622622622623, 62.3623623623624, 62.4624624624625, 
                     62.5625625625626, 62.6626626626627, 62.7627627627628, 62.8628628628629, 
                     62.962962962963, 63.0630630630631, 63.1631631631632, 63.2632632632633, 
                     63.3633633633634, 63.4634634634635, 63.5635635635636, 63.6636636636637, 
                     63.7637637637638, 63.8638638638639, 63.963963963964, 64.0640640640641, 
                     64.1641641641642, 64.2642642642643, 64.3643643643644, 64.4644644644645, 
                     64.5645645645646, 64.6646646646647, 64.7647647647648, 64.8648648648649, 
                     64.964964964965, 65.0650650650651, 65.1651651651652, 65.2652652652653, 
                     65.3653653653654, 65.4654654654655, 65.5655655655656, 65.6656656656657, 
                     65.7657657657658, 65.8658658658659, 65.965965965966, 66.0660660660661, 
                     66.1661661661662, 66.2662662662663, 66.3663663663664, 66.4664664664665, 
                     66.5665665665666, 66.6666666666667, 66.7667667667668, 66.8668668668669, 
                     66.966966966967, 67.0670670670671, 67.1671671671672, 67.2672672672673, 
                     67.3673673673674, 67.4674674674675, 67.5675675675676, 67.6676676676677, 
                     67.7677677677678, 67.8678678678679, 67.967967967968, 68.0680680680681, 
                     68.1681681681682, 68.2682682682683, 68.3683683683684, 68.4684684684685, 
                     68.5685685685686, 68.6686686686687, 68.7687687687688, 68.8688688688689, 
                     68.968968968969, 69.0690690690691, 69.1691691691692, 69.2692692692693, 
                     69.3693693693694, 69.4694694694695, 69.5695695695696, 69.6696696696697, 
                     69.7697697697698, 69.8698698698699, 69.96996996997, 70.0700700700701, 
                     70.1701701701702, 70.2702702702703, 70.3703703703704, 70.4704704704705, 
                     70.5705705705706, 70.6706706706707, 70.7707707707708, 70.8708708708709, 
                     70.970970970971, 71.0710710710711, 71.1711711711712, 71.2712712712713, 
                     71.3713713713714, 71.4714714714715, 71.5715715715716, 71.6716716716717, 
                     71.7717717717718, 71.8718718718719, 71.971971971972, 72.0720720720721, 
                     72.1721721721722, 72.2722722722723, 72.3723723723724, 72.4724724724725, 
                     72.5725725725726, 72.6726726726727, 72.7727727727728, 72.8728728728729, 
                     72.972972972973, 73.0730730730731, 73.1731731731732, 73.2732732732733, 
                     73.3733733733734, 73.4734734734735, 73.5735735735736, 73.6736736736737, 
                     73.7737737737738, 73.8738738738739, 73.973973973974, 74.0740740740741, 
                     74.1741741741742, 74.2742742742743, 74.3743743743744, 74.4744744744745, 
                     74.5745745745746, 74.6746746746747, 74.7747747747748, 74.8748748748749, 
                     74.974974974975, 75.0750750750751, 75.1751751751752, 75.2752752752753, 
                     75.3753753753754, 75.4754754754755, 75.5755755755756, 75.6756756756757, 
                     75.7757757757758, 75.8758758758759, 75.975975975976, 76.0760760760761, 
                     76.1761761761762, 76.2762762762763, 76.3763763763764, 76.4764764764765, 
                     76.5765765765766, 76.6766766766767, 76.7767767767768, 76.8768768768769, 
                     76.976976976977, 77.0770770770771, 77.1771771771772, 77.2772772772773, 
                     77.3773773773774, 77.4774774774775, 77.5775775775776, 77.6776776776777, 
                     77.7777777777778, 77.8778778778779, 77.977977977978, 78.0780780780781, 
                     78.1781781781782, 78.2782782782783, 78.3783783783784, 78.4784784784785, 
                     78.5785785785786, 78.6786786786787, 78.7787787787788, 78.8788788788789, 
                     78.978978978979, 79.0790790790791, 79.1791791791792, 79.2792792792793, 
                     79.3793793793794, 79.4794794794795, 79.5795795795796, 79.6796796796797, 
                     79.7797797797798, 79.8798798798799, 79.97997997998, 80.0800800800801, 
                     80.1801801801802, 80.2802802802803, 80.3803803803804, 80.4804804804805, 
                     80.5805805805806, 80.6806806806807, 80.7807807807808, 80.8808808808809, 
                     80.980980980981, 81.0810810810811, 81.1811811811812, 81.2812812812813, 
                     81.3813813813814, 81.4814814814815, 81.5815815815816, 81.6816816816817, 
                     81.7817817817818, 81.8818818818819, 81.981981981982, 82.0820820820821, 
                     82.1821821821822, 82.2822822822823, 82.3823823823824, 82.4824824824825, 
                     82.5825825825826, 82.6826826826827, 82.7827827827828, 82.8828828828829, 
                     82.982982982983, 83.0830830830831, 83.1831831831832, 83.2832832832833, 
                     83.3833833833834, 83.4834834834835, 83.5835835835836, 83.6836836836837, 
                     83.7837837837838, 83.8838838838839, 83.983983983984, 84.0840840840841, 
                     84.1841841841842, 84.2842842842843, 84.3843843843844, 84.4844844844845, 
                     84.5845845845846, 84.6846846846847, 84.7847847847848, 84.8848848848849, 
                     84.984984984985, 85.0850850850851, 85.1851851851852, 85.2852852852853, 
                     85.3853853853854, 85.4854854854855, 85.5855855855856, 85.6856856856857, 
                     85.7857857857858, 85.8858858858859, 85.985985985986, 86.0860860860861, 
                     86.1861861861862, 86.2862862862863, 86.3863863863864, 86.4864864864865, 
                     86.5865865865866, 86.6866866866867, 86.7867867867868, 86.8868868868869, 
                     86.986986986987, 87.0870870870871, 87.1871871871872, 87.2872872872873, 
                     87.3873873873874, 87.4874874874875, 87.5875875875876, 87.6876876876877, 
                     87.7877877877878, 87.8878878878879, 87.987987987988, 88.0880880880881, 
                     88.1881881881882, 88.2882882882883, 88.3883883883884, 88.4884884884885, 
                     88.5885885885886, 88.6886886886887, 88.7887887887888, 88.8888888888889, 
                     88.988988988989, 89.0890890890891, 89.1891891891892, 89.2892892892893, 
                     89.3893893893894, 89.4894894894895, 89.5895895895896, 89.6896896896897, 
                     89.7897897897898, 89.8898898898899, 89.98998998999, 90.0900900900901, 
                     90.1901901901902, 90.2902902902903, 90.3903903903904, 90.4904904904905, 
                     90.5905905905906, 90.6906906906907, 90.7907907907908, 90.8908908908909, 
                     90.990990990991, 91.0910910910911, 91.1911911911912, 91.2912912912913, 
                     91.3913913913914, 91.4914914914915, 91.5915915915916, 91.6916916916917, 
                     91.7917917917918, 91.8918918918919, 91.991991991992, 92.0920920920921, 
                     92.1921921921922, 92.2922922922923, 92.3923923923924, 92.4924924924925, 
                     92.5925925925926, 92.6926926926927, 92.7927927927928, 92.8928928928929, 
                     92.992992992993, 93.0930930930931, 93.1931931931932, 93.2932932932933, 
                     93.3933933933934, 93.4934934934935, 93.5935935935936, 93.6936936936937, 
                     93.7937937937938, 93.8938938938939, 93.993993993994, 94.0940940940941, 
                     94.1941941941942, 94.2942942942943, 94.3943943943944, 94.4944944944945, 
                     94.5945945945946, 94.6946946946947, 94.7947947947948, 94.8948948948949, 
                     94.994994994995, 95.0950950950951, 95.1951951951952, 95.2952952952953, 
                     95.3953953953954, 95.4954954954955, 95.5955955955956, 95.6956956956957, 
                     95.7957957957958, 95.8958958958959, 95.995995995996, 96.0960960960961, 
                     96.1961961961962, 96.2962962962963, 96.3963963963964, 96.4964964964965, 
                     96.5965965965966, 96.6966966966967, 96.7967967967968, 96.8968968968969, 
                     96.996996996997, 97.0970970970971, 97.1971971971972, 97.2972972972973, 
                     97.3973973973974, 97.4974974974975, 97.5975975975976, 97.6976976976977, 
                     97.7977977977978, 97.8978978978979, 97.997997997998, 98.0980980980981, 
                     98.1981981981982, 98.2982982982983, 98.3983983983984, 98.4984984984985, 
                     98.5985985985986, 98.6986986986987, 98.7987987987988, 98.8988988988989, 
                     98.998998998999, 99.0990990990991, 99.1991991991992, 99.2992992992993, 
                     99.3993993993994, 99.4994994994995, 99.5995995995996, 99.6996996996997, 
                     99.7997997997998, 99.8998998998999, 100)
