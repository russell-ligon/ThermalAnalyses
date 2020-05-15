##############
### ProcessBOP ###
##############



library("reshape2")
#library("magic")
library("data.table")
#ProcessBOP is a collection of functions written for dealing with output from the CowLog behavioral analysis program modified for Birds of Paradise. 
#ProcessBOP was written by Russell Ligon (2016).

#ProcessBOP
#This function extracts values from a BOPLog output file and requires the following variables for input:
#filename - full pathway to desired directory containing BOPLog table (in default, csv format, see example table)
# e.g.  "C:/Users/Rusty/Amazon Drive/BOP/Video_analysis/BOPlog output_oneatatime/455287_2000-01-01T000000_000.csv"
#USE "filename<-file.names[16]" when troubleshooting ProcessBOP
ProcessBOP<-function(filename,location,species.name,specialclass,hb){

#function for pulling right side of text string, for naming purposes    
  right = function (string, char){
    substr(string,nchar(string)-(char-1),nchar(string))
  }
  hangingbird<-hb
######Function for eliminating event behaviors during "OFF" and fixing duration behaviors
  #DURATION BEHAVIORS
  DurationBehaviors<-c("BP1","BP3","SS1","SS2","O3","O4","O5","OFF","O6","O7","O8","OPMH.D","OPMB.D","OPMF.D","OPMW.D","OPMT.D","RepWing","RepTail","RepHead","RepTorso",
                       "RM1","RM2","RM3","RM4","MO2","PU1","SP")
  DurationBehaviors.nosp<-c("BP1","BP3","SS1","SS2","O3","O4","O5","OFF","O6","O7","O8","OPMH.D","OPMB.D","OPMF.D","OPMW.D","OPMT.D","RepWing","RepTail","RepHead","RepTorso",
                            "RM1","RM2","RM3","RM4","MO2","PU1")
  DurationBehaviors.nosp.noOFF<-c("BP1","BP3","SS1","SS2","O3","O4","O5","O6","O7","O8","OPMH.D","OPMB.D","OPMF.D","OPMW.D","OPMT.D","RepWing","RepTail","RepHead","RepTorso",
                            "RM1","RM2","RM3","RM4","MO2","PU1")
  
  EventBehaviors<-c("BP2","O1","O2","OPMH","OPMB","OPMF","OPMW","OPMT","OPAH","OPAC1","OPAC2","OPAC3","OPAW","OPATW","OPAMW","OPAMTT",
                    "OPABH","OPABB","OPABF","OPABW","OPABT","PC1","PC2","PC3","PC4","MO1","SP1","SP2","SP3","SP4")
  
fixOFFs<-function(tabletofix){  
  offstarts<-nrow(subset(tabletofix,code=="OFF(start)"))
  offstops<-nrow(subset(tabletofix,code=="OFF(end)"))	
  tabletofix<-tabletofix[order(tabletofix$time),]
  if(offstarts>0){#If there are "OFFS", this if statement eliminates any behaviors measured during these OFFs and adds "STOPS"
    #At the beginning of the OFF period and adds "STARTS" at the end of the OFF period. Also deletes behaviors that happen while OFF.
   offstartings<-subset(tabletofix,code=="OFF(start)")[,1]
   offstoppins<-subset(tabletofix,code=="OFF(end)")[,1]
    
    for (vvv in 1:offstarts){#offstarts
      tabletofix<-tabletofix[order(tabletofix$time),]
      #stopBegin<-subset(tabletofix,code=="OFF(start)")[vvv,1]
      #stopEnd<-subset(tabletofix,code=="OFF(end)")[vvv,1]
      stopBegin<-offstartings[vvv]
      stopEnd<-offstoppins[vvv]
      
      events<- tabletofix[tabletofix$code %in% EventBehaviors,]
      gap.wise.misses<-events[which(events$time < stopEnd & events$time > stopBegin),]
      
      if(nrow(gap.wise.misses)>0){
        
        tabletofix<-tabletofix[-(which((tabletofix$time%in% gap.wise.misses$time) & (tabletofix$code%in% gap.wise.misses$code))),] #Removes rows of eventbehaviors which happen while birds is OFF screen
        #[-(which(tabletofix$time == oldend$time & tabletofix$code==oldend$code)),]
      }
      #testlist<-"O4"
      for (db in DurationBehaviors.nosp){  ##testlist #DurationBehaviors.nosp
        behav.obs <- tabletofix[grep(db, tabletofix$code), ]
        bouts<-nrow(behav.obs)/2
        if(bouts>0.6){
          starts<-behav.obs[grep("start", behav.obs$code), ]
          ends<-behav.obs[grep("end", behav.obs$code), ]
          
          oldend<-starts[1,][-1,]
          newend<-starts[1,][-1,]
          oldstart<-starts[1,][-1,]
          newstart<-starts[1,][-1,]
          
          if(nrow(starts)==nrow(ends)){ #only do this if starts/ends are equal for a given behavior
            for(nnn in 1:nrow(starts)){
              if (starts[nnn,1]==stopBegin & ends[nnn,1]<stopEnd & ends[nnn,2]!="OFF(end)"){ #S3_ If start occurs as OFF starts, and end occurs while OFF is active, delete
                starttimetoeliminate<-starts[nnn,1]
                stoptimetoeliminate<-ends[nnn,1]
                startremove<-behav.obs[which(behav.obs$time == starttimetoeliminate),] 
                endremove<-behav.obs[which(behav.obs$time == stoptimetoeliminate),] 
                
                gap.wise.misses2<-rbind(startremove,endremove)
                
                tabletofix<-tabletofix[-((match(gap.wise.misses2$time,tabletofix$time))),]
              }
              
              if(starts[nnn,1]<stopBegin & ends[nnn,1]==stopEnd & ends[nnn,2]!="OFF(end)" ){#S2_This part handles all duration behaviors that start before each (nnn) OFF starts and end as each OFF ends.
                oldend<-ends[nnn,]
                
                
                newend<-ends[1,]
                newend$time<-stopBegin
                
                
                tabletofix<-rbind(tabletofix,newend) #adds new end at start of OFF period
                tabletofix<-tabletofix[-(which(tabletofix$time == oldend$time & tabletofix$code==oldend$code)),] #gets rid of old end behavior
              }  
              
              if(starts[nnn,1]<stopBegin & ends[nnn,1]<stopEnd & ends[nnn,1]>stopBegin & ends[nnn,2]!="OFF(end)"){#S1_This part handles all duration behaviors that start before each (nnn) OFF starts and end while OFF is on.
                oldend<-ends[nnn,]
                
                
                newend<-ends[1,]
                newend$time<-stopBegin
                
                
                tabletofix<-rbind(tabletofix,newend) #adds new end at start of OFF period
                tabletofix<-tabletofix[-(which(tabletofix$time == oldend$time & tabletofix$code==oldend$code)),] #gets rid of old end behavior
              }  
              
              if(starts[nnn,1]>stopBegin & starts[nnn,1]<stopEnd & ends[nnn,1]==stopEnd & ends[nnn,2]!="OFF(end)" ){#E4_This part handles all duration behaviors that start while OFF and ends with ON
                
                oldstart<-starts[nnn,]
                
                newstart<-starts[1,]
                newstart$time<-stopEnd
               
                tabletofix<-rbind(tabletofix,newstart) #
                tabletofix<-tabletofix[-(which(tabletofix$time == oldstart$time & tabletofix$code==oldstart$code)),] #gets rid of old start behavior
                
              }  
              
              
              if(starts[nnn,1]==stopBegin & ends[nnn,1]==stopEnd & ends[nnn,2]!="OFF(end)" ){#E3_This part handles all duration behaviors that start and end exactly as the OFF starts and ends
                oldend<-ends[nnn,]
                oldstart<-starts[nnn,]
                

                tabletofix<-tabletofix[-(which(tabletofix$time == oldend$time & tabletofix$code==oldend$code)),] #gets rid of old end behavior
                tabletofix<-tabletofix[-(which(tabletofix$time == oldstart$time & tabletofix$code==oldstart$code)),] #gets rid of old start behavior
                
              }  
              
              if(starts[nnn,1]==stopBegin & ends[nnn,1]>stopEnd & ends[nnn,2]!="OFF(end)" ){#E2_This part handles all duration behaviors that start as OFF starts and end after each OFF ends.
                oldstart<-starts[nnn,]
                
                
                newstart<-starts[1,]
                newstart$time<-stopEnd
                
                
                tabletofix<-rbind(tabletofix,newstart) #adds new end at start of OFF period
                tabletofix<-tabletofix[-(which(tabletofix$time == oldstart$time & tabletofix$code==oldstart$code)),] #gets rid of old end behavior
              } 
              
              if(starts[nnn,1]>stopBegin & starts[nnn,1]<stopEnd & ends[nnn,1]>stopEnd & ends[nnn,2]!="OFF(end)"){#E1_This part handles all duration behaviors that start while OFF is on, then after OFF is off
                oldstart<-starts[nnn,]
                
                
                newstart<-starts[1,]
                newstart$time<-stopEnd
                
                
                tabletofix<-rbind(tabletofix,newstart) #adds new end at start of OFF period
                tabletofix<-tabletofix[-(which(tabletofix$time == oldstart$time & tabletofix$code==oldstart$code)),] #gets rid of old end behavior
              }  
              
              
              if (starts[nnn,1]<stopBegin & ends[nnn,1]>stopEnd & ends[nnn,2]!="OFF(end)"){ #SPAN_if start happens before bird goes OFF AND stop happens after bird comes back (OFF(end)), cut out middle and add new stop and start
                starttimetoadd<-starts[nnn,]
                starttimetoadd$time<-stopEnd
                
                stoptimetoadd<-ends[nnn,]
                stoptimetoadd$time<-stopBegin
                 
                
                addstartstop<-rbind(starttimetoadd,stoptimetoadd)
                
                tabletofix<-rbind(tabletofix,addstartstop)
              }
              
              
              if (starts[nnn,1]>stopBegin & ends[nnn,1]<stopEnd & ends[nnn,2]!="OFF(end)"){ #MID_if start AND stop of a duration behavior occur while bird is OFF, delete
                starttimetoeliminate<-starts[nnn,1]
                stoptimetoeliminate<-ends[nnn,1]
                startremove<-behav.obs[which(behav.obs$time == starttimetoeliminate),] 
                endremove<-behav.obs[which(behav.obs$time == stoptimetoeliminate),] 
                
                gap.wise.misses2<-rbind(startremove,endremove)
                
                tabletofix<-tabletofix[-((match(gap.wise.misses2$time,tabletofix$time))),]
              }
              
              
            }# for every row in 'starts', which has all starts for a behavior, it calculates the difference between end and start and then adds this to 'duration'
          }    
        }
      }
      
    }
    
  } 
  
  #REMOVES behaviors that have a start and stop at the same time
  
  
  for (db2 in DurationBehaviors.nosp){  
    behav.observations <- tabletofix[grep(db2, tabletofix$code), ]
    startings<-behav.observations[grep("start",behav.observations$code),]
    endings<-behav.observations[grep("end",behav.observations$code),]
    
    startstoremove<-startings[(startings$time %in% endings$time),]
    endstoremove<-endings[(endings$time %in% startings$time),]
    kill<-rbind(startstoremove,endstoremove)
    tabletofix<-tabletofix[!(tabletofix$time %in% kill$time & tabletofix$code %in% kill$code),]
    
  }	
  
  tabletofix[order(tabletofix$time),]->table
  
  return(table)
}  




# setting wd, and performing initial naming tasks related to file output 
############################
dir<-getwd()
saveithere<-paste(dir,"/checks/",sep="")
  FileisNamed<-right(filename,46)
  
  name<-paste(FileisNamed,"pdf",sep=".")
  name<-gsub("/", ".", name, fixed = TRUE)
  name<-gsub(".csv", "", name, fixed = TRUE)
  name<-gsub("000000_000", "", name, fixed = TRUE)
  name<-gsub("_2000-01-01T","",name, fixed = TRUE)
  name<-gsub("ut_oneatatime.","",name, fixed = TRUE)
  name<-gsub("eatatime.","",name, fixed = TRUE)
  name<-gsub("sis.FixedCSVs.","",name, fixed = TRUE)
  

  name2<-gsub("/", ".", FileisNamed, fixed = TRUE)
  name2<-gsub(".csv", "", name2, fixed = TRUE)
  name2<-gsub("000000_000", "", name2, fixed = TRUE)
  name2<-gsub("_2000-01-01T","",name2, fixed = TRUE)
  name2<-gsub("ut_oneatatime.","",name2, fixed = TRUE)
  name2<-gsub("eatatime.","",name2, fixed = TRUE)
  name2<-gsub("sis.FixedCSVs.","",name2, fixed = TRUE)
  
  
  fullname<-paste(location,name,sep="/")
###########  
#read in table
  read.table(filename,sep=",",header=TRUE)->table

#Order selections according to begin time
	table[order(table$time),]->table

	

#FIRST STEP PREPROCESSING
	
#Early files used "Mpres" instead of "M pres"...these were incorrect and are removed in this step		
	rowswithMpresnospace<-which(table$code == "Mpres")
	if (length(rowswithMpresnospace)!=0){
	table<-table[-((which(table$code == "Mpres"))),] #Removes 'Mpres'
	}
	
####Remove spaces from 'code' identifiers
	table$code<-gsub(" ", "", table$code, fixed = TRUE)

#Order selections according to begin time	
	table[order(table$time),]->table	
	
#ELIMINATES DUPLICATE ROWS	
	table<-unique(table) 
###########
	#DATA CHECKING LOCATION, BEHAVIOR BY BEHAVIOR
	behaviortoinvestigate<-"OPMF"
	info<-table[grep(behaviortoinvestigate, table$code), ]
	rle(table[grep(behaviortoinvestigate, table$code), ][,2])
######################
	
#Make repetitively pressed buttons (e.g. OPMW, OPAW, OPAC3) into duration of RM behaviors (e.g. RM1)

  flappinglist<-c("OPAW","OPAC3","OPAMW")	
	wingsflapping<-subset(table,code %in% flappinglist)	
	wingsflapping<-subset(wingsflapping,!duplicated(time))
	if(nrow(wingsflapping)>0){
	wingsflapping$timetonext<-c(diff(wingsflapping$time),1000)
	wingsflapping$keep<- ifelse((wingsflapping$code==(shift(wingsflapping$code,1L))) | (wingsflapping$code!=(shift(wingsflapping$code,1L)) & wingsflapping$timetonext>0.2),"yes","no")
	wingsflapping<-subset(wingsflapping,keep=="yes")
	wingsflapping<-wingsflapping[,c(1:3)]
	}
	
	
	taillist<-c("OPMT")
	tailflapping<-subset(table,code %in% taillist)	
	
	headlist<-c("OPAH")
	headflapping<-subset(table,code %in% headlist)
	  
	torsolist<-c("OPAC1","OPAC2","OPATW","OPAMTT")
	torsoflapping<-subset(table, code %in% torsolist)
	torsoflapping<-subset(torsoflapping,!duplicated(time))
	if(nrow(torsoflapping)>0){
	torsoflapping$timetonext<-c(diff(torsoflapping$time),1000)
	torsoflapping$keep<- ifelse((torsoflapping$code==(shift(torsoflapping$code,1L))) | (torsoflapping$code!=(shift(torsoflapping$code,1L)) & torsoflapping$timetonext>0.2),"yes","no")
	torsoflapping<-subset(torsoflapping,keep=="yes")
	torsoflapping<-torsoflapping[,c(1:3)]
	}
	
	
	
repeptiveornamentaccentuations<-list(wingsflapping,tailflapping,headflapping,torsoflapping)	

durationofrepetbehavs<-table[1,]	
durationofrepetbehavs<-durationofrepetbehavs[-1,]
RepBehavLabels<-c("RepWing","RepTail","RepHead","RepTorso")

for (RMs in 1:length(repeptiveornamentaccentuations)){	
#flappers<-tailflapping

  newlabel<-RepBehavLabels[RMs]
  flappers<-repeptiveornamentaccentuations[[RMs]]
  flappers<-flappers[order(flappers$time),]
  
  if(nrow(flappers)>1){
    attach(flappers)
    
    timetonext<-c(diff(time),"FALSE")	
  	TFseq<-timetonext<1 
  	NewVec<-vector(length=length(TFseq))
  	
  	for (gg in 1:length(TFseq)) {
  	  
  	  xnow<-TFseq[gg]
  	  pre<-ifelse(gg>1,TFseq[(gg-1)],NA)
  	  nex<-ifelse(gg==length(TFseq),NA,TFseq[(gg+1)])
  	  
  	  if(is.na(pre)){  
      	  if(xnow==TRUE & is.na(pre)){
      	    newvalue<-"start"
      	  } else {
      	    if (xnow==FALSE & is.na(pre)){
      	      newvalue<-"off"
      	    }
      	  }
  	    } else {
  	        if (is.na(nex)){
  	          if (xnow==FALSE & (pre)==FALSE){
  	            newvalue<-"off"
  	          }
  	          if (xnow==FALSE & (pre)==TRUE){
  	            newvalue<-"end"
  	          }
  	        } else {
        	  if(xnow==TRUE & (nex)==TRUE & (pre)==FALSE){
        	    newvalue<-"start"
        	  } 
        	  if (xnow==TRUE & (nex)==FALSE & (pre)==TRUE){
        	    newvalue<-"on"
        	  }
  	        if(xnow==TRUE & (nex)==FALSE & (pre)==FALSE){
        	    newvalue<-"start"
        	  } 
    	      if (xnow==TRUE & (nex)==TRUE & (pre)==TRUE){
    	        newvalue<-"on"
    	      }
        	  if (xnow==FALSE & (nex)==TRUE & (pre)==TRUE){
        	    newvalue<-"end"
        	  }
        	  if (xnow==FALSE & (nex)==TRUE & (pre)==FALSE){
        	    newvalue<-"off"
        	  }
    	      if (xnow==FALSE & (nex)==FALSE & (pre)==FALSE){
    	        newvalue<-"off"
    	      }
    	      if (xnow==FALSE & (nex)==FALSE & (pre)==TRUE){
    	        newvalue<-"end"
    	      }
  	    }
  	  }
  	  NewVec[gg]<-newvalue
  	}
  	detach(flappers)
  	
  	flappers$times<-NewVec
  	startsofrep<-subset(flappers,times=="start")
  	if(nrow(startsofrep)>0){startsofrep$code<-paste(newlabel,"(start)",sep="")}
  	endsofrep<-subset(flappers,times=="end")
  	if(nrow(endsofrep)>0){endsofrep$code<-paste(newlabel,"(end)",sep="")}
  	
  	thesebehaviors<-rbind(startsofrep,endsofrep)
  	thesebehaviors<-thesebehaviors[,-4]
  
  	
  	durationofrepetbehavs<-rbind(durationofrepetbehavs,thesebehaviors)
  	durationofrepetbehavs<-durationofrepetbehavs[order(durationofrepetbehavs$time),]	
  }
}	

table<-rbind(table,durationofrepetbehavs)	
	
	
#SECOND STEP PROCESSING (handling OFF periods)	
	
	
#############################################################	
	#if last OFF (end) is the same as the END, this cuts off analysis at that point
	
	if(nrow(subset(table,code=="OFF(end)"))>0){
	  if(isTRUE(all.equal(max((subset(table,code=="OFF(end)")[,1]),na.rm=TRUE),subset(table,code=="END")[,1],tolerance = 0.0009))){	  
	    OriginalEndingTime<-	max((subset(table,code=="OFF(end)")[,1]),na.rm=TRUE)
	    TimeBirdFlewOff<-tail(subset(table,code=="OFF(start)")[,1],1)
	    FinishEnd<-data.frame("time"=TimeBirdFlewOff,"code"="END","class"=0)
	    beforeflewaway<-subset(table,time<=TimeBirdFlewOff & code !="OFF(start)")
	    previousOFFs<-subset(table,time<TimeBirdFlewOff & code =="OFF(start)")
	    
	    AtEnd<-subset(table,time==OriginalEndingTime & code !="END" & code != "OFF(end)")
	    if(nrow(AtEnd)>0){
	      AtEnd$time<-TimeBirdFlewOff
	    }
	    
	    
	    tabletable<-rbind(beforeflewaway,previousOFFs,AtEnd,FinishEnd)
	    table<-tabletable[order(tabletable$time),]
	  }
	}	
	
	#Finds END time and uses to calculate TotalTime
	TotalTime<-table[which(table$code == "END"), ][,1]	
	
	######if first OFF (start) is at t=0, this cuts off analysis during this first 'missing' bit
	BirdWasOffatBeginning<-0
	TimeBirdFlewON<-0
	if(nrow(subset(table,code=="OFF(start)"))>0){
	  if(isTRUE(all.equal(min((subset(table,code=="OFF(start)")[,1]),na.rm=TRUE),0,tolerance = 0.0009))){	  
	    BirdWasOffatBeginning<-1
	    OriginalStartingTime<-	0
	    TimeBirdFlewON<-subset(table,code=="OFF(end)")[1,1]
	    #FinishEnd<-data.frame("time"=TimeBirdFlewOff,"code"="END","class"=0)
	    afterflewON<-subset(table,time>=TimeBirdFlewON & code !="OFF(end)")
	    subsequentOFFs<-subset(table,time>TimeBirdFlewON & code =="OFF(end)")
	    
	    AtStart<-subset(table,time==0 &  code != "OFF(start)")
	    if(nrow(AtStart)>0){
	      AtStart$time<-TimeBirdFlewON
	    }
	    
	    
	    tabletabletable<-rbind(afterflewON,subsequentOFFs,AtStart)
	    table<-tabletabletable[order(tabletabletable$time),]
	    
	    
	  }
	}	
	
	
	OFFtimes<-table[grep("OFF", table$code), ][,1]
	OFFstartTimes<-subset(table,code=="OFF(start)")[,1]  
	OFFendTimes<-subset(table,code=="OFF(end)")[,1]  
	
#############################################################	
	
#PLOT A	
	

table<-table[order(table$time),]	
table<-fixOFFs(table)	

#Adds O2 for every changed in O3
O2additions<-rbind(table[which(table$code == "O3(start)"),],table[which(table$code == "O3(end)"),])
O2additions$code<-(gsub("O3(end)", "O2", gsub("O3(start)", "O2", O2additions$code, fixed = TRUE), fixed = TRUE)) #substituting
table<-rbind(table,O2additions)	


#THIRD STEP PROCESSING (handling body orientation issues, e.g. O5s and O6s, O4 contingencies)	
#############################################################			
	#############################################################	
  	#############################################################	
		#Finds END time and uses to calculate TotalTime
	TotalTime<-table[which(table$code == "END"), ][,1]
	
	
	#If there are no 05s, then the bird was upright the entire time
	if(nrow(table[which(table$code == "O5(start)"), ])!=0){ 
	  firstO5<-min(table[which(table$code == "O5(start)"), ][,1],na.rm = TRUE) #Takes times of all O5 starts and finds first (minimum)
	  lastO5<-max(table[which(table$code == "O5(end)"), ][,1],na.rm = TRUE) #Takes times of all O5 starts and finds last (maximum)
	 # O5starttable<-rbind(table[which(table$code == "O5(end)"),],table[which(table$code == "O5(start)"),],table[which(table$code == "END"),])#binds all O5 starts, ends, and the Trial End
	  O5starttable<-table[grep("O5", table$code), ]
	  endtime<-table[which(table$code == "END"),][,1]
	  
	  O6starttable<-O5starttable
	 # O6starttable<-O6starttable[-(which(O6starttable$code=="O5(start)" & O6starttable$time %in% OFFendTimes)),]
	 # O6starttable<-O6starttable[-(which(O6starttable$code=="O5(end)" & O6starttable$time %in% OFFstartTimes)),]
	  
	  if (BirdWasOffatBeginning>0){
	    if (isTRUE(all.equal(firstO5,TimeBirdFlewON,tolerance = 0.0009))){
	      if (isTRUE(all.equal(lastO5,endtime,tolerance = 0.0009))) { #If O5 begins and ends a trial, substiutions should not be made on those obs
	        O6starttable<-O6starttable[!(O6starttable$code == "O5(start)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, TimeBirdFlewON, tolerance = 0.0009)))), ]
	        O6starttable<-O6starttable[!(O6starttable$code == "O5(end)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]
	        O6starttable$code<-(gsub("O5(end)", "O6(start)", O6starttable$code, fixed = TRUE)) 
	        O6starttable$code<-(gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE))#substituting O6starts for O5 ends, etc
	        
	        O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
	         
	        
	       
	      } else { #If O5 begins a trial, substitutions should not be made on that first obs, and O6(end) should be at the end
	        O6starttable<-O6starttable[!(O6starttable$code == "O5(start)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, TimeBirdFlewON, tolerance = 0.0009)))), ]
	        O6starttable$code<-(gsub("O5(end)", "O6(start)", O6starttable$code, fixed = TRUE)) 
	        O6starttable$code<-(gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE))#substituting O6starts for O5 ends, etc
	        O6starttable$code<-(gsub("END", "O6(end)", O6starttable$code, fixed = TRUE)) #substituting O6starts for O5 ends, etc
	        
	        O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
	        
	        O6ends<-O6starttable[1,]
	        O6ends$time<-endtime
	        O6ends$code<-"O6(end)"
	        
	        O6starttable<-rbind(O6starttable,O6ends)
	         
	      }
	    } else { #IF O5 doesn't start the trial, O6 must
	      if (isTRUE(all.equal(lastO5,endtime,tolerance = 0.0009))) {
	        O6startsatzero<-O6starttable[1,]
	        O6startsatzero$time<-TimeBirdFlewON
	        O6startsatzero$code<-"O6(start)" 
	        O6starttable<-rbind(O6startsatzero,O6starttable)#This (+3 lines above) adds an O6(start) at t=0
	        
	        O6starttable<-O6starttable[!(O6starttable$code == "O5(end)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]#Don't replace O5(end) if it's at the END
	        O6starttable$code<-(gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE))
	        O6starttable$code<-(gsub("O5(end)", "O6(start)", O6starttable$code, fixed = TRUE))
	        
	        O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
	         
	        
	      } else { #When O5 neither starts the trial nor ends the trial, all subs can be made, END can be turned into O6(end), and O6(start) can be added at t=0
	        O6startsatzero<-O6starttable[1,]
	        O6startsatzero$time<-TimeBirdFlewON
	        O6startsatzero$code<-"O6(start)" 
	        O6starttable<-rbind(O6startsatzero,O6starttable)#This (+3 lines above) adds an O6(start) at t=0
	        
	        O6starttable$code<-(gsub("O5(end)", "O6(start)", gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE), fixed = TRUE)) #substituting O6starts for O5 ends, etc
	        O6starttable$code<-(gsub("END", "O6(end)", O6starttable$code, fixed = TRUE)) #If 05 doesn't end the trial, O6 must
	        
	        O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
	        
	        O6ends<-O6starttable[1,]
	        O6ends$time<-endtime
	        O6ends$code<-"O6(end)"
	        
	        O6starttable<-rbind(O6starttable,O6ends)
	         
	        
	      }
	    }
	  } else {
	  if (isTRUE(all.equal(firstO5,0,tolerance = 0.0009))){
	    if (isTRUE(all.equal(lastO5,endtime,tolerance = 0.0009))) { #If O5 begins and ends a trial, substiutions should not be made on those obs
	      O6starttable<-O6starttable[!(O6starttable$code == "O5(start)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, 0, tolerance = 0.0009)))), ]
	      O6starttable<-O6starttable[!(O6starttable$code == "O5(end)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]
	      O6starttable$code<-(gsub("O5(end)", "O6(start)", O6starttable$code, fixed = TRUE)) 
	      O6starttable$code<-(gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE))#substituting O6starts for O5 ends, etc
	      
	      O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
	       
	      
	    } else { #If O5 begins a trial, substitutions should not be made on that first obs, and O6(end) should be at the end
	      O6starttable<-O6starttable[!(O6starttable$code == "O5(start)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, 0, tolerance = 0.0009)))), ]
	      O6starttable$code<-(gsub("O5(end)", "O6(start)", O6starttable$code, fixed = TRUE)) 
	      O6starttable$code<-(gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE))#substituting O6starts for O5 ends, etc
	      O6starttable$code<-(gsub("END", "O6(end)", O6starttable$code, fixed = TRUE)) #substituting O6starts for O5 ends, etc
	      
	      O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
	      
	      O6ends<-O6starttable[1,]
	      O6ends$time<-endtime
	      O6ends$code<-"O6(end)"
	      
	      O6starttable<-rbind(O6starttable,O6ends)
	       
	    }
	  } else { #IF O5 doesn't start the trial, O6 must
	    if (isTRUE(all.equal(lastO5,endtime,tolerance = 0.0009))) {
	      O6startsatzero<-O6starttable[1,]
	      O6startsatzero$time<-0
	      O6startsatzero$code<-"O6(start)" 
	      O6starttable<-rbind(O6startsatzero,O6starttable)#This (+3 lines above) adds an O6(start) at t=0
	      
	      O6starttable<-O6starttable[!(O6starttable$code == "O5(end)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]#Don't replace O5(end) if it's at the END
	      O6starttable$code<-(gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE))
	      O6starttable$code<-(gsub("O5(end)", "O6(start)", O6starttable$code, fixed = TRUE))
	      
	      O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
	       
	      
	    } else { #When O5 neither starts the trial nor ends the trial, all subs can be made, END can be turned into O6(end), and O6(start) can be added at t=0
	      O6startsatzero<-O6starttable[1,]
	      O6startsatzero$time<-0
	      O6startsatzero$code<-"O6(start)" 
	      O6starttable<-rbind(O6startsatzero,O6starttable)
	      #This (+3 lines above) adds an O6(start) at t=0
	      
	      O6ends<-O6starttable[1,]
	      O6ends$time<-endtime
	      O6ends$code<-"O6(end)"
	      
	      O6starttable$code<-(gsub("O5(end)", "O6(start)", gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE), fixed = TRUE)) #substituting O6starts for O5 ends, etc
	      O6starttable$code<-(gsub("END", "O6(end)", O6starttable$code, fixed = TRUE)) #If 05 doesn't end the trial, O6 must
	      
	      O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
	      
	      O6starttable<-rbind(O6starttable,O6ends)
	       
	      
	      }
	    }
	  }
	} else {
	  if(BirdWasOffatBeginning>0){
	      O6starttable<-data.frame("time"=c(TimeBirdFlewON,TotalTime),"code"=c("O6(start)","O6(end)"),"class"=c(5,5))
	  } else {
	      O6starttable<-data.frame("time"=c(0,TotalTime),"code"=c("O6(start)","O6(end)"),"class"=c(5,5))
	  }
	} 
	
  table2<-rbind(table,O6starttable)
  

	table<-table2 #keeps table our working, modified dataframe
	table<-unique(table) #ELIMINATES DUPLICATE ROWS
	table<-table[order(table$time),]	
	table<-fixOFFs(table)	#Fixes OFF issues

#PLOT D
#PLOTB<-table

	
	#Order selections according to begin time
	table<-table[order(table$time),]
	
	O4active<-rbind(table[which(table$code == "O4(start)"),],table[which(table$code == "O4(end)"),])
	O4active<-O4active[order(O4active$time),]#Order selections according to begin time
	O4starts<-table[which(table$code == "O4(start)"),]
	O4ends<-table[which(table$code == "O4(end)"),]
	O4bouts<-length(table[which(table$code == "O4(start)"),][,1])
	

	compareO5<-table[which(table$code == "O5(start)"),]
	compareO5end<-table[which(table$code == "O5(end)"),]
	compareO6<-table[which(table$code == "O6(start)"),]
	compareO6end<-table[which(table$code == "O6(end)"),]
	
	#This if addresses all postural (O5, O6,O7,O8) issues (if there are at least one or more O4s)
	if(nrow(O4starts)==nrow(O4ends) & O4bouts>0){   #IF BRANCH IS NOT VERTICAL THE ENTIRE TIME, WE NEED THIS CONTINGENCY (ELSE FOR CASES WHEN WE NEVER HAVE ANY O4 PERIODS/BOUTS)

	    #
	    if(nrow(compareO5)==nrow(compareO5end) & nrow(compareO5)!=0 & nrow(compareO5end)!=0) { #ONLY DO THIS IF EQUAL STARTS/STOPS FOR O5 (and !=0)
	      for(j in 1:length(compareO5[,1])){
	        
	        for(i in 1:O4bouts){
	          start<-O4starts[i,1]
	          stop<-O4ends[i,1]
	          startrow<-O4starts[i,]
	          stoprow<-O4ends[i,]
	        
	        if(compareO5$time[j]==start & compareO5end$time[j]<stop ){#S3_This part handles O5s that start exactly as O4 starts and end before each O4 ends.
	          oldendO5<-compareO5end[j,]
	          oldstartO5<-compareO5[j,]
	          
	          O7start<-oldstartO5
	          O7start$code="O7(start)"
	          O7end<-oldendO5
	          O7end$code<-"O7(end)"
	          newO7s<-rbind(O7start,O7end)
	          newO7s$class<-paste(i,j)
	          
	          O7O8s<-newO7s
	          O7O8s$class<-paste(i,j)
	          table<-rbind(table,newO7s) 
	          
	        }  
	        
	        if(compareO5$time[j]<start & compareO5end$time[j]==stop ){#S2_This part handles O5s that start before each O4 starts and end as each O4 ends.
	          oldendO5<-compareO5end[j,]
	          oldstartO5<-compareO5[j,]
	          
	          O8start<-oldstartO5
	          O8start$code="O8(start)"
	          O8end<-startrow
	          O8end$code<-"O8(end)"
	          newO8s<-rbind(O8start,O8end)
	          
	          O7start<-startrow
	          O7start$code="O7(start)"
	          O7end<-stoprow
	          O7end$code<-"O7(end)"
	          newO7s<-rbind(O7start,O7end)
	          newO7s$class<-paste(i,j)
	          
	          O7O8s<-rbind(newO7s,newO8s)
	          O7O8s$class<-paste(i,j)
	          table<-rbind(table,newO7s) #adds new O7s and O8s
	         
	        }  
	        
	        if(compareO5$time[j]<start & compareO5end$time[j]<stop & compareO5end$time[j]>start ){#S1_This part handles all O5s that start before each O4 starts and end while O4 is on
	          oldendO5<-compareO5end[j,]
	          oldstartO5<-compareO5[j,]
	          
	          O8start<-oldstartO5
	          O8start$code="O8(start)"
	          O8end<-startrow
	          O8end$code<-"O8(end)"
	          newO8s<-rbind(O8start,O8end)
	          
	          O7start<-startrow
	          O7start$code="O7(start)"
	          O7end<-oldendO5
	          O7end$code<-"O7(end)"
	          newO7s<-rbind(O7start,O7end)
	          newO7s$class<-paste(i,j)
	          
	          O7O8s<-rbind(newO7s,newO8s)
	          O7O8s$class<-paste(i,j)
	          table<-rbind(table,newO7s) 
	        }  
	        
	        if(compareO5$time[j]>start & compareO5$time[j]<stop & compareO5end$time[j]==stop ){#E4_This part handles all O5s that start while O4 and end exactly as O4 ends
	          oldendO5<-compareO5end[j,]
	          oldstartO5<-compareO5[j,]
	          
	          Nope<-oldstartO5[-1,]
	          newO8s<-Nope
	          
	          O7start<-oldstartO5
	          O7start$code="O7(start)"
	          O7end<-oldendO5
	          O7end$code<-"O7(end)"
	          newO7s<-rbind(O7start,O7end)
	          newO7s$class<-paste(i,j)
	          
	          O7O8s<-rbind(newO7s,newO8s)
	          O7O8s$class<-paste(i,j)
	          table<-rbind(table,newO7s) 
	          
	        }  
	        
	        if(compareO5$time[j]==start & compareO5end$time[j]==stop ){#E3_This part handles all O5s that start and end exactly as O4 starts and ends
	          oldendO5<-compareO5end[j,]
	          oldstartO5<-compareO5[j,]
	          
	          Nope<-oldstartO5[-1,]
	          newO8s<-Nope
	          
	          O7start<-oldstartO5
	          O7start$code="O7(start)"
	          O7end<-oldendO5
	          O7end$code<-"O7(end)"
	          newO7s<-rbind(O7start,O7end)
	          newO7s$class<-paste(i,j)
	          
	          O7O8s<-rbind(newO7s,newO8s)
	          O7O8s$class<-paste(i,j)
	          table<-rbind(table,newO7s) 
	          
	        }  
	        
	        if(compareO5$time[j]==start & compareO5end$time[j]>stop ){#E2_This part handles O5s that start as O4 starts and end after O4 ends.
	          oldendO5<-compareO5end[j,]
	          oldstartO5<-compareO5[j,]
	          
	          O8start<-stoprow
	          O8start$code="O8(start)"
	          O8end<-oldendO5
	          O8end$code<-"O8(end)"
	          newO8s<-rbind(O8start,O8end)
	          
	          O7start<-oldstartO5
	          O7start$code="O7(start)"
	          O7end<-stoprow
	          O7end$code<-"O7(end)"
	          newO7s<-rbind(O7start,O7end)
	          newO7s$class<-paste(i,j)
	          
	          O7O8s<-rbind(newO7s,newO8s)
	          O7O8s$class<-paste(i,j)
	          table<-rbind(table,newO7s) 
	        } 
	        
	        if(compareO5$time[j]>start & compareO5$time[j]<stop & compareO5end$time[j]>stop ){#E1_This part handles O5s that start while O4 is on, then end after O4 is off
	          oldendO5<-compareO5end[j,]
	          oldstartO5<-compareO5[j,]
	          
	          O8start<-stoprow
	          O8start$code="O8(start)"
	          O8end<-oldendO5
	          O8end$code<-"O8(end)"
	          newO8s<-rbind(O8start,O8end)
	          
	          O7start<-oldstartO5
	          O7start$code="O7(start)"
	          O7end<-stoprow
	          O7end$code<-"O7(end)"
	          newO7s<-rbind(O7start,O7end)
	          newO7s$class<-paste(i,j)
	          
	          O7O8s<-rbind(newO7s,newO8s)
	          O7O8s$class<-paste(i,j)
	          table<-rbind(table,newO7s) 
	          
	        }  
	        
	        if (compareO5$time[j]<start & compareO5end$time[j]>stop){ #SPAN_if start happens before bird goes OFF AND stop happens after bird comes back (OFF(end)), cut out middle and add new stop and start
	          oldendO5<-compareO5end[j,]
	          oldstartO5<-compareO5[j,]
	          
	          O8start<-rbind(oldstartO5,stoprow)
	          O8start$code="O8(start)"
	          O8end<-rbind(startrow,oldendO5)
	          O8end$code<-"O8(end)"
	          newO8s<-rbind(O8start,O8end)
	          
	          O7start<-startrow
	          O7start$code="O7(start)"
	          O7end<-stoprow
	          O7end$code<-"O7(end)"
	          newO7s<-rbind(O7start,O7end)
	          newO7s$class<-paste(i,j)
	          
	          O7O8s<-rbind(newO7s,newO8s)
	          O7O8s$class<-paste(i,j)
	          table<-rbind(table,newO7s) 
	          
	        }
	        
	        if (compareO5$time[j]>start & compareO5end$time[j]<stop){ #MID_if start AND stop of a duration behavior occur while bird is OFF, delete
	          oldendO5<-compareO5end[j,]
	          oldstartO5<-compareO5[j,]
	          
	          O7start<-oldstartO5
	          O7start$code="O7(start)"
	          O7end<-oldendO5
	          O7end$code<-"O7(end)"
	          newO7s<-rbind(O7start,O7end)
	          newO7s$class<-paste(i,j)
	          
	          O7O8s<-newO7s
	          O7O8s$class<-paste(i,j)
	          table<-rbind(table,newO7s) 
	        }
	        
	        if ((compareO5$time[j]<start & compareO5end$time[j]<start) | (compareO5$time[j]>stop & compareO5end$time[j]>stop)){ #PRE OR POST_if start AND stop of a duration behavior occur while bird is not on O4
	          oldendO6<-compareO6end[j,]
	          oldstartO6<-compareO6[j,]
	          O7start<-rbind(oldstartO6,oldendO6)
	          O7start$code=c("O7(start)","O7(end)")
	          newO7s<-O7start
	          
	          table<-rbind(table,newO7s)
	        }   
	          
	      }
	    }
	   }
	    if(nrow(compareO6)==nrow(compareO6end) & nrow(compareO6)!=0 & nrow(compareO6end)!=0) { #ONLY DO THIS IF EQUAL STARTS/STOPS FOR O6 (and !=0)
	      for(j in 1:length(compareO6[,1])){
	        
	        for(i in 1:O4bouts){
	          start<-O4starts[i,1]
	          stop<-O4ends[i,1]
	          startrow<-O4starts[i,]
	          stoprow<-O4ends[i,]
	        
	        
	        
	        if(compareO6$time[j]==start & compareO6end$time[j]<stop ){#S3_This part handles O6s that start exactly as O4 starts and end before each O4 ends.
	          oldendO6<-compareO6end[j,]
	          oldstartO6<-compareO6[j,]
	          

	          
	        }  
	        if(compareO6$time[j]<start & compareO6end$time[j]==stop ){#S2_This part handles O6s that start before each O4 starts and end as each O4 ends.
	          oldendO6<-compareO6end[j,]
	          oldstartO6<-compareO6[j,]
	          
	          O7start<-oldstartO6
	          O7start$code="O7(start)"
	          O7end<-startrow
	          O7end$code<-"O7(end)"
	          newO7s<-rbind(O7start,O7end)
	          newO7s$class<-paste(i,j)
	          
	          O8start<-startrow
	          O8start$code="O8(start)"
	          O8end<-stoprow
	          O8end$code<-"O8(end)"
	          newO8s<-rbind(O8start,O8end)
	          
	          O8O7s<-rbind(newO8s,newO7s)
	          O8O7s$class<-paste(i,j)
	          table<-rbind(table,newO7s) #adds new O8s and O7s
	          
	        }  
	        
	        if(compareO6$time[j]<start & compareO6end$time[j]<stop & compareO6end$time[j]>start ){#S1_This part handles all O6s that start before each O4 starts and end while O4 is on
	          oldendO6<-compareO6end[j,]
	          oldstartO6<-compareO6[j,]
	          
	          O7start<-oldstartO6
	          O7start$code="O7(start)"
	          O7end<-startrow
	          O7end$code<-"O7(end)"
	          newO7s<-rbind(O7start,O7end)
	          newO7s$class<-paste(i,j)
	          
	          O8start<-startrow
	          O8start$code="O8(start)"
	          O8end<-oldendO6
	          O8end$code<-"O8(end)"
	          newO8s<-rbind(O8start,O8end)
	          
	          O8O7s<-rbind(newO8s,newO7s)
	          O8O7s$class<-paste(i,j)
	          table<-rbind(table,newO7s) 
	        }  
	        
	        if(compareO6$time[j]>start & compareO6$time[j]<stop & compareO6end$time[j]==stop ){#E4_This part handles all O6s that start while O4 and end exactly as O4 ends
	          oldendO6<-compareO6end[j,]
	          oldstartO6<-compareO6[j,]
	          

	          
	        }  
	        
	        if(compareO6$time[j]==start & compareO6end$time[j]==stop ){#E3_This part handles all O6s that start and end exactly as O4 starts and ends
	          oldendO6<-compareO6end[j,]
	          oldstartO6<-compareO6[j,]
	          
	        
	          
	        }  
	        
	        if(compareO6$time[j]==start & compareO6end$time[j]>stop ){#E2_This part handles O6s that start as O4 starts and end after O4 ends.
	          oldendO6<-compareO6end[j,]
	          oldstartO6<-compareO6[j,]
	          
	          O7start<-stoprow
	          O7start$code="O7(start)"
	          O7end<-oldendO6
	          O7end$code<-"O7(end)"
	          newO7s<-rbind(O7start,O7end)
	          newO7s$class<-paste(i,j)
	          
	          O8start<-oldstartO6
	          O8start$code="O8(start)"
	          O8end<-stoprow
	          O8end$code<-"O8(end)"
	          newO8s<-rbind(O8start,O8end)
	          
	          O8O7s<-rbind(newO8s,newO7s)
	          O8O7s$class<-paste(i,j)
	          table<-rbind(table,newO7s) 
	        } 
	        
	        if(compareO6$time[j]>start & compareO6$time[j]<stop & compareO6end$time[j]>stop ){#E1_This part handles O6s that start while O4 is on, then end after O4 is off
	          oldendO6<-compareO6end[j,]
	          oldstartO6<-compareO6[j,]
	          
	          O7start<-stoprow
	          O7start$code="O7(start)"
	          O7end<-oldendO6
	          O7end$code<-"O7(end)"
	          newO7s<-rbind(O7start,O7end)
	          newO7s$class<-paste(i,j)
	          
	          O8start<-oldstartO6
	          O8start$code="O8(start)"
	          O8end<-stoprow
	          O8end$code<-"O8(end)"
	          newO8s<-rbind(O8start,O8end)
	          
	          O8O7s<-rbind(newO8s,newO7s)
	          O8O7s$class<-paste(i,j)
	          table<-rbind(table,newO7s) 
	          
	        }  
	        
	        if (compareO6$time[j]<start & compareO6end$time[j]>stop){ #SPAN_if start happens before bird goes O4 AND stop happens after O4 goes OFF, cut out middle and add new stop and start
	          oldendO6<-compareO6end[j,]
	          oldstartO6<-compareO6[j,]
	          
	          O7start<-rbind(oldstartO6,stoprow)
	          O7start$code="O7(start)"
	          O7end<-rbind(startrow,oldendO6)
	          O7end$code<-"O7(end)"
	          newO7s<-rbind(O7start,O7end)
	          newO7s$class<-paste(i,j)
	          
	          O8start<-startrow
	          O8start$code="O8(start)"
	          O8end<-stoprow
	          O8end$code<-"O8(end)"
	          newO8s<-rbind(O8start,O8end)
	          
	          O8O7s<-rbind(newO8s,newO7s)
	          O8O7s$class<-paste(i,j)
	          table<-rbind(table,newO7s) 
	          
	        }
	        
	        if (compareO6$time[j]>start & compareO6end$time[j]<stop){ #MID_if start AND stop of a duration behavior occur while bird is OFF, delete
	          oldendO6<-compareO6end[j,]
	          oldstartO6<-compareO6[j,]
	          
	        
	        }
	        
	        if ((compareO6$time[j]<start & compareO6end$time[j]<start) | (compareO6$time[j]>stop & compareO6end$time[j]>stop)){ #PRE OR POST_if start AND stop of O6 occurs completely before or after an O4 bout
	          oldendO6<-compareO6end[j,]
	          oldstartO6<-compareO6[j,]
	          O7start<-rbind(oldstartO6,oldendO6)
	          O7start$code=c("O7(start)","O7(end)")
	          newO7s<-O7start
	          
	          table<-rbind(table,newO7s)
	          }
	        }
	      }
	    } 
	
}	    
	 ################################################
	#This if addresses all postural (O5, O6,O7,O8) issues (if there are no O4s)
	if(nrow(O4starts)==nrow(O4ends) & O4bouts==0){
	  oldendO5<-compareO5end
	  oldstartO5<-compareO5
	  
	  oldendO6<-compareO6end
	  oldstartO6<-compareO6
	  
	  if(nrow(oldstartO5)>0){
	  O8start<-oldstartO5
	  O8start$code="O8(start)"
	  O8end<-oldendO5
	  O8end$code<-"O8(end)"
	  newO8s<-rbind(O8start,O8end)
	  newO8s$class<-5
	  } else {
	    newO8s<-oldendO5
	    if(nrow(newO8s)>0){
	    newO8s$class<-5
	    }
	  }
	  
	  if(nrow(oldstartO6)>0){
	  O7start<-oldstartO6
	  O7start$code="O7(start)"
	  O7end<-compareO6end
	  O7end$code<-"O7(end)"
	  newO7s<-rbind(O7start,O7end)
	  newO7s$class<-5
	  } else {
	    newO7s<-oldstartO6
	    if(nrow(newO7s)>0){
	      newO7s$class<-5
	    }

	  }
	  
	  O7O8s<-rbind(newO7s,newO8s)
	  O7O8s$class<-5
	  table<-rbind(table,O7O8s) 
	}

	midcheck<-table
	table<-unique(table) #ELIMINATES DUPLICATE ROWS
	table<-table[order(table$time),]	
	
	
	#ORDER (start/stop) FIXER FOR DURATION BEHAVIOR
  sevens<-table[grep("O7",table$code),] #O7 values
  if(length(sevens[,1])>0) {sevens$class<-"77"} #Standardize O7 class values
  sevens<-unique(sevens)
  
  
  sevens$keep<- ifelse((sevens$code==(data.table::shift(sevens$code,n=1L,type="lead"))) ,"no","yes")
  sevens<-subset(sevens,keep=="yes" | is.na(keep))[,-4]
  sevens<-subset(sevens,!(is.na(time)))
  
	table<-table[!(grepl("O7", table$code)),] #get rid of O7 valuse
	

  

	  

	  
	
	
	table<-rbind(table,sevens)
	
	
	
	
	table<-fixOFFs(table)	#Fixes OFF issues
	
# 
# 	AllO7<-table[which(table$code=="O7(start)" | table$code=="O7(end)" ),]
# 	AllO7<-AllO7[order(AllO7$time),]	
# 	
# 	AllO8<-table[which(table$code=="O8(start)" | table$code=="O8(end)"),]
# 	AllO8<-AllO8[order(AllO8$time),]	
# 	
# 	noO7O8<-table[-which(table$code=="O7(start)" | table$code=="O7(end)" | table$code=="O8(start)" | table$code=="O8(end)"),]
# 	table<-rbind(noO7O8,AllO7,AllO8)
#########################################################################	
	#If there are no O7s, then the bird was O8 the entire time
	TTolerance <- 0.0009
	
	
	if(nrow(table[which(table$code == "O7(start)"), ])!=0){ 
	  firstO7<-min(table[which(table$code == "O7(start)"), ][,1],na.rm = TRUE) #Takes times of all O7 starts and finds first (minimum)
	  lastO7<-max(table[which(table$code == "O7(end)"), ][,1],na.rm = TRUE) #Takes times of all O7 starts and finds last (maximum)
	  # O7starttable<-rbind(table[which(table$code == "O7(end)"),],table[which(table$code == "O7(start)"),],table[which(table$code == "END"),])#binds all O7 starts, ends, and the Trial End
	  O7starttable<-table[grep("O7", table$code), ]
	  endtime<-table[which(table$code == "END"),][,1]
	  
	  O8starttable<-O7starttable
	  # O8starttable<-O8starttable[-(which(O8starttable$code=="O7(start)" & O8starttable$time %in% OFFendTimes)),]
	  # O8starttable<-O8starttable[-(which(O8starttable$code=="O7(end)" & O8starttable$time %in% OFFstartTimes)),]
	  
	  if (BirdWasOffatBeginning>0){
	    if (isTRUE(all.equal(firstO7,TimeBirdFlewON,tolerance = 0.0009))){
	      if (isTRUE(all.equal(lastO7,endtime,tolerance = 0.0009))) { #If O7 begins and ends a trial, substiutions should not be made on those obs
	        O8starttable<-O8starttable[!(O8starttable$code == "O7(start)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, TimeBirdFlewON, tolerance = 0.0009)))), ]
	        O8starttable<-O8starttable[!(O8starttable$code == "O7(end)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]
	        O8starttable$code<-(gsub("O7(end)", "O8(start)", O8starttable$code, fixed = TRUE)) 
	        O8starttable$code<-(gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE))#substituting O8starts for O7 ends, etc
	        
	        O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
	        
	        
	        
	      } else { #If O7 begins a trial, substitutions should not be made on that first obs, and O8(end) should be at the end
	        O8starttable<-O8starttable[!(O8starttable$code == "O7(start)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, TimeBirdFlewON, tolerance = 0.0009)))), ]
	        O8starttable$code<-(gsub("O7(end)", "O8(start)", O8starttable$code, fixed = TRUE)) 
	        O8starttable$code<-(gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE))#substituting O8starts for O7 ends, etc
	        O8starttable$code<-(gsub("END", "O8(end)", O8starttable$code, fixed = TRUE)) #substituting O8starts for O7 ends, etc
	        
	        O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
	        
	        O8ends<-O8starttable[1,]
	        O8ends$time<-endtime
	        O8ends$code<-"O8(end)"
	        
	        O8starttable<-rbind(O8starttable,O8ends)
	        
	      }
	    } else { #IF O7 doesn't start the trial, O8 must
	      if (isTRUE(all.equal(lastO7,endtime,tolerance = 0.0009))) {
	        O8startsatzero<-O8starttable[1,]
	        O8startsatzero$time<-TimeBirdFlewON
	        O8startsatzero$code<-"O8(start)" 
	        O8starttable<-rbind(O8startsatzero,O8starttable)#This (+3 lines above) adds an O8(start) at t=0
	        
	        O8starttable<-O8starttable[!(O8starttable$code == "O7(end)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]#Don't replace O7(end) if it's at the END
	        O8starttable$code<-(gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE))
	        O8starttable$code<-(gsub("O7(end)", "O8(start)", O8starttable$code, fixed = TRUE))
	        
	        O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
	        
	        
	      } else { #When O7 neither starts the trial nor ends the trial, all subs can be made, END can be turned into O8(end), and O8(start) can be added at t=0
	        O8startsatzero<-O8starttable[1,]
	        O8startsatzero$time<-TimeBirdFlewON
	        O8startsatzero$code<-"O8(start)" 
	        O8starttable<-rbind(O8startsatzero,O8starttable)#This (+3 lines above) adds an O8(start) at t=0
	        
	        O8starttable$code<-(gsub("O7(end)", "O8(start)", gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE), fixed = TRUE)) #substituting O8starts for O7 ends, etc
	        O8starttable$code<-(gsub("END", "O8(end)", O8starttable$code, fixed = TRUE)) #If 05 doesn't end the trial, O8 must
	        
	        O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
	        
	        O8ends<-O8starttable[1,]
	        O8ends$time<-endtime
	        O8ends$code<-"O8(end)"
	        
	        O8starttable<-rbind(O8starttable,O8ends)
	        
	        
	      }
	    }
	  } else {
	    if (isTRUE(all.equal(firstO7,0,tolerance = 0.0009))){
	      if (isTRUE(all.equal(lastO7,endtime,tolerance = 0.0009))) { #If O7 begins and ends a trial, substiutions should not be made on those obs
	        O8starttable<-O8starttable[!(O8starttable$code == "O7(start)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, 0, tolerance = 0.0009)))), ]
	        O8starttable<-O8starttable[!(O8starttable$code == "O7(end)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]
	        O8starttable$code<-(gsub("O7(end)", "O8(start)", O8starttable$code, fixed = TRUE)) 
	        O8starttable$code<-(gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE))#substituting O8starts for O7 ends, etc
	        
	        O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
	        
	        
	      } else { #If O7 begins a trial, substitutions should not be made on that first obs, and O8(end) should be at the end
	        O8starttable<-O8starttable[!((O8starttable$code == "O7(end)" & abs(O8starttable$time - endtime) < TTolerance)), ]
	        O8starttable$code<-(gsub("O7(end)", "O8(start)", O8starttable$code, fixed = TRUE)) 
	        O8starttable$code<-(gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE))#substituting O8starts for O7 ends, etc
	        O8starttable$code<-(gsub("END", "O8(end)", O8starttable$code, fixed = TRUE)) #substituting O8starts for O7 ends, etc
	        
	        O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
	        
	        O8ends<-O8starttable[1,]
	        O8ends$time<-endtime
	        O8ends$code<-"O8(end)"
	        
	        O8starttable<-rbind(O8starttable,O8ends)
	        
	      }
	    } else { #IF O7 doesn't start the trial, O8 must
	      if (isTRUE(all.equal(lastO7,endtime,tolerance = 0.0009))) {
	        O8startsatzero<-O8starttable[1,]
	        O8startsatzero$time<-0
	        O8startsatzero$code<-"O8(start)" 
	        O8starttable<-rbind(O8startsatzero,O8starttable)#This (+3 lines above) adds an O8(start) at t=0
	        
	        O8starttable<-O8starttable[!(O8starttable$code == "O7(end)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]#Don't replace O7(end) if it's at the END
	        O8starttable$code<-(gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE))
	        O8starttable$code<-(gsub("O7(end)", "O8(start)", O8starttable$code, fixed = TRUE))
	        
	        O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
	        
	        
	      } else { #When O7 neither starts the trial nor ends the trial, all subs can be made, END can be turned into O8(end), and O8(start) can be added at t=0
	        O8startsatzero<-O8starttable[1,]
	        O8startsatzero$time<-0
	        O8startsatzero$code<-"O8(start)" 
	        O8starttable<-rbind(O8startsatzero,O8starttable)
	        #This (+3 lines above) adds an O8(start) at t=0
	        
	        O8ends<-O8starttable[1,]
	        O8ends$time<-endtime
	        O8ends$code<-"O8(end)"
	        
	        O8starttable$code<-(gsub("O7(end)", "O8(start)", gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE), fixed = TRUE)) #substituting O8starts for O7 ends, etc
	        O8starttable$code<-(gsub("END", "O8(end)", O8starttable$code, fixed = TRUE)) #If 05 doesn't end the trial, O8 must
	        
	        O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
	        
	        O8starttable<-rbind(O8starttable,O8ends)
	        
	        
	      }
	    }
	  }
	} else {
	  if(BirdWasOffatBeginning>0){
	    O8starttable<-data.frame("time"=c(TimeBirdFlewON,TotalTime),"code"=c("O8(start)","O8(end)"),"class"=c(5,5))
	  } else {
	    O8starttable<-data.frame("time"=c(0,TotalTime),"code"=c("O8(start)","O8(end)"),"class"=c(5,5))
	  }
	} 
	
	
	
	table5<-rbind(table,O8starttable)	

	table<-table5 #keeps table our working, modified dataframe
	table<-unique(table) #ELIMINATES DUPLICATE ROWS
	table<-table[order(table$time),]	
	
	
	#table<-fixOFFs(table)
	
	
	
	#############################################################	
	
##########REMOVES instances where there are multiple *starts* in a row
# 	OrientationDurations<-c("O3","O4","O5","O6","O7","O8")	
# 	
# 	for (OD2 in OrientationDurations){  
# 	  behav.observations <- table[grep(OD2, table$code), ]
# 	  goo<-rle(behav.observations$code)
# 	  
# 	  firstone2<-behav.observations[(which(rep(goo$lengths>1,times=goo$lengths))),]
# 	  firstone<-firstone2[1,]
# 	  alternators<-behav.observations[-(which(rep(goo$lengths>1,times=goo$lengths))),]
# 	  totalsequence<-rbind(firstone,alternators)
# 	  totalsequence<-totalsequence[order(totalsequence$time),]    
# 	  
# 	
# 	  kill<-firstone2[-1,]
# 	  table<-table[!(table$time %in% kill$time & table$code %in% kill$code),]
# 	  
# 	}	
	######################################################
	####################################################
	#ELIMINATE FAULTILY CREATED DEPENDENT VARIABLES (e.g. O7(end) when 08(start) IF O8(start) = 0 )
	EndList<-c("BP1(end)","BP3(end)",
	           "SS1(end)","SS2(end)",
	           "O3(end)","O4(end)",
	           "O5(end)","OFF(end)",
	           "O6(end)","O7(end)","O8(end)",
	           "OPMH.D(end)","OPMB.D(end)","OPMF.D(end)","OPMW.D(end)","OPMT.D(end)",
	           "RM1(end)","RM2(end)","RM3(end)","RM4(end)",
	           "MO2(end)","PU1(end)")
	
	for(xx in EndList){ #Nice little function that eliminates any "end" variables if they occur at t=0
	  times2<-subset(table,code==xx)[,1]
	  if(length(times2)>0){ 
	    if(min(times2)==0){
	      table<-table[-which(table$code == xx & table$time==0),]
	    }
	  }
	}
	
	for(xx in EndList){ #Nice little function that eliminates any "end" variables if they occur at the moment the bird flies on screen for the first time
	  times<-subset(table,code==xx)[,1]
	  if(length(times)>0){    
	    if(min(times)==TimeBirdFlewON){
	      table<-table[-which(table$code == xx & table$time==TimeBirdFlewON),]
	    }
	  }
	}	

StartList2<-c("BP1(start)","BP3(start)",
	           "SS1(start)","SS2(start)",
	           "O3(start)","O4(start)",
	           "O5(start)","OFF(end)",
	           "O6(start)","O7(start)","O8(start)",
	           "OPMH.D(start)","OPMB.D(start)","OPMF.D(start)","OPMW.D(start)","OPMT.D(start)",
	           "RM1(start)","RM2(start)","RM3(start)","RM4(start)",
	           "MO2(start)","PU1(start)")

StartList<-c("BP1(start)","BP3(start)",
              "SS1(start)","SS2(start)",
              "O3(start)","O4(start)",
              "O5(start)",
              "O6(start)","O7(start)","O8(start)",
              "OPMH.D(start)","OPMB.D(start)","OPMF.D(start)","OPMW.D(start)","OPMT.D(start)",
              "RM1(start)","RM2(start)","RM3(start)","RM4(start)",
              "MO2(start)","PU1(start)")
		
table<-fixOFFs(table)	#Fixes OFF issues


for(xj in StartList){ #Nice little function that eliminates any "start" variables if they occur at the moment the bird flies off screen
	  starttimes<-subset(table,code==xj)[,1]
	  flyofftimes<-subset(table,code=="OFF(start)")[,1]
	  if(length(starttimes)>0 & length(flyofftimes)){   
	    table<-table[!(table$time %in% flyofftimes & table$code==xj),]        
	 	  }
	}
	############################################
	###############################
secondcheck<-table

  	
	#############################################################  
	#OFF PERIODS---HOW MANY AND HOW MANY BEHAVIORS WERE MEASURED WHILE "OFF"
	
	happenedwhileoff<-data.frame()
	offstarts<-nrow(subset(table,code=="OFF(start)"))
	offstops<-nrow(subset(table,code=="OFF(end)"))
	off.ss.eq<-(offstarts==offstops)
	for (v in 1:nrow(subset(table,code=="OFF(start)"))){
	  stopBegin<-subset(table,code=="OFF(start)")[v,1]
	  stopEnd<-subset(table,code=="OFF(end)")[v,1]
	  
	  gap.wise.misses<-table[which(table$time < stopEnd & table$time > stopBegin),]
	  happenedwhileoff<-rbind(happenedwhileoff,gap.wise.misses)
	}
	
	happenedwhileoff<-paste(happenedwhileoff$time,collapse="|")
	# 	
	#############################################################	
	table<-unique(table) #ELIMINATES DUPLICATE ROWS	
	


	
	
	
	#FOURTH DATA PROCESSING STEP (removes starts/stops of same behavior, at the same time)
	#These are induced if a duration-based orienatation spans an O4 switch, but if they do not change when O4 changes, they can be eliminated here	
	#############################################################	
	#############################################################	
	






#############################################################################
###################################################################


#DATA SUMMARIZING STAGE BEGINS HERE
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	

#Summary measurements extracted from BOPLog output	
is.odd <- function(x) x %% 2 == 1 #defines is.odd function used to ensure even number of observations for state behaviors





#creates empty dataframe that duration behavior summary info goes into
summarydurationbehaviors<-data.frame(matrix(vector(),0,4,dimnames=list(c(),c("Behavior","Frequency","Duration","Bouts"))),stringsAsFactors=F)
#Creates variable to be filled with duration behaviors with odd numbers
odd.behaviors<-c()

#Loops through all duration behaviors identified above and calculates duration, n, and mean time/bout
for (m in DurationBehaviors){      
  behav.obs <- table[grep(m, table$code), ]
  behav.obs$time<-as.numeric(behav.obs$time)
   if(m=="SP"){
     
     SP1<- subset(behav.obs,code=="SP1")
     SP2<- subset(behav.obs,code=="SP2")   
     SP3<- subset(behav.obs,code=="SP3")  
     SP4<- subset(behav.obs,code=="SP4")
     
     if(nrow(SP1)!=0){
     SP1.SP2.duration<-0  
     for(i in 1:nrow(SP1)){
       SP1.SP2.duration<-SP1.SP2.duration+(SP2$time[i]-SP1$time[i])
     }
     SP1.SP2.freq<-nrow(SP1)
     SP1.SP2.bouts<-  SP1.SP2.duration/ SP1.SP2.freq
     SP1.SP2<-cbind("SP1.SP2",SP1.SP2.freq,SP1.SP2.duration,SP1.SP2.bouts)
     } else {
        SP1.SP2<-data.frame("SP1.SP2",0,0,0)
       }
     colnames(SP1.SP2)<-c("Behavior","Frequency","Duration","Bouts")
     
     if(nrow(SP2)!=0){
     SP2.SP3.duration<-0  
     for(i in 1:nrow(SP2)){
       SP2.SP3.duration<-SP2.SP3.duration+(SP3$time[i]-SP2$time[i])
     }
     SP2.SP3.freq<-nrow(SP2)
     SP2.SP3.bouts<-  SP2.SP3.duration/ SP2.SP3.freq
     SP2.SP3<-cbind("SP2.SP3",SP2.SP3.freq,SP2.SP3.duration,SP2.SP3.bouts)
     } else {
       SP2.SP3<-data.frame("SP2.SP3",0,0,0)
       }
     colnames(SP2.SP3)<-c("Behavior","Frequency","Duration","Bouts")
     
     if(nrow(SP3)!=0){
     SP3.SP4.duration<-0  
     for(i in 1:nrow(SP3)){
       SP3.SP4.duration<-SP3.SP4.duration+(SP4$time[i]-SP3$time[i])
     }
     SP3.SP4.freq<-nrow(SP3)
     SP3.SP4.bouts<-  SP3.SP4.duration/ SP3.SP4.freq
     SP3.SP4<-cbind("SP3.SP4",SP3.SP4.freq,SP3.SP4.duration,SP3.SP4.bouts)
     } else {
       SP3.SP4<-data.frame("SP3.SP4",0,0,0)
        }
     colnames(SP3.SP4)<-c("Behavior","Frequency","Duration","Bouts")
     
   } else {
    if(is.odd(nrow(behav.obs))) {
       behav.obs ##prints table of behaviors if n is odd, suggesting an unclosed duration. 
       odd.behaviors<-c(odd.behaviors,m)
       behav.freq<-"ODD"           #Assigns ODD if behav is mising
       behav.duration<-"ODD"
       behav.bouts<- "ODD" 
    }  else {
        if(nrow(behav.obs)==0) {#Assigns zeros if behav is mising,
          behav.freq<-0
          behav.duration<-0
          behav.bouts<-  0 
      } else {#Otherwise returns durations (behav.duration), frequency (behav.freq), and avg bout length (behav.bouts)
        starts<-behav.obs[grep("start", behav.obs$code), ]
        ends<-behav.obs[grep("end", behav.obs$code), ]
        duration<-0
        for(n in 1:nrow(starts)){
          duration<-duration+(ends$time[n]-starts$time[n])
        }# for every row in 'starts', which has all starts for a behavior, it calculates the difference between end and start and then adds this to 'duration'
        behav.freq<-nrow(starts)
        behav.duration<-duration
        behav.bouts<-  behav.duration/behav.freq
            }
    }
   }
  summary<-cbind(m,behav.freq,behav.duration,behav.bouts)
  colnames(summary)<-c("Behavior","Frequency","Duration","Bouts")
  summarydurationbehaviors<-rbind(summarydurationbehaviors,summary)
  
}
summarydurationbehaviors<-rbind(summarydurationbehaviors,SP1.SP2)
summarydurationbehaviors<-rbind(summarydurationbehaviors,SP2.SP3)
summarydurationbehaviors<-rbind(summarydurationbehaviors,SP3.SP4)

rownames(summarydurationbehaviors)<-summarydurationbehaviors$Behavior
summarydurationbehaviors<-summarydurationbehaviors[,-1]
summarydurationbehaviors$Frequency<-(as.numeric(as.character(summarydurationbehaviors$Frequency)))
summarydurationbehaviors$Duration<-(as.numeric(as.character(summarydurationbehaviors$Duration)))
summarydurationbehaviors$Bouts<-(as.numeric(as.character(summarydurationbehaviors$Bouts)))
summarydurationbehaviors$ids<-rownames(summarydurationbehaviors)

rowwisedurations<-dcast(melt(summarydurationbehaviors, id.var="ids"), 1~ids+variable)
rowwisedurations<-rowwisedurations[,-1]
######################
#EVENT BEHAVIORS
summaryeventehaviors<-data.frame(matrix(vector(),0,2,dimnames=list(c(),c("Behavior","Number"))),stringsAsFactors=F)
for (m in EventBehaviors){       #Loop through each event behavior and counts its occurences
  eventtable<-table[table$code %in% EventBehaviors,]
  #use<-eventtable[jj %in% eventtable$code, ]
  behav.obs <- subset(eventtable,code==m)
  behav.obs$time<-as.numeric(behav.obs$time)
  number<-nrow(behav.obs)
  summary<-cbind(m,number)
  summaryeventehaviors<-rbind(summaryeventehaviors,summary)
}
rownames(summaryeventehaviors)<-summaryeventehaviors$m #applies rownames of behaviors
summaryeventehaviors$number<-as.numeric(as.character(summaryeventehaviors$number)) #
rowwiseevents<-as.data.frame(t(summaryeventehaviors)) #transpose function makes numbers lose their numeric classification
rowwiseevents<-rowwiseevents[-1,]
rowwiseevents<-as.data.frame(t(as.data.frame(sapply(rowwiseevents, function(x) as.numeric(as.character(x))))))  #this function, applied to every column, returns numeric values
#############################################
#COMBINED STATE AND EVENT BEHAVIOR SUMMARY DATA

stateANDduration<-cbind(rowwisedurations,rowwiseevents)
rownames(stateANDduration)<-NULL

#################
#Additional summary measures
#####################################
if (BirdWasOffatBeginning>0){
  accuratetime<-table
  accuratetime$time<-accuratetime$time-TimeBirdFlewON
  TotalTime<-accuratetime[which(accuratetime$code == "END"), ][,1]
  
  table$time<-table$time-TimeBirdFlewON+0.1
}


  TotalTimeWatched<-TotalTime-stateANDduration$OFF_Duration
  PropTimeMoving<-stateANDduration$BP1_Duration/TotalTimeWatched

  Nfemales<-length(table[which(table$code == "Fpres"),][,1])
	Nfemales.disp<-length(table[which(table$code == "Fdisp"),][,1])
	Nmales<-length(table[which(table$code == "Mpres"),][,1])
	
	Time.upright <- TotalTimeWatched-stateANDduration$O3_Duration
	Time.inverted<- stateANDduration$O3_Duration
	
	
	O1.obs <- subset(table,code=="O1")
	O2.obs <- subset(table,code=="O2")
	timespent<-0
	femaleorientedbouts<-nrow(O1.obs)#Number of bouts oriented towards females
	if(nrow(O1.obs)!=0){ #IF there are is at least 1 obs of O1, then we can calculate time oriented towards female
	  
	 
      if(nrow(O2.obs)!=0){ #IF there are O1s and O2s, we can use O2s to calculate end of O1s
	  
	  
	      for(j in 1:femaleorientedbouts){
	        lookingstart<-as.numeric(as.character(O1.obs[j,"time"]))
	        alltimesafterorientingtowardsfemale<-subset(O2.obs,time>lookingstart)
	        lookingend<-0
	    
	      if(nrow(alltimesafterorientingtowardsfemale)>0){
	        lookingend<-as.numeric(as.character(min(alltimesafterorientingtowardsfemale$time,na.rm=TRUE)))
	        } else {	lookingend<-TotalTime
	                }
	        
	      timespent<-timespent+(lookingend-lookingstart)
	       }
      } else {#If there is an O1 and no subsequent O2s, then we use the end time (TotalTime) to calc duration spent oriented towards female
        timespent<-TotalTime-as.numeric(as.character(O1.obs[1,"time"]))
      }
	}

	Time.orientedtofemale <- timespent
	

	lengthodds<-length(odd.behaviors)

	
	odd.behaviors<-paste(odd.behaviors,collapse="-")
	
	behaviorsscored<-nrow(table)
	
	belongsto<-species.name
	
	unique.behaviors.scored<-as.data.frame(unique(table$code))
	colnames(unique.behaviors.scored)<-"uniques"
	n.unique.behavs<-length(unique.behaviors.scored[ grep("end", unique.behaviors.scored$uniques, invert = TRUE,ignore.case = TRUE) , ])
	
	
#  printme<-data.frame(FileisNamed,belongsto,TotalTimeWatched,behaviorsscored,n.unique.behavs,Nmales,Nfemales,Nfemales.disp,Time.upright,Time.inverted,Time.orientedtofemale,femaleorientedbouts,stateANDduration,odd.behaviors,off.ss.eq,happenedwhileoff)
	
  origsourcefolder<-getwd()
######################################################################################################################################
#####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  
  
  compositebehavioralmetrictable<-table[table$code!="END" & table$code!="Mpres" & table$code!="Fpres"&
                                         # table$code!="OFF(start)"& table$code!="OFF(end)" &
                                          table$code!="Fdisp",]
  
  compositebehavioralmetrictable$time<-round(compositebehavioralmetrictable$time,digit=1)
  
  #adds 0.1 sec to OFF(starts) to avoid overlap with behaviors that end as birds goes OFF
  compositebehavioralmetrictable[compositebehavioralmetrictable$code=="OFF(start)","time"]<-compositebehavioralmetrictable[compositebehavioralmetrictable$code=="OFF(start)","time"]+0.1
  
  #subtracts 0.1 sec to OFF(end)s to avoid overlap with behaviors that start as bird comes back on
  compositebehavioralmetrictable[compositebehavioralmetrictable$code=="OFF(end)","time"]<-compositebehavioralmetrictable[compositebehavioralmetrictable$code=="OFF(end)","time"]-0.1
  
  newlength<-round(TotalTime/0.1)+2 #Complete # of 0.1 sec intervals for the video (adds 1, b/c we add 0.1 to the first behavior so it didn't occur at 0 = 0th column)
#  allbehaviorstoconverttenthsecondintervals<-cbind(StartList2,EndList)
  
  
  AllBehaviors<-c("O5","O6","O3","O4",
                  "BP1","BP2","BP3",
                  "SS1","SS2",
                  "OFF", "O2",
                  
                  #"O1","O2",
                  "OPMH","OPMB","OPMF","OPMW","OPMT",
                  "OPAH","OPAC1","OPAC2","OPAC3","OPAW","OPATW","OPAMW","OPAMTT",
                  "OPABH","OPABB","OPABF","OPABW","OPABT",
                  "PC1","PC2","PC3","PC4",
                  #"RepWing","RepTail","RepHead","RepTorso",
                  "RM1","RM2","RM3","RM4",
                  "MO1","MO2","PU1",
                  "SP1","SP2","SP3","SP4")
  
  
  
  allbehaviorstolatercombine<-length(AllBehaviors)
  completebehaviormatrix<-matrix(data=NA,nrow=allbehaviorstolatercombine,ncol=newlength)
  row.names(completebehaviormatrix)<-AllBehaviors
  
  DurationBehaviors.nosp.no078<-c("O5","O6","O3","O4",
                                  "BP1","BP3","OFF","SS1","SS2",
                                  "OPMH.D","OPMB.D","OPMF.D","OPMW.D","OPMT.D",
                                  #"RepWing","RepTail","RepHead","RepTorso",
                                  #"RM1","RM2","RM3","RM4",
                                  "MO2","PU1")
  
  
  
#Special class (0,1,2,3) determines which specials are treated as duration behaviors  
  if(specialclass==0){
    DurationBehaviors.toUse<-DurationBehaviors.nosp.no078
  }
  if(specialclass==1){
    DurationBehaviors.toUse<-c(DurationBehaviors.nosp.no078,"SP3")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP3"]<-"SP3(start)")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP4"]<-"SP3(end)")
  }
  if(specialclass==2){
    DurationBehaviors.toUse<-c(DurationBehaviors.nosp.no078,"SP1")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP1"]<-"SP1(start)")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP2"]<-"SP1(end)")
  }
  if(specialclass==3){
    DurationBehaviors.toUse<-c(DurationBehaviors.nosp.no078,c("SP1","SP3"))
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP1"]<-"SP1(start)")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP2"]<-"SP1(end)")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP3"]<-"SP3(start)")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP4"]<-"SP3(end)")
    

  }
  
  for (cb in 1:length(DurationBehaviors.toUse)){ #loop through all initially scored behaviors, but not Specials
      
      thisone<-DurationBehaviors.toUse[cb] #gives new name, "thisone" to the current behavior
      
      start.code<-paste(thisone,"(start)",sep="")
      end.code<-paste(thisone,"(end)",sep="")
      
      starts<-compositebehavioralmetrictable[compositebehavioralmetrictable$code==start.code, ] #df of starts
      ends<-compositebehavioralmetrictable[compositebehavioralmetrictable$code==end.code, ] #df of ends
      
      number.of.starts<-length(starts[,2])
    
    if(number.of.starts!=0){  #Only run the "fill-in loop" if the behavior was recorded at all    
      for (ooo in 1:number.of.starts)  { #This loop fills the matrix with 1s for every 1/10 sec where the behavior is occuring
       
       columnstofill<-10* (seq(starts[ooo,1],ends[ooo,1],by=0.1)   )

           if(grepl(".D",thisone)){
             thisone<-substr(thisone,1,nchar(thisone)-2)
           } 
       
       completebehaviormatrix[rownames(completebehaviormatrix)==thisone,columnstofill]<-1 #fills columns while behavior is active
       
       
      }
    }
  }  
  
  
  
  EventBehaviors.Fill<-c("BP2","O1","O2",
                         "OPMH","OPMB","OPMF","OPMW","OPMT","OPAH",
                         "OPAC1","OPAC2","OPAC3","OPAW","OPATW","OPAMW","OPAMTT",
                         "OPABH","OPABB","OPABF","OPABW","OPABT",
                         "PC1","PC2","PC3","PC4",
                         "MO1","SP1","SP2","SP3","SP4")
  
  event.length<-length(EventBehaviors)
  
  for (ebs in 1:event.length){
    
    behaviortopull<-EventBehaviors.Fill[ebs]
    behavs.to.fill<-compositebehavioralmetrictable[compositebehavioralmetrictable$code==behaviortopull, ]
    column.to.add<-10*behavs.to.fill$time
    column.to.add<-ifelse(column.to.add==0,1,column.to.add)
    
    completebehaviormatrix[rownames(completebehaviormatrix)==behaviortopull,column.to.add]<-1 #fills columns while behavior is active
    
  }
  

  
  compositebehaviorcategorizer<-function(dataframe,x){
    x2<-dataframe[,x,drop=FALSE]
    yyy<-as.data.frame(x2[complete.cases(x2),])
    behav<-paste(rownames(yyy),sep=".",collapse = ".")
    
    }  
  
temporalbehaviorcategories<-matrix(data=NA,nrow=1,ncol=newlength)
    
 for(pdq in 1:newlength){
   temporalbehaviorcategories[1,pdq]<-compositebehaviorcategorizer(completebehaviormatrix,pdq)
 }
  

  
as.data.frame(unique(temporalbehaviorcategories[1,]))
  

no.OFFs<-as.matrix(temporalbehaviorcategories[,!grepl("OFF",temporalbehaviorcategories[1,])])#gets rid of behavior combos with OFF
no.OFFs[no.OFFs == ""] <- NA
no.OFFs<-na.omit(no.OFFs)



###########################################################################################################
#Finds behaviors missing a body-axis (O5/O6) and adds axis from nearest (time-wise) axis
deficient.behaviors<-no.OFFs[((!grepl("O5", no.OFFs[,1], fixed=TRUE)) & (!grepl("O6", no.OFFs[,1], fixed=TRUE)) ),]
containing.behaviors<-no.OFFs[((grepl("O5", no.OFFs[,1], fixed=TRUE)) | (grepl("O6", no.OFFs[,1], fixed=TRUE)) ),]

index.of.deficient<-which((!grepl("O5", no.OFFs[,1], fixed=TRUE)) & (!grepl("O6", no.OFFs[,1], fixed=TRUE)) )
index.of.containing<-which((grepl("O5", no.OFFs[,1], fixed=TRUE)) | (grepl("O6", no.OFFs[,1], fixed=TRUE)) )

####
#Little loop through all deficient (posture-absent) behaviors, adding posture from nearest time
if(length(index.of.deficient)>0){
  for(defic in 1:length(index.of.deficient)){
    defbehav<-deficient.behaviors[defic]
    defind<-index.of.deficient[defic]
    
    replacementbehavior<-no.OFFs[which.min(abs(index.of.containing-defind)),]
    replacementposture<-ifelse(grepl("O5",replacementbehavior),"O5","O6")
    #replacementpostureDOT<-paste(replacementposture,".", sep ="")
    
    
    
    newbehaviorvalue<- paste(replacementposture,defbehav,sep=".")
   
    
  

    no.OFFs[defind,1]<-newbehaviorvalue

  }
}  
####################################################
#######################################################
#no.OFFs<-as.matrix(gsub("O5.O6","O5",no.OFFs[,1]))#swaps out impossible O5/O6

#need to deal with O5.O6s (maybe randomly choose one?)
no.OFFs[((grepl("O5.O6", no.OFFs[,1], fixed=TRUE)) ),]
simultaneousO5O6index<-grep("O5.O6", no.OFFs[,1], fixed=TRUE,value=FALSE)
fixerO5O6er<-c("O5","O6")

if(length(simultaneousO5O6index)>0){
  for(abc in 1:length(simultaneousO5O6index)){
    no.OFFs[simultaneousO5O6index[abc],1]<-as.matrix(gsub("O5.O6",fixerO5O6er[sample(1:2,1)],no.OFFs[simultaneousO5O6index[abc],1]))
  }
} #randomly replaces any simultatneous O5.O6 with either O5 or O6



#no.OFFs<-as.matrix(gsub("O5.O6",fixerO5O6er[sample(1:2,1)],no.OFFs[,1]))#swaps out impossible O5/O6


#Counts occurences for each ID type
summarytable<-table(no.OFFs)
behaviorcounts<-table(no.OFFs)
proportiontable<-prop.table(summarytable)
uniquebehaviors<-length(summarytable)

six<-no.OFFs
#ALL BEHAVIORS and TRANSITIONS
if(uniquebehaviors>1){

##This loop creates the denominator for the calculation of the inverse Simpson Index for behavioral diversity
val2<-0
for (u in 1:uniquebehaviors){
  proptosquare<-proportiontable[u]
  val <-proptosquare^2.0
  val2<-val+val2
}			
#inverse Simpson Index for behavioral diversity
iSimpson<-1.0/as.numeric(val2)
relativebehavioraldiversity<-iSimpson/(uniquebehaviors) #relative color diversity, accounting for # of color classes

##########################################################################################################
#uses "createSequenceMatrix" function (from *markovchain*) to calculate transition matrix
c.no.OFFs<-as.character(no.OFFs)


TransitionMatrix.markov<-createSequenceMatrix(c.no.OFFs,toRowProbs = TRUE,sanitize=FALSE)
TransitionMatrix.cum<-createSequenceMatrix(c.no.OFFs,sanitize=FALSE)

#TransitionMatrix<-TransitionMatrix[-ncol(TransitionMatrix),-ncol(TransitionMatrix)]


time.value<-sum(TransitionMatrix.cum)/10 #seconds watched


noselfs<-TransitionMatrix.cum
diag(noselfs)<-NA

prop.table.excludeNAs<- function(x) {x/sum(x, na.rm=TRUE)}
proportiontable.trans.Simpson<-prop.table.excludeNAs(noselfs) #get rid of diagonal/self-transitions for Simpson Analyses

proportiontable.trans<-prop.table.excludeNAs(TransitionMatrix.cum) 


if(uniquebehaviors>1){
  val3<-matrix(nrow=(((uniquebehaviors^2)-uniquebehaviors)/2),ncol=2)
  rownum<-1
  for (u in 1:uniquebehaviors){
    for(v in u:uniquebehaviors){#iteratively reduces columns analyzed in next loop to avoid double counting transitions 
      if (u!=v){
        val3[rownum,2]<-as.numeric(proportiontable.trans.Simpson[u,v]+proportiontable.trans.Simpson[v,u])
        val3[rownum,1]<-paste(row.names(proportiontable.trans.Simpson)[u],colnames(proportiontable.trans.Simpson)[v],sep="-")
        rownum<-rownum+1
      }
    }
  }
  val4<-as.data.frame(val3)
  val4[,2]<-as.numeric(as.character(val4[,2])) 

  simpsonT<-0

  for(z in 1:nrow(val4)){
  propsq<-val4[z,2]^2.0
  simpsonT<-simpsonT+propsq
  }


iSimpsonT<-1.0/simpsonT #Calculates inverse Simpson index for horizontal transitions
relativetransitiondiversity<-iSimpsonT/(uniquebehaviors*(uniquebehaviors-1)/2)  
} else {
 
  iSimpsonT <- NA
  relativetransitiondiversity<- NA
}


PropSelfTransitioning<-sum(diag(TransitionMatrix.cum))/sum(TransitionMatrix.cum)

network.TransitionMatrix.cum<-network(proportiontable.trans,directed=TRUE,loops=TRUE,hyper=FALSE,matrix.type="adjacency")



showme<-graph_from_adjacency_matrix(proportiontable.trans,mode="directed",weighted=TRUE,diag=FALSE,add.colnames=TRUE)
showme.undirected<-graph_from_adjacency_matrix(proportiontable.trans,mode="undirected",weighted=TRUE,diag=FALSE,add.colnames=TRUE)

showme.looped<- graph_from_adjacency_matrix(proportiontable.trans,mode="directed",weighted=TRUE,diag=TRUE,add.colnames=TRUE)
showme.looped.undirected<- graph_from_adjacency_matrix(proportiontable.trans,mode="undirected",weighted=TRUE,diag=TRUE,add.colnames=TRUE)



E(showme)$width <- 1+E(showme)$weight/50
V(showme)$nombres<-colnames(proportiontable.trans)

E(showme.looped)$width <- 1+E(showme.looped)$weight/50
V(showme.looped)$nombres<-colnames(proportiontable.trans)

# Compute approx self-transition and use that to set node size:
selfsizesscaled<-(diag(proportiontable.trans))*100+1
V(showme)$selfsize<-sqrt(selfsizesscaled)
V(showme)$size <- ((selfsizesscaled)^2)*9

V(showme.looped)$selfsize<-log(selfsizesscaled)+1



#####################################################
# Network Clustering
# l <- layout.fruchterman.reingold(showme.looped)
# l <- layout_with_kk(showme.looped)
# l <- layout_with_lgl(showme.looped)
l <- layout_with_graphopt(showme.looped)
l2 <-layout_with_graphopt(showme)
#l <- layout_nicely(showme.looped)
#behaviorplotlist


#clp <- cluster_label_prop(showme.looped)
clp <-cluster_label_prop(showme.looped.undirected) # input graph should be undirected to make sense.
clp2<-cluster_label_prop(showme.undirected)
# setwd("C:/Users/Rusty/Amazon Drive/BOP/Video_analysis")
# source("PictoGrammerMega_Dec13ab.R", chdir = F)
# #PictoGrammer<-function(proportiontable,hangingbird)
#




################################################################################################################################
#PictoGrammer(proportiontable.trans,hangingbird,uniquebehaviors)

origwd<-dir
setwd(dir<-"C:/Users/Rusty/Amazon Drive/BOP/Video_analysis/Pictograms")

pictogram.names<-list.files(dir,full.names=FALSE,pattern=".png")

proportiontable<-proportiontable
hangingbird<-hangingbird


selflooppoints<-read.table("SelfLoop.csv",sep=",",header=F)

#Creates empty list, then loops through pictogram directory, pulling in each PG, adding it to the list, and naming it
piclist<-list()
for (pg in 1:56){
  img <- readPNG(pictogram.names[pg])  
  
  if(pg>23 && pg <42) {
    w <- matrix(rgb(img[,,1],img[,,2],img[,,3], img[,,4] * 0.75), nrow=dim(img)[1]) #0.5 is alpha
    rimg <- as.raster(w) # raster multilayer object
  } else {
    rimg <- as.raster(img) # raster multilayer object  
  }
  pgname<-gsub('.{4}$', '',pictogram.names[pg])
  piclist[[pg]]<-rimg
  names(piclist)[[pg]]<-pgname
}

#Loops through each unique behavior (from proportiontable) and creates a new, unique, behavioral pictogram
behaviorplotlist<-list()
par(bg="transparent")
plot.new()
goo=1
#pdf("Behav30-trans-fix.pdf",width= 12, height= 12,family="NimbusRom")
#op <- par(mfrow=c(6,6))# rows, columns
#par(mar=c(1,1,1,1)) 


colorweights<-E(showme.looped)$weight/min(E(showme.looped)$weight) 
colorweights<-log(colorweights+0.01)+1

Lab.palette<-colorRampPalette(c("blue",  "red"),space = "Lab")
#  Lab.palette<-colorRampPalette(c("#ffeda0","#feb24c","#f03b20"),space = "Lab")


E(showme.looped)$color<-Lab.palette(max(colorweights))[colorweights]

col.tr <- grDevices::adjustcolor(E(showme.looped)$color, alpha=1) #adds translucency to arrows

SELFCOLORS<-log(diag(proportiontable.trans)/min(E(showme.looped)$weight))+1 
SELFCOLORS.colors<-SELFCOLORS
#Little loop that adds same heatmap scaled colors to self-transitions, for use in behavior boxes
for(ggg in 1:length(SELFCOLORS)){
  selfvalue<-SELFCOLORS[ggg]
  if(selfvalue>0){
    SELFCOLORS.colors[ggg]<-Lab.palette(max(colorweights))[selfvalue]
  } else {
    SELFCOLORS.colors[ggg]<-NA
  }
  
}


for (ub in 1:uniquebehaviors){
  par(bg="transparent")
  graphics::plot(goo, type ="l", xlab="", ylab="", xlim=c(0, 100), ylim=c(0, 100),axes=FALSE)
  
  #53
  if( (length(grep("SS1", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`SS1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #54
  if( (length(grep("SS2", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`SS2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  ##########################################################  
  #2
  if( (length(grep("BP1", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`BP1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #1
  if( (length(grep("BP1.BP2", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`BP1-BP2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #3
  if( (length(grep("BP3", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`BP3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #4
  if( (length(grep("O1", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`O1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #5
  if( (length(grep("O2", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
    rasterImage(piclist$`O2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  ###############################################
  
  #If hangingbird  
  if(hangingbird==1){  
    #6
    if( (length(grep("O5.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`O5-O3-O4-MO1or2a`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #7
    if( (length(grep("O5.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
      rasterImage(piclist$`O5-O3-O4a`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    
    #10
    if( (length(grep("O6.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`O6-O3-O4-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #11
    if( (length(grep("O6.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
      rasterImage(piclist$`O6-O3-O4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }  
  }
  
  #If non-hangingbird
  if(hangingbird==0){  
    #8
    if( (length(grep("O5.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`O5-O3-O4-MO1or2b`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #9
    if( (length(grep("O5.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
      rasterImage(piclist$`O5-O3-O4b`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    
    #10.a
    if( (length(grep("O6.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`O6-O3-O4-MO1or2a`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #11.a
    if( (length(grep("O6.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
      rasterImage(piclist$`O6-O3-O4a`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }  
  }

  #12
  if( (length(grep("O5.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1 && (length(grep("O3.O4.O5", names(proportiontable)[ub])))==0){
    rasterImage(piclist$`O5-O4-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #13
  if( (length(grep("O5.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0 && (length(grep("O3.O4.O5", names(proportiontable)[ub])))==0){
    rasterImage(piclist$`O5-O4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }  
  #14
  if( (length(grep("O6.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1 && (length(grep("O3.O4.O6", names(proportiontable)[ub])))==0){
    rasterImage(piclist$`O6-O4-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #15
  if( (length(grep("O6.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0 && (length(grep("O3.O4.O6", names(proportiontable)[ub])))==0){
    rasterImage(piclist$`O6-O4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  
  #18
  if( (length(grep("O5", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1 && (length(grep("O4", names(proportiontable)[ub])))==0){
    rasterImage(piclist$`O5-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #19
  if( (length(grep("O5", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0 && (length(grep("O4", names(proportiontable)[ub])))==0){
    rasterImage(piclist$`O5`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #20
  if( (length(grep("O6", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1 && (length(grep("O4", names(proportiontable)[ub])))==0){
    rasterImage(piclist$`O6-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #21
  if( (length(grep("O6", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0 && (length(grep("O4", names(proportiontable)[ub])))==0){
    rasterImage(piclist$`O6`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  ###############################################  
  #22
  if( (length(grep("OPABB", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`OPABB`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }  
  #23
  if( (length(grep("OPABF", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`OPABF`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  } 
  #24
  if( (length(grep("OPABH", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`OPABH`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #25
  if( (length(grep("OPABT", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`OPABT`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #26
  if( (length(grep("OPABW", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`OPABW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #27
  if( (length(grep("OPAC1", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`OPAC1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #28
  if( (length(grep("OPAC2", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`OPAC2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #29
  if( (length(grep("OPAC3", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`OPAC3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #30
  if( (length(grep("OPAH", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`OPAH`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #31
  if( (length(grep("OPAMTT", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`OPAMTT`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #32
  if( (length(grep("OPAMW", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`OPAMW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #33
  if( (length(grep("OPATW", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`OPATW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #34
  if( (length(grep("OPAW", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`OPAW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #35
  if( (length(grep("OPMB", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`OPMB`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #36
  if( (length(grep("OPMF", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`OPMF`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #37
  if( (length(grep("OPMH", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`OPMH`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #38
  if( (length(grep("OPMT", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`OPMT`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #39
  if( (length(grep("OPMW", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`OPMW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  ######################################
  #40
  if( (length(grep("PC1", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`PC1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #41
  if( (length(grep("PC2", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`PC2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #42
  if( (length(grep("PC3", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`PC3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #43
  if( (length(grep("PC4", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`PC4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #44
  if( (length(grep("PU1", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`PU1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  ######################################  
  #45
  if( (length(grep("RM1", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`RM1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #46
  if( (length(grep("RM2", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`RM2`, xleft = 1, ybottom=1, xright=100,ytop=90)    
  }
  #47
  if( (length(grep("RM3", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`RM3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #48
  if( (length(grep("RM4", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`RM4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  ######################################  
  #49
  if( (length(grep("SP1", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`SP1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }
  #50
  if( (length(grep("SP2", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`SP2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }  
  #51
  if( (length(grep("SP3", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`SP3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  } 
  #52
  if( (length(grep("SP4", names(proportiontable)[ub])))==1){
    rasterImage(piclist$`SP4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
  }  
  
  ############  
  text(50,88,ub,cex=1.5) #ADD BEHAVIOR ID NUMBER
  ########################################### 
  ############  
 # box(which="plot",lty="solid",lwd=30,col=SELFCOLORS.colors[ub]) #ADD HEATMAPPED COLORED BOX AROUND PLOT, CORRESPONDING TO LOOPS
  arrow.head.color<-ifelse(is.na(SELFCOLORS.colors[ub]),NA,SELFCOLORS.colors[ub])
  
#   diagram::curvedarrow(from=c(45,80), to=c(55,80), lwd = 7, lty = 1, 
#                        # lcol = "black", arr.col = "black", 
#                        lcol = SELFCOLORS.colors[ub], arr.col = NA,
#                        arr.pos = 1, curve = -2, dr = 0.01, arr.type="triangle",endhead=TRUE,
#                        segment = c(0, 1))
#   
 # polygon(x=c(58,52,55),y=c(86,86,75),col=arrow.head.color,border=arrow.head.color)
  if(!(is.na(arrow.head.color))) {
  polygon(x=30+0.4*selflooppoints[,1],y=79+0.25*selflooppoints[,2],col=arrow.head.color)
  }
  ########################################### 
  setwd("C:/Users/Rusty/Amazon Drive/BOP/Video_analysis/Pictograms/temp")
  pngplot <- recordPlot() 
  xyz<-names(proportiontable)[ub]
  individualbehaviorname<-paste(xyz,"png",sep=".")

  png(file = individualbehaviorname, bg = "transparent")
  replayPlot(pngplot)
  dev.off()
  
  completepicture<-readPNG(individualbehaviorname)
  rimg <- as.raster(completepicture) # raster multilayer object
  
  newpgname<-gsub('.{4}$', '',xyz)
  behaviorplotlist[[ub]]<-rimg
  names(behaviorplotlist)[[ub]]<-xyz
  
  setwd(dir)
}

#dev.off()



#  behaviorplotlist is a named list corresponding to all unique behaviors from a given trial



setwd(origwd)









# NETWORK SUMMARY VARIABLES

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#6.1 Density

#The proportion of present edges from all possible edges in the network.

#densityscore<-ecount(showme.looped)/(vcount(showme.looped)*(vcount(showme.looped)-1)) #for a directed network
densityscore<-ecount(showme.looped)/(vcount(showme.looped)^2)#for a directed network WITH self-loops

#Mean degree
deg <- degree(showme.looped,loops=FALSE)
mean.degree<-mean(deg)

#Mean path length
mean_path_length<-mean_distance(showme.looped, directed=T)


#network diameter is the longest geodesic distance (length of the shortest path between two nodes) in the network.
network.diameter<-diameter(showme.looped, directed=T,weights=NA)


#average clustering coefficient
avg.cluster.coeff<-transitivity(showme.looped)

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

printme<-data.frame(FileisNamed,belongsto,TotalTimeWatched,behaviorsscored,n.unique.behavs,
                    uniquebehaviors,densityscore,mean.degree,mean_path_length,network.diameter,avg.cluster.coeff,
                    iSimpson,relativebehavioraldiversity,iSimpsonT,relativetransitiondiversity,
                    PropSelfTransitioning,Smallworldness,
                    Nmales,Nfemales,Nfemales.disp,Time.upright,Time.inverted,Time.orientedtofemale,
                    femaleorientedbouts,stateANDduration,odd.behaviors,off.ss.eq,happenedwhileoff)






#CONNECTIVITY measures reflect degree to which network differs from complete network
#Edge density: % of edges compared to maximum
#Average degree: Avg. number of links
#Average path length: avg of shortest pathes between reachable nodes
#Network diameter: longest of shortest paths

#CENTRALITY measures quantify heterogeneity in network structure
#Average clustering coefficient:
#Components











#V(showme.looped)$size<-80
netname<-paste(name2,"Network.pdf",sep='_')

setwd(location)
pdf(netname,width= 12, height= 12,family="NimbusRom")
  op <- par(mfrow=c(1,1))# rows, columns
  par(mar=c(0,0,0,0)) 
  # plot(clp, showme.looped,edge.curved=.1,edge.arrow.size = 0.45,vertex.size=18,layout=l,rescale=FALSE,
  #      xlim=range(l[,1]),ylim=range(l[,2]))
  
  
  # plot(clp, showme.looped,edge.curved=.1,edge.arrow.size = 0.45,vertex.size=18,layout=l,rescale=TRUE,
  #      xlim=c(-1,1),ylim=c(-1,1))
  
  V(showme.looped)$raster <- behaviorplotlist
  V(showme)$raster <- behaviorplotlist
  

  
  
# print(  plot(clp, showme.looped,edge.curved=.1,layout=l,rescale=TRUE,vertex.shape="raster",
#        vertex.label=NA, 
#        #vertex.size=15*(V(showme.looped)$selfsize), # resizes vertices based on self-self transitions
#        vertex.size=25,
#        edge.color=col.tr,
#        edge.arrow.size=0.35,
#        
#        edge.width=3,
#        edge.arrow.width=1.6,
#        edge.loop.angle=1.5,
#        
#        xlim=c(-1,1),ylim=c(-1,1)))
  
print(  plot(clp, showme,edge.curved=.1,layout=l,rescale=TRUE,vertex.shape="raster",
             vertex.label=NA, 
             #vertex.size=15*(V(showme.looped)$selfsize), # resizes vertices based on self-self transitions
             vertex.size=25,
             edge.color=col.tr,
             edge.arrow.size=0.35,
             
             edge.width=3,
             edge.arrow.width=1.6,
             edge.loop.angle=1.5,
             
             xlim=c(-1,1),ylim=c(-1,1)))

print(text(0.9,-0.75,paste("Timescored(s) = ",time.value,sep=""),col="black",cex=0.75,adj=(0)))
print(text(0.9,-0.8,paste("behaviors = ",uniquebehaviors,sep=""),col="black",cex=0.75,adj=(0)))  
print(text(0.9,-0.85,paste("mean deg = ",round(mean.degree,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
print(text(0.9,-0.9,paste("density score = ",round(densityscore,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
print(text(0.9,-0.95,paste("mean path L = ",round(mean_path_length,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
print(text(0.9,-1,paste("network.diameter = ",network.diameter,sep=""),col="black",cex=0.75,adj=(0)))  
print(text(0.9,-1.05,paste("Smallworldness = ",Smallworldness,sep=""),col="black",cex=0.75,adj=(0)))


  #tkplot(clp, showme.looped)
#####  par("usr")
  
  # library(plotrix)
  # 
  # scalelist<-list()
  # gettinsmaller<-seq(from=1,to=0.52, by=-0.02)
  # for (vc in 1:25){
  # 
  #   
  #   gettingsmaller<-gettinsmaller[vc]
  #   pointlocations<-cbind(rescale(l[,1],gettingsmaller*c(-1,1)),rescale(l[,2],gettingsmaller*c(-1,1)))
  #   par(new=TRUE)
  #   plot(pointlocations,xlab="", ylab="",xlim=c(-1,1),ylim=c(-1,1),pch=vc)
  # 
  # }
  # dev.off()
  
  
  
  # for(bbb in 1:uniquebehaviors){
  #   rasterImage(behaviorplotlist[[bbb]], xleft = pointlocations[bbb,1], xright=pointlocations[bbb,1]+0.15,
  #               ybottom=pointlocations[bbb,2], ytop=pointlocations[bbb,2]+0.15)  
  #   par(new=TRUE)
  #   plot(pointlocations,xlab="", ylab="",xlim=c(-1,1),ylim=c(-1,1))
  # }
  # 
  # pointlocations<-cbind(rescale(l[,1],0.92*c(-1,1)),rescale(l[,2],0.92*c(-1,1)))
  # pointlocations<-pointlocations
  # 
  # for(bbb in 1:uniquebehaviors){
  #   rasterImage(behaviorplotlist[[bbb]], xleft = pointlocations[bbb,1], xright=pointlocations[bbb,1]+0.15,
  #               ybottom=pointlocations[bbb,2], ytop=pointlocations[bbb,2]+0.15)  
  #   par(new=TRUE)
  #   
  # }

 dev.off()

#Another network program
#geph

}

#FILTERED BEHAVIORS (minus behaviors that only happened once)
no.OFFs.filtered1<-no.OFFs[which(duplicated(no.OFFs)),,drop=FALSE]
  
  summarytable.filtered1<-table(no.OFFs.filtered1)
  behaviorcounts.filtered1<-table(no.OFFs.filtered1)
  proportiontable.filtered1<-prop.table(summarytable.filtered1)
  uniquebehaviors.filtered1<-length(summarytable.filtered1)  

seven<-no.OFFs.filtered1
  
if(uniquebehaviors.filtered1>1){  
  ##This loop creates the denominator for the calculation of the inverse Simpson Index for behavioral diversity
  val2.filtered1<-0
  for (u in 1:uniquebehaviors.filtered1){
    proptosquare<-proportiontable[u]
    val <-proptosquare^2.0
    val2.filtered1<-val+val2.filtered1
  }			
  #inverse Simpson Index for behavioral diversity
  iSimpson.filtered1<-1.0/as.numeric(val2.filtered1)
  relativebehavioraldiversity.filtered1<-iSimpson.filtered1/(uniquebehaviors.filtered1) #relative color diversity, accounting for # of color classes
  
  ##########################################################################################################
  #uses "createSequenceMatrix" function (from *markovchain*) to calculate transition matrix
  c.no.OFFs<-as.character(no.OFFs.filtered1)
  
  
  TransitionMatrix.markov<-createSequenceMatrix(c.no.OFFs,toRowProbs = TRUE,sanitize=FALSE)
  TransitionMatrix.cum.filtered1<-createSequenceMatrix(c.no.OFFs,sanitize=FALSE)
  
  #TransitionMatrix<-TransitionMatrix[-ncol(TransitionMatrix),-ncol(TransitionMatrix)]
  
  time.value<-sum(TransitionMatrix.cum.filtered1)/10 #seconds watched
  
  
  
  noselfs<-TransitionMatrix.cum.filtered1
  diag(noselfs)<-NA
  
  prop.table.excludeNAs<- function(x) {x/sum(x, na.rm=TRUE)}
  proportiontable.trans.Simpson<-prop.table.excludeNAs(noselfs) #get rid of diagonal/self-transitions for Simpson Analyses
  
  proportiontable.trans<-prop.table.excludeNAs(TransitionMatrix.cum.filtered1) 
  
  
  if(uniquebehaviors.filtered1>1){
    val3<-matrix(nrow=(((uniquebehaviors.filtered1^2)-uniquebehaviors.filtered1)/2),ncol=2)
    rownum<-1
    for (u in 1:uniquebehaviors.filtered1){
      for(v in u:uniquebehaviors.filtered1){#iteratively reduces columns analyzed in next loop to avoid double counting transitions 
        if (u!=v){
          val3[rownum,2]<-as.numeric(proportiontable.trans.Simpson[u,v]+proportiontable.trans.Simpson[v,u])
          val3[rownum,1]<-paste(row.names(proportiontable.trans.Simpson)[u],colnames(proportiontable.trans.Simpson)[v],sep="-")
          rownum<-rownum+1
        }
      }
    }
    val4<-as.data.frame(val3)
    val4[,2]<-as.numeric(as.character(val4[,2])) 
    
    simpsonT<-0
    
    for(z in 1:nrow(val4)){
      propsq<-val4[z,2]^2.0
      simpsonT<-simpsonT+propsq
    }
    
    
    iSimpsonT.filtered1<-1.0/simpsonT #Calculates inverse Simpson index for horizontal transitions
    relativetransitiondiversity.filtered1<-iSimpsonT.filtered1/(uniquebehaviors.filtered1*(uniquebehaviors.filtered1-1)/2)  
  } else {
    
    iSimpsonT.filtered1 <- NA
    relativetransitiondiversity.filtered1<- NA
  }
  
  
  PropSelfTransitioning.filtered1<-sum(diag(TransitionMatrix.cum.filtered1))/sum(TransitionMatrix.cum.filtered1)
  
  network.TransitionMatrix.cum.filtered1<-network(proportiontable.trans,directed=TRUE,loops=TRUE,hyper=FALSE,matrix.type="adjacency")
  
  
  
  showme<-graph_from_adjacency_matrix(proportiontable.trans,mode="directed",weighted=TRUE,diag=FALSE,add.colnames=TRUE)
  showme.undirected<-graph_from_adjacency_matrix(proportiontable.trans,mode="undirected",weighted=TRUE,diag=FALSE,add.colnames=TRUE)
  
  showme.looped<- graph_from_adjacency_matrix(proportiontable.trans,mode="directed",weighted=TRUE,diag=TRUE,add.colnames=TRUE)
  showme.looped.undirected<- graph_from_adjacency_matrix(proportiontable.trans,mode="undirected",weighted=TRUE,diag=TRUE,add.colnames=TRUE)
  
  
  
  E(showme)$width <- 1+E(showme)$weight/50
  V(showme)$nombres<-colnames(proportiontable.trans)
  
  E(showme.looped)$width <- 1+E(showme.looped)$weight/50
  V(showme.looped)$nombres<-colnames(proportiontable.trans)
  
  # Compute approx self-transition and use that to set node size:
  selfsizesscaled<-(diag(proportiontable.trans))*100+1
  V(showme)$selfsize<-sqrt(selfsizesscaled)
  V(showme)$size <- ((selfsizesscaled)^2)*9
  
  V(showme.looped)$selfsize<-log(selfsizesscaled)+1
  
  
  
  #####################################################
  # Network Clustering
  # l <- layout.fruchterman.reingold(showme.looped)
  # l <- layout_with_kk(showme.looped)
  # l <- layout_with_lgl(showme.looped)
  l <- layout_with_graphopt(showme.looped)
  l2 <-layout_with_graphopt(showme)
  #l <- layout_nicely(showme.looped)
  #behaviorplotlist
  
  
  #clp <- cluster_label_prop(showme.looped)
  clp <-cluster_label_prop(showme.looped.undirected) # input graph should be undirected to make sense.
  clp2<-cluster_label_prop(showme.undirected)
  # setwd("C:/Users/Rusty/Amazon Drive/BOP/Video_analysis")
  # source("PictoGrammerMega_Dec13ab.R", chdir = F)
  # #PictoGrammer<-function(proportiontable,hangingbird)
  #
  
  
  
  
  ################################################################################################################################
  #PictoGrammer(proportiontable.trans,hangingbird,uniquebehaviors)
  
  origwd<-dir
  setwd(dir<-"C:/Users/Rusty/Amazon Drive/BOP/Video_analysis/Pictograms")
  
  pictogram.names<-list.files(dir,full.names=FALSE,pattern=".png")
  
  proportiontable<-proportiontable
  hangingbird<-hangingbird
  
  
  selflooppoints<-read.table("SelfLoop.csv",sep=",",header=F)
  
  #Creates empty list, then loops through pictogram directory, pulling in each PG, adding it to the list, and naming it
  piclist<-list()
  for (pg in 1:56){
    img <- readPNG(pictogram.names[pg])  
    
    if(pg>23 && pg <42) {
      w <- matrix(rgb(img[,,1],img[,,2],img[,,3], img[,,4] * 0.75), nrow=dim(img)[1]) #0.5 is alpha
      rimg <- as.raster(w) # raster multilayer object
    } else {
      rimg <- as.raster(img) # raster multilayer object  
    }
    pgname<-gsub('.{4}$', '',pictogram.names[pg])
    piclist[[pg]]<-rimg
    names(piclist)[[pg]]<-pgname
  }
  
  #Loops through each unique behavior (from proportiontable) and creates a new, unique, behavioral pictogram
  behaviorplotlist<-list()
  par(bg="transparent")
  plot.new()
  goo=1
  #pdf("Behav30-trans-fix.pdf",width= 12, height= 12,family="NimbusRom")
  #op <- par(mfrow=c(6,6))# rows, columns
  #par(mar=c(1,1,1,1)) 
  
  
  colorweights<-E(showme.looped)$weight/min(E(showme.looped)$weight) 
  colorweights<-log(colorweights+0.01)+1
  
  Lab.palette<-colorRampPalette(c("blue",  "red"),space = "Lab")
  #  Lab.palette<-colorRampPalette(c("#ffeda0","#feb24c","#f03b20"),space = "Lab")
  
  
  E(showme.looped)$color<-Lab.palette(max(colorweights))[colorweights]
  
  col.tr <- grDevices::adjustcolor(E(showme.looped)$color, alpha=1) #adds translucency to arrows
  
  SELFCOLORS<-log(diag(proportiontable.trans)/min(E(showme.looped)$weight))+1 
  SELFCOLORS.colors<-SELFCOLORS
  #Little loop that adds same heatmap scaled colors to self-transitions, for use in behavior boxes
  for(ggg in 1:length(SELFCOLORS)){
    selfvalue<-SELFCOLORS[ggg]
    if(selfvalue>0){
      SELFCOLORS.colors[ggg]<-Lab.palette(max(colorweights))[selfvalue]
    } else {
      SELFCOLORS.colors[ggg]<-NA
    }
    
  }
  
  
  for (ub in 1:uniquebehaviors.filtered1){
    par(bg="transparent")
    graphics::plot(goo, type ="l", xlab="", ylab="", xlim=c(0, 100), ylim=c(0, 100),axes=FALSE)
    
    #53
    if( (length(grep("SS1", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`SS1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #54
    if( (length(grep("SS2", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`SS2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    ##########################################################  
    #2
    if( (length(grep("BP1", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`BP1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #1
    if( (length(grep("BP1.BP2", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`BP1-BP2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #3
    if( (length(grep("BP3", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`BP3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #4
    if( (length(grep("O1", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`O1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #5
    if( (length(grep("O2", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
      rasterImage(piclist$`O2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    ###############################################
    
    #If hangingbird  
    if(hangingbird==1){  
      #6
      if( (length(grep("O5.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`O5-O3-O4-MO1or2a`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #7
      if( (length(grep("O5.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O5-O3-O4a`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      
      #10
      if( (length(grep("O6.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`O6-O3-O4-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #11
      if( (length(grep("O6.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O6-O3-O4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }  
    }
    
    #If non-hangingbird
    if(hangingbird==0){  
      #8
      if( (length(grep("O5.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`O5-O3-O4-MO1or2b`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #9
      if( (length(grep("O5.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O5-O3-O4b`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      
      #10.a
      if( (length(grep("O6.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`O6-O3-O4-MO1or2a`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #11.a
      if( (length(grep("O6.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O6-O3-O4a`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }  
    }
    
    #12
    if( (length(grep("O5.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1 && (length(grep("O3.O4.O5", names(proportiontable)[ub])))==0){
      rasterImage(piclist$`O5-O4-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #13
    if( (length(grep("O5.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0 && (length(grep("O3.O4.O5", names(proportiontable)[ub])))==0){
      rasterImage(piclist$`O5-O4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }  
    #14
    if( (length(grep("O6.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1 && (length(grep("O3.O4.O6", names(proportiontable)[ub])))==0){
      rasterImage(piclist$`O6-O4-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #15
    if( (length(grep("O6.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0 && (length(grep("O3.O4.O6", names(proportiontable)[ub])))==0){
      rasterImage(piclist$`O6-O4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    
    #18
    if( (length(grep("O5", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1 && (length(grep("O4", names(proportiontable)[ub])))==0){
      rasterImage(piclist$`O5-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #19
    if( (length(grep("O5", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0 && (length(grep("O4", names(proportiontable)[ub])))==0){
      rasterImage(piclist$`O5`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #20
    if( (length(grep("O6", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1 && (length(grep("O4", names(proportiontable)[ub])))==0){
      rasterImage(piclist$`O6-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #21
    if( (length(grep("O6", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0 && (length(grep("O4", names(proportiontable)[ub])))==0){
      rasterImage(piclist$`O6`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    ###############################################  
    #22
    if( (length(grep("OPABB", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPABB`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }  
    #23
    if( (length(grep("OPABF", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPABF`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    } 
    #24
    if( (length(grep("OPABH", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPABH`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #25
    if( (length(grep("OPABT", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPABT`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #26
    if( (length(grep("OPABW", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPABW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #27
    if( (length(grep("OPAC1", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPAC1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #28
    if( (length(grep("OPAC2", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPAC2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #29
    if( (length(grep("OPAC3", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPAC3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #30
    if( (length(grep("OPAH", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPAH`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #31
    if( (length(grep("OPAMTT", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPAMTT`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #32
    if( (length(grep("OPAMW", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPAMW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #33
    if( (length(grep("OPATW", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPATW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #34
    if( (length(grep("OPAW", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPAW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #35
    if( (length(grep("OPMB", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPMB`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #36
    if( (length(grep("OPMF", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPMF`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #37
    if( (length(grep("OPMH", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPMH`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #38
    if( (length(grep("OPMT", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPMT`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #39
    if( (length(grep("OPMW", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPMW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    ######################################
    #40
    if( (length(grep("PC1", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`PC1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #41
    if( (length(grep("PC2", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`PC2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #42
    if( (length(grep("PC3", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`PC3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #43
    if( (length(grep("PC4", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`PC4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #44
    if( (length(grep("PU1", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`PU1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    ######################################  
    #45
    if( (length(grep("RM1", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`RM1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #46
    if( (length(grep("RM2", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`RM2`, xleft = 1, ybottom=1, xright=100,ytop=90)    
    }
    #47
    if( (length(grep("RM3", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`RM3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #48
    if( (length(grep("RM4", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`RM4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    ######################################  
    #49
    if( (length(grep("SP1", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`SP1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }
    #50
    if( (length(grep("SP2", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`SP2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }  
    #51
    if( (length(grep("SP3", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`SP3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    } 
    #52
    if( (length(grep("SP4", names(proportiontable)[ub])))==1){
      rasterImage(piclist$`SP4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
    }  
    
    ############  
    text(50,88,ub,cex=1.5) #ADD BEHAVIOR ID NUMBER
    ########################################### 
    ############  
    # box(which="plot",lty="solid",lwd=30,col=SELFCOLORS.colors[ub]) #ADD HEATMAPPED COLORED BOX AROUND PLOT, CORRESPONDING TO LOOPS
    arrow.head.color<-ifelse(is.na(SELFCOLORS.colors[ub]),NA,SELFCOLORS.colors[ub])
    
    #   diagram::curvedarrow(from=c(45,80), to=c(55,80), lwd = 7, lty = 1, 
    #                        # lcol = "black", arr.col = "black", 
    #                        lcol = SELFCOLORS.colors[ub], arr.col = NA,
    #                        arr.pos = 1, curve = -2, dr = 0.01, arr.type="triangle",endhead=TRUE,
    #                        segment = c(0, 1))
    #   
    # polygon(x=c(58,52,55),y=c(86,86,75),col=arrow.head.color,border=arrow.head.color)
    if(!(is.na(arrow.head.color))) {
      polygon(x=30+0.4*selflooppoints[,1],y=79+0.25*selflooppoints[,2],col=arrow.head.color)
    }
    ########################################### 
    setwd("C:/Users/Rusty/Amazon Drive/BOP/Video_analysis/Pictograms/temp")
    pngplot <- recordPlot() 
    xyz<-names(proportiontable)[ub]
    individualbehaviorname<-paste(xyz,"png",sep=".")
    
    png(file = individualbehaviorname, bg = "transparent")
    replayPlot(pngplot)
    dev.off()
    
    completepicture<-readPNG(individualbehaviorname)
    rimg <- as.raster(completepicture) # raster multilayer object
    
    newpgname<-gsub('.{4}$', '',xyz)
    behaviorplotlist[[ub]]<-rimg
    names(behaviorplotlist)[[ub]]<-xyz
    
    setwd(dir)
  }
  
  #dev.off()
  
  
  
  #  behaviorplotlist is a named list corresponding to all unique behaviors from a given trial
  
  
  
  setwd(origwd)
  
  
  
  
  
  
  
  
  
  # NETWORK SUMMARY VARIABLES
  
  #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
  #6.1 Density
  
  #The proportion of present edges from all possible edges in the network.
  
  densityscore.filtered1<-ecount(showme.looped)/(vcount(showme.looped)*(vcount(showme.looped)-1)) #for a directed network
  
  #Mean degree
  deg <- degree(showme.looped,loops=FALSE)
  mean.degree.filtered1<-mean(deg)
  
  #Mean path length
  mean_path_length.filtered1<-mean_distance(showme.looped, directed=T)
  
  
  #network diameter is the longest geodesic distance (length of the shortest path between two nodes) in the network.
  network.diameter.filtered1<-diameter(showme.looped, directed=T,weights=NA)
  
  
  #average clustering coefficient
  avg.cluster.coeff.filtered1<-transitivity(showme.looped)
  
  ############################################
  #SMALL WORLDNESS --- video-wise, filtered
  #number of nodes/vertices in graph
  vertices<- uniquebehaviors.filtered1
  #number of edges in G(n,m) graph
  edges<- sum(TransitionMatrix.cum.filtered1!=0)
  
  rando.network.filtered1<-sample_gnm(n=vertices, m=edges, directed = TRUE, loops = TRUE)
  
  Trobserved<-avg.cluster.coeff.filtered1
  mean.Trrandom<-transitivity(rando.network.filtered1)
  
  SPobserved<-mean_path_length.filtered1
  mean.SPrandom<-mean_distance(rando.network.filtered1, directed=T)  
  
  Smallworldness.filtered1<- (Trobserved/mean.Trrandom)/(SPobserved/mean.SPrandom)
  ############################################
  
  
  printme<-data.frame(FileisNamed,belongsto,TotalTimeWatched,behaviorsscored,n.unique.behavs,
                      uniquebehaviors,densityscore,mean.degree,mean_path_length,network.diameter,avg.cluster.coeff,
                      iSimpson,relativebehavioraldiversity,iSimpsonT,relativetransitiondiversity,
                      PropSelfTransitioning,Smallworldness,
                      
                      uniquebehaviors.filtered1,densityscore.filtered1,mean.degree.filtered1,mean_path_length.filtered1,network.diameter.filtered1,avg.cluster.coeff.filtered1,
                      iSimpson.filtered1,relativebehavioraldiversity.filtered1,iSimpsonT.filtered1,relativetransitiondiversity.filtered1,
                      PropSelfTransitioning.filtered1,Smallworldness.filtered1,
                      Nmales,Nfemales,Nfemales.disp,Time.upright,Time.inverted,Time.orientedtofemale,
                      femaleorientedbouts,stateANDduration,odd.behaviors,off.ss.eq,happenedwhileoff)
  
  
  
  
  
  
  #CONNECTIVITY measures reflect degree to which network differs from complete network
  #Edge density: % of edges compared to maximum
  #Average degree: Avg. number of links
  #Average path length: avg of shortest pathes between reachable nodes
  #Network diameter: longest of shortest paths
  
  #CENTRALITY measures quantify heterogeneity in network structure
  #Average clustering coefficient:
  #Components
  
  
  
  
  
  
  
  
  
  
  
  #V(showme.looped)$size<-80
  netname<-paste(name2,"Network_filtered1.pdf",sep='_')
  
  setwd(location)
  pdf(netname,width= 12, height= 12,family="NimbusRom")
  op <- par(mfrow=c(1,1))# rows, columns
  par(mar=c(0,0,0,0)) 
  # plot(clp, showme.looped,edge.curved=.1,edge.arrow.size = 0.45,vertex.size=18,layout=l,rescale=FALSE,
  #      xlim=range(l[,1]),ylim=range(l[,2]))
  
  
  # plot(clp, showme.looped,edge.curved=.1,edge.arrow.size = 0.45,vertex.size=18,layout=l,rescale=TRUE,
  #      xlim=c(-1,1),ylim=c(-1,1))
  
  V(showme.looped)$raster <- behaviorplotlist
  V(showme)$raster <- behaviorplotlist
  
  
  
  
  # print(  plot(clp, showme.looped,edge.curved=.1,layout=l,rescale=TRUE,vertex.shape="raster",
  #        vertex.label=NA, 
  #        #vertex.size=15*(V(showme.looped)$selfsize), # resizes vertices based on self-self transitions
  #        vertex.size=25,
  #        edge.color=col.tr,
  #        edge.arrow.size=0.35,
  #        
  #        edge.width=3,
  #        edge.arrow.width=1.6,
  #        edge.loop.angle=1.5,
  #        
  #        xlim=c(-1,1),ylim=c(-1,1)))
  
  print(  plot(clp, showme,edge.curved=.1,layout=l,rescale=TRUE,vertex.shape="raster",
               vertex.label=NA, 
               #vertex.size=15*(V(showme.looped)$selfsize), # resizes vertices based on self-self transitions
               vertex.size=25,
               edge.color=col.tr,
               edge.arrow.size=0.35,
               
               edge.width=3,
               edge.arrow.width=1.6,
               edge.loop.angle=1.5,
               
               xlim=c(-1,1),ylim=c(-1,1)))
  
  print(text(0.9,-0.75,paste("Timescored(s) = ",time.value,sep=""),col="black",cex=0.75,adj=(0)))
  print(text(0.9,-0.8,paste("behaviors = ",uniquebehaviors.filtered1,sep=""),col="black",cex=0.75,adj=(0)))  
  print(text(0.9,-0.85,paste("mean deg = ",round(mean.degree.filtered1,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
  print(text(0.9,-0.9,paste("density score = ",round(densityscore.filtered1,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
  print(text(0.9,-0.95,paste("mean path L = ",round(mean_path_length.filtered1,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
  print(text(0.9,-1,paste("network.diameter = ",network.diameter.filtered1,sep=""),col="black",cex=0.75,adj=(0)))  
  print(text(0.9,-1.05,paste("Smallworldness = ",Smallworldness.filtered1,sep=""),col="black",cex=0.75,adj=(0)))
  
  #tkplot(clp, showme.looped)
  #####  par("usr")
  
  # library(plotrix)
  # 
  # scalelist<-list()
  # gettinsmaller<-seq(from=1,to=0.52, by=-0.02)
  # for (vc in 1:25){
  # 
  #   
  #   gettingsmaller<-gettinsmaller[vc]
  #   pointlocations<-cbind(rescale(l[,1],gettingsmaller*c(-1,1)),rescale(l[,2],gettingsmaller*c(-1,1)))
  #   par(new=TRUE)
  #   plot(pointlocations,xlab="", ylab="",xlim=c(-1,1),ylim=c(-1,1),pch=vc)
  # 
  # }
  # dev.off()
  
  
  
  # for(bbb in 1:uniquebehaviors){
  #   rasterImage(behaviorplotlist[[bbb]], xleft = pointlocations[bbb,1], xright=pointlocations[bbb,1]+0.15,
  #               ybottom=pointlocations[bbb,2], ytop=pointlocations[bbb,2]+0.15)  
  #   par(new=TRUE)
  #   plot(pointlocations,xlab="", ylab="",xlim=c(-1,1),ylim=c(-1,1))
  # }
  # 
  # pointlocations<-cbind(rescale(l[,1],0.92*c(-1,1)),rescale(l[,2],0.92*c(-1,1)))
  # pointlocations<-pointlocations
  # 
  # for(bbb in 1:uniquebehaviors){
  #   rasterImage(behaviorplotlist[[bbb]], xleft = pointlocations[bbb,1], xright=pointlocations[bbb,1]+0.15,
  #               ybottom=pointlocations[bbb,2], ytop=pointlocations[bbb,2]+0.15)  
  #   par(new=TRUE)
  #   
  # }
  
  dev.off()
  
  #Another network program
  #geph
} else {
  c.no.OFFs<-as.character(no.OFFs.filtered1)
  TransitionMatrix.cum.filtered1<-  createSequenceMatrix(c.no.OFFs,sanitize=FALSE)
  
  densityscore.filtered1<-mean.degree.filtered1<-mean_path_length.filtered1<-network.diameter.filtered1<-avg.cluster.coeff.filtered1<-iSimpson.filtered1<-relativebehavioraldiversity.filtered1<-iSimpsonT.filtered1<-relativetransitiondiversity.filtered1<-PropSelfTransitioning.filtered1<-Smallworldness.filtered1<-NA
  
  printme<-data.frame(FileisNamed,belongsto,TotalTimeWatched,behaviorsscored,n.unique.behavs,
                      uniquebehaviors,densityscore,mean.degree,mean_path_length,network.diameter,avg.cluster.coeff,
                      iSimpson,relativebehavioraldiversity,iSimpsonT,relativetransitiondiversity,
                      PropSelfTransitioning,Smallworldness,
                      
                      uniquebehaviors.filtered1,densityscore.filtered1,mean.degree.filtered1,mean_path_length.filtered1,network.diameter.filtered1,avg.cluster.coeff.filtered1,
                      iSimpson.filtered1,relativebehavioraldiversity.filtered1,iSimpsonT.filtered1,relativetransitiondiversity.filtered1,
                      PropSelfTransitioning.filtered1,Smallworldness.filtered1,
                      Nmales,Nfemales,Nfemales.disp,Time.upright,Time.inverted,Time.orientedtofemale,
                      femaleorientedbouts,stateANDduration,odd.behaviors,off.ss.eq,happenedwhileoff)
}





  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  
###PLOTTING PORTION OF SCRIPT
#####################################################################################################################################
  EventBehaviors<-c("BP2","O1","O2","OPMH","OPMB","OPMF","OPMW","OPMT","OPAH","OPAC1","OPAC2","OPAC3","OPAW","OPATW","OPAMW","OPAMTT",
                    "OPABH","OPABB","OPABF","OPABW","OPABT","PC1","PC2","PC3","PC4","MO1","SP1","SP2","SP3","SP4")
  
  
  DurationBehaviors2<-c("BP1","BP3","SS1","SS2","O3","O4","O5","OFF","O6","O7","O8","OPMH.D","OPMB.D","OPMF.D","OPMW.D","OPMT.D","RepWing","RepTail","RepHead","RepTorso",
                        "RM1","RM2","RM3","RM4","MO2","PU1")
  SPDurations.start<-c("SP1","SP2","SP3")
  SPDurations.stop<-c("SP2","SP3","SP4")
  
  forplot<-matrix(nrow=(TotalTime*10+3),ncol=length(c(DurationBehaviors2,EventBehaviors))) #Creates empty matrix with 1 column/behavior and 10 rows/sec watched (0.1s resolution) (+2 accounts for 0 and buffer at end)
  colnames(forplot)<-c(DurationBehaviors2,EventBehaviors)
  
  fulltime<-as.character(seq(0,TotalTime+0.2,0.1))
  rownames(forplot)<-fulltime
  forplot<-data.frame(forplot)
  
  q.time<-as.character(round(table$time/0.1)*0.1)
  table<-cbind(table,q.time)
  
  #eventdata<-table[which(table$code %in% EventBehaviors),]
  
  for (jj in EventBehaviors){
    eventtable<-table[which(table$code %in% EventBehaviors),]
    #behav.obs <- eventtable[grep(m, eventtable$code), ]
    use<-subset(eventtable,code==jj)
    use<-unique(use)
    rowlist<-match(use$q.time,rownames(forplot))
    forplot[rowlist,jj]<-1 #Adds 1 to every 0.1 sec time window that matches for each behavior at each time
    
  }
  
  
  for (kk in DurationBehaviors2){
    use<-table[grep(kk, table$code), ]
    use<-unique(use)
    bouts<-nrow(use)/2
    if(bouts>0.6){#Only if there are >0 bouts, do the following (adding in 1s during the time periods when each duration behavior is active)
      for (ll in 1:bouts){
        start<-(ll-1)*2+1
        end<-start+1
        subuse<-use[start:end,]
        listotimes<-as.numeric(as.character(subuse$q.time))
        timeon<-seq(listotimes[1],listotimes[2],0.1)
        timeon2<-as.character(timeon)
        
        matchingrows<-match(timeon2,rownames(forplot)) #Can return non-matches for true matches based on floating point values
        forplot[matchingrows,kk]<-1 #Adds 1 to every 0.1 sec time window that matches for each behavior at each time
        
      }
    }
    
  }
  
  
  for (mm in 1:length(SPDurations.start)){
    use.starts<-table[grep(SPDurations.start[mm], table$code), ]
    use.ends<-table[grep(SPDurations.stop[mm], table$code), ]
  }
  
  
  DistinctColors<-c("#ae7b4a","#6264de","#48cd60","#9744bf","#50a428","#da5ac9","#91bf3c","#be77e6","#5ebd63","#a62d8a","#4dcc93","#d84498","#368939","#7753b5","#b4b43b","#587ee4",
                    "#e1aa34","#695faf","#e8792d","#4da9d9","#d04228","#4ed6c7","#db3973","#319a6b","#d3384c","#3bbac6","#e56d5c","#57b597","#e482d5","#69882b","#a756a7","#82bc75",
                    "#dd6a93","#2d7350","#b197e3","#a68222","#7399dc","#c5792b","#4167a5","#d5ae69","#7f5e9c","#aab56e","#984167","#5c894d","#e293c1","#5a6624","#a86597","#8b8239",
                    "#a74851","#e7986b","#7f5724","#dd8383","#a84e28")
  
  
  notimes<-forplot
  
  offperiod<-notimes[,c("OFF")]
  offperiod[is.na(offperiod)]<-0 #REPLACE NAs with ZERO
  
  notimes[,c("OFF")]<-NULL #removes fulltime and  OFF columns
  #remove columns with NO data
  columnstoremove<-c()
  for(rr in 1:ncol(notimes)){
    if(sum(notimes[,rr],na.rm = TRUE)==0){
      columnstoremove<-c(columnstoremove,rr)
    }
  }
  cleaned<-notimes[,-c(columnstoremove)]
  cleaned[is.na(cleaned)]<-0 #REPLACE NAs with ZERO
  
  str(cleaned)
  nnn.pred <- nrow(cleaned)
  behavs<-ncol(cleaned)
  date.seq <- seq(from=0, to=1, length=nnn.pred)
  rowtot <- apply(cleaned, 1, sum, na.rm=T )
  coltot <- apply(cleaned, 2, sum, na.rm=T )
  
  
  modifiedtime<-file.info(filename)[5]
  
  
  shadowtext <- function(x, y=NULL, labels, col='white', bg='black', 
                         theta= seq(0, 2*pi, length.out=50), r=0.1, ... ) {
    
    xy <- xy.coords(x,y)
    xo <- r*strwidth('A')
    yo <- r*strheight('A')
    
    # draw background text with small shift in x and y in background colour
    for (i in theta) {
      text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
    }
    # draw actual text in exact xy position in foreground colour
    text(xy$x, xy$y, labels, col=col, ... )
  }
  
  
  
  
  height<-behavs/50
  

 
  
  #pdf("C",width= 8, height= 12,family="NimbusRom")
  prefix<-paste("Unequal",odd.behaviors,sep="_")
  newunequalname<-paste(prefix,name,sep="_")
  newunequalname2<-paste(location,newunequalname,sep="/")
  
  newname<-ifelse(lengthodds<1,fullname,newunequalname2)
    

   
  pdf(newname,width= 8, height= 12,family="NimbusRom")
  # Cummulative Value Plot
  par(mar=c(1,0,0,0)) 
  op <- par(mfrow=c(behavs+1,1))# rows, columns (one more than behavs in 'cleaned', to include OFF periods)
  
  #Plot OFF periods at top
  yyy.old <- offperiod
  yyy <- c( yyy.old, rev(rep(0, length(yyy.old))) )
  
  behaviorname<-"OffScreen"
  
  print(plot( seq(0, 1, nnn.pred), seq(0, 1, nnn.pred),type="n",xlim = c(-.05,1.0),ylim = c(0,1.0),xlab="", ylab="",cex.main = 1.0,main = "",axes = F,col.main = "black" ))
  # Below 
  print(axis(1,at = c(0.00,1),labels=c("","")))
 
  print(axis(2, at = c(0, 0.25, 0.5, 0.75, 1.0), labels = c("0", "", "0.5","", "1")))	
  
  print(box(lwd=2))
  xxx <- c( date.seq, rev(date.seq))
  print(polygon(xxx, yyy, col="black", border=NA))
  OFFxxs<-xxx
  OFFyys<-yyy
  #text(0,0.5,behaviorname,,col="black",cex=1.75,font=4)
  #text(0,0.5,behaviorname,bg="black",col=DistinctColors[iii.nlcd],cex=1.5,font=4)
  print(shadowtext(-0.09,height,behaviorname,col="white",cex=2,adj=(0)))
  print(text(-0.09,.95,modifiedtime[1,1],col="black",cex=0.75,adj=(0)))
  print(text(-0.09,.65,name,col="black",cex=0.75,adj=(0)))

  summarytextsize<-15/behavs
  print(shadowtext(1,0.95,paste("Duration",TotalTime,sep=" "),col="red",cex=summarytextsize,adj=(1)))
  print(shadowtext(1,0.75,paste("Unequal",odd.behaviors,sep="="),col="red",cex=summarytextsize,adj=(1)))
  print(shadowtext(1,0.55,paste("SSequal",off.ss.eq,sep="="),col="red",cex=summarytextsize,adj=(1)))
  print(shadowtext(1,0.35,paste("WhileOff",happenedwhileoff,sep="="),col="red",cex=summarytextsize,adj=(1)))
  print(shadowtext(1,0.15,paste("Nmales",Nmales,sep="="),col="red",cex=summarytextsize,adj=(1)))
  
  
  print(axis(1,at = c(0.5),labels=""))
 # print(text(0.489,0.10,round(TotalTime/2,digits=1)))
  
  
  for (iii.nlcd in 1:behavs){
    yyy.old <- cleaned[,iii.nlcd]
    yyy <- c( yyy.old, rev(rep(0, length(yyy.old))) )
    
    behaviorname<-colnames(cleaned)[iii.nlcd]
    
    print(plot( seq(0, 1, nnn.pred), seq(0, 1, nnn.pred),type="n",xlim = c(-.05,1.0),ylim = c(0,1.0),xlab="", ylab="",cex.main = 1.0,main = "",axes = F,col.main = "black" ))
    # Below 
    print(axis(1,at = c(0.00,1),labels=c("","")))
    print(axis(2, at = c(0, 0.25, 0.5, 0.75, 1.0), labels = c("0", "", "0.5","", "1")))	
    
    print(box(lwd=2))
    xxx <- c( date.seq, rev(date.seq))
      
      print(polygon(xxx, yyy, col=DistinctColors[iii.nlcd], border=NA))
  
    print(polygon(OFFxxs, OFFyys, col=rgb(.2,.2,.2,alpha=.2), border=NA))
    #text(0,0.5,behaviorname,,col="black",cex=1.75,font=4)
    #text(0,0.5,behaviorname,bg="black",col=DistinctColors[iii.nlcd],cex=1.5,font=4)
    print(shadowtext(-0.09,height,behaviorname,col=DistinctColors[iii.nlcd],cex=2,adj=(0)))
    print(axis(1,at = c(0.5),labels=""))
    #print(text(0.48,0.10,round(TotalTime/2,digits=1)))
    
    
  } 
  
  dev.off()
  
  
  
######################################################################################################################################  
######################################################################################################################################
  
  
setwd(origsourcefolder)
	
sevencomponents<-list(printme,TransitionMatrix.cum,TransitionMatrix.cum.filtered1,behaviorcounts,behaviorcounts.filtered1,six,seven)
    
return(sevencomponents)
	}
ProcessBOP.cust<-function(filename,location,species.name,specialclass,hb,dimensions){
  
  #function for pulling right side of text string, for naming purposes    
  right = function (string, char){
    substr(string,nchar(string)-(char-1),nchar(string))
  }
  hangingbird<-hb
  ######Function for eliminating event behaviors during "OFF" and fixing duration behaviors
  #DURATION BEHAVIORS
  DurationBehaviors<-c("BP1","BP3","SS1","SS2","O3","O4","O5","OFF","O6","O7","O8","OPMH.D","OPMB.D","OPMF.D","OPMW.D","OPMT.D","RepWing","RepTail","RepHead","RepTorso",
                       "RM1","RM2","RM3","RM4","MO2","PU1","SP")
  DurationBehaviors.nosp<-c("BP1","BP3","SS1","SS2","O3","O4","O5","OFF","O6","O7","O8","OPMH.D","OPMB.D","OPMF.D","OPMW.D","OPMT.D","RepWing","RepTail","RepHead","RepTorso",
                            "RM1","RM2","RM3","RM4","MO2","PU1")
  DurationBehaviors.nosp.noOFF<-c("BP1","BP3","SS1","SS2","O3","O4","O5","O6","O7","O8","OPMH.D","OPMB.D","OPMF.D","OPMW.D","OPMT.D","RepWing","RepTail","RepHead","RepTorso",
                                  "RM1","RM2","RM3","RM4","MO2","PU1")
  
  EventBehaviors<-c("BP2","O1","O2","OPMH","OPMB","OPMF","OPMW","OPMT","OPAH","OPAC1","OPAC2","OPAC3","OPAW","OPATW","OPAMW","OPAMTT",
                    "OPABH","OPABB","OPABF","OPABW","OPABT","PC1","PC2","PC3","PC4","MO1","SP1","SP2","SP3","SP4")
  
  fixOFFs<-function(tabletofix){  
    offstarts<-nrow(subset(tabletofix,code=="OFF(start)"))
    offstops<-nrow(subset(tabletofix,code=="OFF(end)"))	
    tabletofix<-tabletofix[order(tabletofix$time),]
    if(offstarts>0){#If there are "OFFS", this if statement eliminates any behaviors measured during these OFFs and adds "STOPS"
      #At the beginning of the OFF period and adds "STARTS" at the end of the OFF period. Also deletes behaviors that happen while OFF.
      offstartings<-subset(tabletofix,code=="OFF(start)")[,1]
      offstoppins<-subset(tabletofix,code=="OFF(end)")[,1]
      
      for (vvv in 1:offstarts){#offstarts
        tabletofix<-tabletofix[order(tabletofix$time),]
        #stopBegin<-subset(tabletofix,code=="OFF(start)")[vvv,1]
        #stopEnd<-subset(tabletofix,code=="OFF(end)")[vvv,1]
        stopBegin<-offstartings[vvv]
        stopEnd<-offstoppins[vvv]
        
        events<- tabletofix[tabletofix$code %in% EventBehaviors,]
        gap.wise.misses<-events[which(events$time < stopEnd & events$time > stopBegin),]
        
        if(nrow(gap.wise.misses)>0){
          
          tabletofix<-tabletofix[-(which((tabletofix$time%in% gap.wise.misses$time) & (tabletofix$code%in% gap.wise.misses$code))),] #Removes rows of eventbehaviors which happen while birds is OFF screen
          #[-(which(tabletofix$time == oldend$time & tabletofix$code==oldend$code)),]
        }
        #testlist<-"O4"
        for (db in DurationBehaviors.nosp){  ##testlist #DurationBehaviors.nosp
          behav.obs <- tabletofix[grep(db, tabletofix$code), ]
          bouts<-nrow(behav.obs)/2
          if(bouts>0.6){
            starts<-behav.obs[grep("start", behav.obs$code), ]
            ends<-behav.obs[grep("end", behav.obs$code), ]
            
            oldend<-starts[1,][-1,]
            newend<-starts[1,][-1,]
            oldstart<-starts[1,][-1,]
            newstart<-starts[1,][-1,]
            
            if(nrow(starts)==nrow(ends)){ #only do this if starts/ends are equal for a given behavior
              for(nnn in 1:nrow(starts)){
                if (starts[nnn,1]==stopBegin & ends[nnn,1]<stopEnd & ends[nnn,2]!="OFF(end)"){ #S3_ If start occurs as OFF starts, and end occurs while OFF is active, delete
                  starttimetoeliminate<-starts[nnn,1]
                  stoptimetoeliminate<-ends[nnn,1]
                  startremove<-behav.obs[which(behav.obs$time == starttimetoeliminate),] 
                  endremove<-behav.obs[which(behav.obs$time == stoptimetoeliminate),] 
                  
                  gap.wise.misses2<-rbind(startremove,endremove)
                  
                  tabletofix<-tabletofix[-((match(gap.wise.misses2$time,tabletofix$time))),]
                }
                
                if(starts[nnn,1]<stopBegin & ends[nnn,1]==stopEnd & ends[nnn,2]!="OFF(end)" ){#S2_This part handles all duration behaviors that start before each (nnn) OFF starts and end as each OFF ends.
                  oldend<-ends[nnn,]
                  
                  
                  newend<-ends[1,]
                  newend$time<-stopBegin
                  
                  
                  tabletofix<-rbind(tabletofix,newend) #adds new end at start of OFF period
                  tabletofix<-tabletofix[-(which(tabletofix$time == oldend$time & tabletofix$code==oldend$code)),] #gets rid of old end behavior
                }  
                
                if(starts[nnn,1]<stopBegin & ends[nnn,1]<stopEnd & ends[nnn,1]>stopBegin & ends[nnn,2]!="OFF(end)"){#S1_This part handles all duration behaviors that start before each (nnn) OFF starts and end while OFF is on.
                  oldend<-ends[nnn,]
                  
                  
                  newend<-ends[1,]
                  newend$time<-stopBegin
                  
                  
                  tabletofix<-rbind(tabletofix,newend) #adds new end at start of OFF period
                  tabletofix<-tabletofix[-(which(tabletofix$time == oldend$time & tabletofix$code==oldend$code)),] #gets rid of old end behavior
                }  
                
                if(starts[nnn,1]>stopBegin & starts[nnn,1]<stopEnd & ends[nnn,1]==stopEnd & ends[nnn,2]!="OFF(end)" ){#E4_This part handles all duration behaviors that start while OFF and ends with ON
                  
                  oldstart<-starts[nnn,]
                  
                  newstart<-starts[1,]
                  newstart$time<-stopEnd
                  
                  tabletofix<-rbind(tabletofix,newstart) #
                  tabletofix<-tabletofix[-(which(tabletofix$time == oldstart$time & tabletofix$code==oldstart$code)),] #gets rid of old start behavior
                  
                }  
                
                
                if(starts[nnn,1]==stopBegin & ends[nnn,1]==stopEnd & ends[nnn,2]!="OFF(end)" ){#E3_This part handles all duration behaviors that start and end exactly as the OFF starts and ends
                  oldend<-ends[nnn,]
                  oldstart<-starts[nnn,]
                  
                  
                  tabletofix<-tabletofix[-(which(tabletofix$time == oldend$time & tabletofix$code==oldend$code)),] #gets rid of old end behavior
                  tabletofix<-tabletofix[-(which(tabletofix$time == oldstart$time & tabletofix$code==oldstart$code)),] #gets rid of old start behavior
                  
                }  
                
                if(starts[nnn,1]==stopBegin & ends[nnn,1]>stopEnd & ends[nnn,2]!="OFF(end)" ){#E2_This part handles all duration behaviors that start as OFF starts and end after each OFF ends.
                  oldstart<-starts[nnn,]
                  
                  
                  newstart<-starts[1,]
                  newstart$time<-stopEnd
                  
                  
                  tabletofix<-rbind(tabletofix,newstart) #adds new end at start of OFF period
                  tabletofix<-tabletofix[-(which(tabletofix$time == oldstart$time & tabletofix$code==oldstart$code)),] #gets rid of old end behavior
                } 
                
                if(starts[nnn,1]>stopBegin & starts[nnn,1]<stopEnd & ends[nnn,1]>stopEnd & ends[nnn,2]!="OFF(end)"){#E1_This part handles all duration behaviors that start while OFF is on, then after OFF is off
                  oldstart<-starts[nnn,]
                  
                  
                  newstart<-starts[1,]
                  newstart$time<-stopEnd
                  
                  
                  tabletofix<-rbind(tabletofix,newstart) #adds new end at start of OFF period
                  tabletofix<-tabletofix[-(which(tabletofix$time == oldstart$time & tabletofix$code==oldstart$code)),] #gets rid of old end behavior
                }  
                
                
                if (starts[nnn,1]<stopBegin & ends[nnn,1]>stopEnd & ends[nnn,2]!="OFF(end)"){ #SPAN_if start happens before bird goes OFF AND stop happens after bird comes back (OFF(end)), cut out middle and add new stop and start
                  starttimetoadd<-starts[nnn,]
                  starttimetoadd$time<-stopEnd
                  
                  stoptimetoadd<-ends[nnn,]
                  stoptimetoadd$time<-stopBegin
                  
                  
                  addstartstop<-rbind(starttimetoadd,stoptimetoadd)
                  
                  tabletofix<-rbind(tabletofix,addstartstop)
                }
                
                
                if (starts[nnn,1]>stopBegin & ends[nnn,1]<stopEnd & ends[nnn,2]!="OFF(end)"){ #MID_if start AND stop of a duration behavior occur while bird is OFF, delete
                  starttimetoeliminate<-starts[nnn,1]
                  stoptimetoeliminate<-ends[nnn,1]
                  startremove<-behav.obs[which(behav.obs$time == starttimetoeliminate),] 
                  endremove<-behav.obs[which(behav.obs$time == stoptimetoeliminate),] 
                  
                  gap.wise.misses2<-rbind(startremove,endremove)
                  
                  tabletofix<-tabletofix[-((match(gap.wise.misses2$time,tabletofix$time))),]
                }
                
                
              }# for every row in 'starts', which has all starts for a behavior, it calculates the difference between end and start and then adds this to 'duration'
            }    
          }
        }
        
      }
      
    } 
    
    #REMOVES behaviors that have a start and stop at the same time
    
    
    for (db2 in DurationBehaviors.nosp){  
      behav.observations <- tabletofix[grep(db2, tabletofix$code), ]
      startings<-behav.observations[grep("start",behav.observations$code),]
      endings<-behav.observations[grep("end",behav.observations$code),]
      
      startstoremove<-startings[(startings$time %in% endings$time),]
      endstoremove<-endings[(endings$time %in% startings$time),]
      kill<-rbind(startstoremove,endstoremove)
      tabletofix<-tabletofix[!(tabletofix$time %in% kill$time & tabletofix$code %in% kill$code),]
      
    }	
    
    tabletofix[order(tabletofix$time),]->table
    
    return(table)
  }  
  
  
  
  
  # setting wd, and performing initial naming tasks related to file output 
  ############################
  dir<-getwd()
  saveithere<-paste(dir,"/checks/",sep="")
  FileisNamed<-right(filename,46)
  
  name<-paste(FileisNamed,"pdf",sep=".")
  name<-gsub("/", ".", name, fixed = TRUE)
  name<-gsub(".csv", "", name, fixed = TRUE)
  name<-gsub("000000_000", "", name, fixed = TRUE)
  name<-gsub("_2000-01-01T","",name, fixed = TRUE)
  name<-gsub("ut_oneatatime.","",name, fixed = TRUE)
  name<-gsub("eatatime.","",name, fixed = TRUE)
  name<-gsub("sis.FixedCSVs.","",name, fixed = TRUE)
  
  
  name2<-gsub("/", ".", FileisNamed, fixed = TRUE)
  name2<-gsub(".csv", "", name2, fixed = TRUE)
  name2<-gsub("000000_000", "", name2, fixed = TRUE)
  name2<-gsub("_2000-01-01T","",name2, fixed = TRUE)
  name2<-gsub("ut_oneatatime.","",name2, fixed = TRUE)
  name2<-gsub("eatatime.","",name2, fixed = TRUE)
  name2<-gsub("sis.FixedCSVs.","",name2, fixed = TRUE)
  
  
  fullname<-paste(location,name,sep="/")
  ###########  
  #read in table
  read.table(filename,sep=",",header=TRUE)->table
  
  #Order selections according to begin time
  table[order(table$time),]->table
  
  
  
  #FIRST STEP PREPROCESSING
  
  #Early files used "Mpres" instead of "M pres"...these were incorrect and are removed in this step		
  rowswithMpresnospace<-which(table$code == "Mpres")
  if (length(rowswithMpresnospace)!=0){
    table<-table[-((which(table$code == "Mpres"))),] #Removes 'Mpres'
  }
  
  ####Remove spaces from 'code' identifiers
  table$code<-gsub(" ", "", table$code, fixed = TRUE)
  
  #Order selections according to begin time	
  table[order(table$time),]->table	
  
  #ELIMINATES DUPLICATE ROWS	
  table<-unique(table) 
  ###########
  #DATA CHECKING LOCATION, BEHAVIOR BY BEHAVIOR
  behaviortoinvestigate<-"OPMF"
  info<-table[grep(behaviortoinvestigate, table$code), ]
  rle(table[grep(behaviortoinvestigate, table$code), ][,2])
  ######################
  
  #Make repetitively pressed buttons (e.g. OPMW, OPAW, OPAC3) into duration of RM behaviors (e.g. RM1)
  
  flappinglist<-c("OPAW","OPAC3","OPAMW")	
  wingsflapping<-subset(table,code %in% flappinglist)	
  wingsflapping<-subset(wingsflapping,!duplicated(time))
  if(nrow(wingsflapping)>0){
    wingsflapping$timetonext<-c(diff(wingsflapping$time),1000)
    wingsflapping$keep<- ifelse((wingsflapping$code==(shift(wingsflapping$code,1L))) | (wingsflapping$code!=(shift(wingsflapping$code,1L)) & wingsflapping$timetonext>0.2),"yes","no")
    wingsflapping<-subset(wingsflapping,keep=="yes")
    wingsflapping<-wingsflapping[,c(1:3)]
  }
  
  
  taillist<-c("OPMT")
  tailflapping<-subset(table,code %in% taillist)	
  
  headlist<-c("OPAH")
  headflapping<-subset(table,code %in% headlist)
  
  torsolist<-c("OPAC1","OPAC2","OPATW","OPAMTT")
  torsoflapping<-subset(table, code %in% torsolist)
  torsoflapping<-subset(torsoflapping,!duplicated(time))
  if(nrow(torsoflapping)>0){
    torsoflapping$timetonext<-c(diff(torsoflapping$time),1000)
    torsoflapping$keep<- ifelse((torsoflapping$code==(shift(torsoflapping$code,1L))) | (torsoflapping$code!=(shift(torsoflapping$code,1L)) & torsoflapping$timetonext>0.2),"yes","no")
    torsoflapping<-subset(torsoflapping,keep=="yes")
    torsoflapping<-torsoflapping[,c(1:3)]
  }
  
  
  
  repeptiveornamentaccentuations<-list(wingsflapping,tailflapping,headflapping,torsoflapping)	
  
  durationofrepetbehavs<-table[1,]	
  durationofrepetbehavs<-durationofrepetbehavs[-1,]
  RepBehavLabels<-c("RepWing","RepTail","RepHead","RepTorso")
  
  for (RMs in 1:length(repeptiveornamentaccentuations)){	
    #flappers<-tailflapping
    
    newlabel<-RepBehavLabels[RMs]
    flappers<-repeptiveornamentaccentuations[[RMs]]
    flappers<-flappers[order(flappers$time),]
    
    if(nrow(flappers)>1){
      attach(flappers)
      
      timetonext<-c(diff(time),"FALSE")	
      TFseq<-timetonext<1 
      NewVec<-vector(length=length(TFseq))
      
      for (gg in 1:length(TFseq)) {
        
        xnow<-TFseq[gg]
        pre<-ifelse(gg>1,TFseq[(gg-1)],NA)
        nex<-ifelse(gg==length(TFseq),NA,TFseq[(gg+1)])
        
        if(is.na(pre)){  
          if(xnow==TRUE & is.na(pre)){
            newvalue<-"start"
          } else {
            if (xnow==FALSE & is.na(pre)){
              newvalue<-"off"
            }
          }
        } else {
          if (is.na(nex)){
            if (xnow==FALSE & (pre)==FALSE){
              newvalue<-"off"
            }
            if (xnow==FALSE & (pre)==TRUE){
              newvalue<-"end"
            }
          } else {
            if(xnow==TRUE & (nex)==TRUE & (pre)==FALSE){
              newvalue<-"start"
            } 
            if (xnow==TRUE & (nex)==FALSE & (pre)==TRUE){
              newvalue<-"on"
            }
            if(xnow==TRUE & (nex)==FALSE & (pre)==FALSE){
              newvalue<-"start"
            } 
            if (xnow==TRUE & (nex)==TRUE & (pre)==TRUE){
              newvalue<-"on"
            }
            if (xnow==FALSE & (nex)==TRUE & (pre)==TRUE){
              newvalue<-"end"
            }
            if (xnow==FALSE & (nex)==TRUE & (pre)==FALSE){
              newvalue<-"off"
            }
            if (xnow==FALSE & (nex)==FALSE & (pre)==FALSE){
              newvalue<-"off"
            }
            if (xnow==FALSE & (nex)==FALSE & (pre)==TRUE){
              newvalue<-"end"
            }
          }
        }
        NewVec[gg]<-newvalue
      }
      detach(flappers)
      
      flappers$times<-NewVec
      startsofrep<-subset(flappers,times=="start")
      if(nrow(startsofrep)>0){startsofrep$code<-paste(newlabel,"(start)",sep="")}
      endsofrep<-subset(flappers,times=="end")
      if(nrow(endsofrep)>0){endsofrep$code<-paste(newlabel,"(end)",sep="")}
      
      thesebehaviors<-rbind(startsofrep,endsofrep)
      thesebehaviors<-thesebehaviors[,-4]
      
      
      durationofrepetbehavs<-rbind(durationofrepetbehavs,thesebehaviors)
      durationofrepetbehavs<-durationofrepetbehavs[order(durationofrepetbehavs$time),]	
    }
  }	
  
  table<-rbind(table,durationofrepetbehavs)	
  
  
  #SECOND STEP PROCESSING (handling OFF periods)	
  
  
  #############################################################	
  #if last OFF (end) is the same as the END, this cuts off analysis at that point
  
  if(nrow(subset(table,code=="OFF(end)"))>0){
    if(isTRUE(all.equal(max((subset(table,code=="OFF(end)")[,1]),na.rm=TRUE),subset(table,code=="END")[,1],tolerance = 0.0009))){	  
      OriginalEndingTime<-	max((subset(table,code=="OFF(end)")[,1]),na.rm=TRUE)
      TimeBirdFlewOff<-tail(subset(table,code=="OFF(start)")[,1],1)
      FinishEnd<-data.frame("time"=TimeBirdFlewOff,"code"="END","class"=0)
      beforeflewaway<-subset(table,time<=TimeBirdFlewOff & code !="OFF(start)")
      previousOFFs<-subset(table,time<TimeBirdFlewOff & code =="OFF(start)")
      
      AtEnd<-subset(table,time==OriginalEndingTime & code !="END" & code != "OFF(end)")
      if(nrow(AtEnd)>0){
        AtEnd$time<-TimeBirdFlewOff
      }
      
      
      tabletable<-rbind(beforeflewaway,previousOFFs,AtEnd,FinishEnd)
      table<-tabletable[order(tabletable$time),]
    }
  }	
  
  #Finds END time and uses to calculate TotalTime
  TotalTime<-table[which(table$code == "END"), ][,1]	
  
  ######if first OFF (start) is at t=0, this cuts off analysis during this first 'missing' bit
  BirdWasOffatBeginning<-0
  TimeBirdFlewON<-0
  if(nrow(subset(table,code=="OFF(start)"))>0){
    if(isTRUE(all.equal(min((subset(table,code=="OFF(start)")[,1]),na.rm=TRUE),0,tolerance = 0.0009))){	  
      BirdWasOffatBeginning<-1
      OriginalStartingTime<-	0
      TimeBirdFlewON<-subset(table,code=="OFF(end)")[1,1]
      #FinishEnd<-data.frame("time"=TimeBirdFlewOff,"code"="END","class"=0)
      afterflewON<-subset(table,time>=TimeBirdFlewON & code !="OFF(end)")
      subsequentOFFs<-subset(table,time>TimeBirdFlewON & code =="OFF(end)")
      
      AtStart<-subset(table,time==0 &  code != "OFF(start)")
      if(nrow(AtStart)>0){
        AtStart$time<-TimeBirdFlewON
      }
      
      
      tabletabletable<-rbind(afterflewON,subsequentOFFs,AtStart)
      table<-tabletabletable[order(tabletabletable$time),]
      
      
    }
  }	
  
  
  OFFtimes<-table[grep("OFF", table$code), ][,1]
  OFFstartTimes<-subset(table,code=="OFF(start)")[,1]  
  OFFendTimes<-subset(table,code=="OFF(end)")[,1]  
  
  #############################################################	
  
  #PLOT A	
  
  
  table<-table[order(table$time),]	
  table<-fixOFFs(table)	
  
  #Adds O2 for every changed in O3
  O2additions<-rbind(table[which(table$code == "O3(start)"),],table[which(table$code == "O3(end)"),])
  O2additions$code<-(gsub("O3(end)", "O2", gsub("O3(start)", "O2", O2additions$code, fixed = TRUE), fixed = TRUE)) #substituting
  table<-rbind(table,O2additions)	
  
  
  #THIRD STEP PROCESSING (handling body orientation issues, e.g. O5s and O6s, O4 contingencies)	
  #############################################################			
  #############################################################	
  #############################################################	
  #Finds END time and uses to calculate TotalTime
  TotalTime<-table[which(table$code == "END"), ][,1]
  
  
  #If there are no 05s, then the bird was upright the entire time
  if(nrow(table[which(table$code == "O5(start)"), ])!=0){ 
    firstO5<-min(table[which(table$code == "O5(start)"), ][,1],na.rm = TRUE) #Takes times of all O5 starts and finds first (minimum)
    lastO5<-max(table[which(table$code == "O5(end)"), ][,1],na.rm = TRUE) #Takes times of all O5 starts and finds last (maximum)
    # O5starttable<-rbind(table[which(table$code == "O5(end)"),],table[which(table$code == "O5(start)"),],table[which(table$code == "END"),])#binds all O5 starts, ends, and the Trial End
    O5starttable<-table[grep("O5", table$code), ]
    endtime<-table[which(table$code == "END"),][,1]
    
    O6starttable<-O5starttable
    # O6starttable<-O6starttable[-(which(O6starttable$code=="O5(start)" & O6starttable$time %in% OFFendTimes)),]
    # O6starttable<-O6starttable[-(which(O6starttable$code=="O5(end)" & O6starttable$time %in% OFFstartTimes)),]
    
    if (BirdWasOffatBeginning>0){
      if (isTRUE(all.equal(firstO5,TimeBirdFlewON,tolerance = 0.0009))){
        if (isTRUE(all.equal(lastO5,endtime,tolerance = 0.0009))) { #If O5 begins and ends a trial, substiutions should not be made on those obs
          O6starttable<-O6starttable[!(O6starttable$code == "O5(start)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, TimeBirdFlewON, tolerance = 0.0009)))), ]
          O6starttable<-O6starttable[!(O6starttable$code == "O5(end)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]
          O6starttable$code<-(gsub("O5(end)", "O6(start)", O6starttable$code, fixed = TRUE)) 
          O6starttable$code<-(gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE))#substituting O6starts for O5 ends, etc
          
          O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
          
          
          
        } else { #If O5 begins a trial, substitutions should not be made on that first obs, and O6(end) should be at the end
          O6starttable<-O6starttable[!(O6starttable$code == "O5(start)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, TimeBirdFlewON, tolerance = 0.0009)))), ]
          O6starttable$code<-(gsub("O5(end)", "O6(start)", O6starttable$code, fixed = TRUE)) 
          O6starttable$code<-(gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE))#substituting O6starts for O5 ends, etc
          O6starttable$code<-(gsub("END", "O6(end)", O6starttable$code, fixed = TRUE)) #substituting O6starts for O5 ends, etc
          
          O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
          
          O6ends<-O6starttable[1,]
          O6ends$time<-endtime
          O6ends$code<-"O6(end)"
          
          O6starttable<-rbind(O6starttable,O6ends)
          
        }
      } else { #IF O5 doesn't start the trial, O6 must
        if (isTRUE(all.equal(lastO5,endtime,tolerance = 0.0009))) {
          O6startsatzero<-O6starttable[1,]
          O6startsatzero$time<-TimeBirdFlewON
          O6startsatzero$code<-"O6(start)" 
          O6starttable<-rbind(O6startsatzero,O6starttable)#This (+3 lines above) adds an O6(start) at t=0
          
          O6starttable<-O6starttable[!(O6starttable$code == "O5(end)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]#Don't replace O5(end) if it's at the END
          O6starttable$code<-(gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE))
          O6starttable$code<-(gsub("O5(end)", "O6(start)", O6starttable$code, fixed = TRUE))
          
          O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
          
          
        } else { #When O5 neither starts the trial nor ends the trial, all subs can be made, END can be turned into O6(end), and O6(start) can be added at t=0
          O6startsatzero<-O6starttable[1,]
          O6startsatzero$time<-TimeBirdFlewON
          O6startsatzero$code<-"O6(start)" 
          O6starttable<-rbind(O6startsatzero,O6starttable)#This (+3 lines above) adds an O6(start) at t=0
          
          O6starttable$code<-(gsub("O5(end)", "O6(start)", gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE), fixed = TRUE)) #substituting O6starts for O5 ends, etc
          O6starttable$code<-(gsub("END", "O6(end)", O6starttable$code, fixed = TRUE)) #If 05 doesn't end the trial, O6 must
          
          O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
          
          O6ends<-O6starttable[1,]
          O6ends$time<-endtime
          O6ends$code<-"O6(end)"
          
          O6starttable<-rbind(O6starttable,O6ends)
          
          
        }
      }
    } else {
      if (isTRUE(all.equal(firstO5,0,tolerance = 0.0009))){
        if (isTRUE(all.equal(lastO5,endtime,tolerance = 0.0009))) { #If O5 begins and ends a trial, substiutions should not be made on those obs
          O6starttable<-O6starttable[!(O6starttable$code == "O5(start)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, 0, tolerance = 0.0009)))), ]
          O6starttable<-O6starttable[!(O6starttable$code == "O5(end)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]
          O6starttable$code<-(gsub("O5(end)", "O6(start)", O6starttable$code, fixed = TRUE)) 
          O6starttable$code<-(gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE))#substituting O6starts for O5 ends, etc
          
          O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
          
          
        } else { #If O5 begins a trial, substitutions should not be made on that first obs, and O6(end) should be at the end
          O6starttable<-O6starttable[!(O6starttable$code == "O5(start)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, 0, tolerance = 0.0009)))), ]
          O6starttable$code<-(gsub("O5(end)", "O6(start)", O6starttable$code, fixed = TRUE)) 
          O6starttable$code<-(gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE))#substituting O6starts for O5 ends, etc
          O6starttable$code<-(gsub("END", "O6(end)", O6starttable$code, fixed = TRUE)) #substituting O6starts for O5 ends, etc
          
          O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
          
          O6ends<-O6starttable[1,]
          O6ends$time<-endtime
          O6ends$code<-"O6(end)"
          
          O6starttable<-rbind(O6starttable,O6ends)
          
        }
      } else { #IF O5 doesn't start the trial, O6 must
        if (isTRUE(all.equal(lastO5,endtime,tolerance = 0.0009))) {
          O6startsatzero<-O6starttable[1,]
          O6startsatzero$time<-0
          O6startsatzero$code<-"O6(start)" 
          O6starttable<-rbind(O6startsatzero,O6starttable)#This (+3 lines above) adds an O6(start) at t=0
          
          O6starttable<-O6starttable[!(O6starttable$code == "O5(end)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]#Don't replace O5(end) if it's at the END
          O6starttable$code<-(gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE))
          O6starttable$code<-(gsub("O5(end)", "O6(start)", O6starttable$code, fixed = TRUE))
          
          O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
          
          
        } else { #When O5 neither starts the trial nor ends the trial, all subs can be made, END can be turned into O6(end), and O6(start) can be added at t=0
          O6startsatzero<-O6starttable[1,]
          O6startsatzero$time<-0
          O6startsatzero$code<-"O6(start)" 
          O6starttable<-rbind(O6startsatzero,O6starttable)
          #This (+3 lines above) adds an O6(start) at t=0
          
          O6ends<-O6starttable[1,]
          O6ends$time<-endtime
          O6ends$code<-"O6(end)"
          
          O6starttable$code<-(gsub("O5(end)", "O6(start)", gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE), fixed = TRUE)) #substituting O6starts for O5 ends, etc
          O6starttable$code<-(gsub("END", "O6(end)", O6starttable$code, fixed = TRUE)) #If 05 doesn't end the trial, O6 must
          
          O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
          
          O6starttable<-rbind(O6starttable,O6ends)
          
          
        }
      }
    }
  } else {
    if(BirdWasOffatBeginning>0){
      O6starttable<-data.frame("time"=c(TimeBirdFlewON,TotalTime),"code"=c("O6(start)","O6(end)"),"class"=c(5,5))
    } else {
      O6starttable<-data.frame("time"=c(0,TotalTime),"code"=c("O6(start)","O6(end)"),"class"=c(5,5))
    }
  } 
  
  table2<-rbind(table,O6starttable)
  
  
  table<-table2 #keeps table our working, modified dataframe
  table<-unique(table) #ELIMINATES DUPLICATE ROWS
  table<-table[order(table$time),]	
  table<-fixOFFs(table)	#Fixes OFF issues
  
  #PLOT D
  #PLOTB<-table
  
  
  #Order selections according to begin time
  table<-table[order(table$time),]
  
  O4active<-rbind(table[which(table$code == "O4(start)"),],table[which(table$code == "O4(end)"),])
  O4active<-O4active[order(O4active$time),]#Order selections according to begin time
  O4starts<-table[which(table$code == "O4(start)"),]
  O4ends<-table[which(table$code == "O4(end)"),]
  O4bouts<-length(table[which(table$code == "O4(start)"),][,1])
  
  
  compareO5<-table[which(table$code == "O5(start)"),]
  compareO5end<-table[which(table$code == "O5(end)"),]
  compareO6<-table[which(table$code == "O6(start)"),]
  compareO6end<-table[which(table$code == "O6(end)"),]
  
  #This if addresses all postural (O5, O6,O7,O8) issues (if there are at least one or more O4s)
  if(nrow(O4starts)==nrow(O4ends) & O4bouts>0){   #IF BRANCH IS NOT VERTICAL THE ENTIRE TIME, WE NEED THIS CONTINGENCY (ELSE FOR CASES WHEN WE NEVER HAVE ANY O4 PERIODS/BOUTS)
    
    #
    if(nrow(compareO5)==nrow(compareO5end) & nrow(compareO5)!=0 & nrow(compareO5end)!=0) { #ONLY DO THIS IF EQUAL STARTS/STOPS FOR O5 (and !=0)
      for(j in 1:length(compareO5[,1])){
        
        for(i in 1:O4bouts){
          start<-O4starts[i,1]
          stop<-O4ends[i,1]
          startrow<-O4starts[i,]
          stoprow<-O4ends[i,]
          
          if(compareO5$time[j]==start & compareO5end$time[j]<stop ){#S3_This part handles O5s that start exactly as O4 starts and end before each O4 ends.
            oldendO5<-compareO5end[j,]
            oldstartO5<-compareO5[j,]
            
            O7start<-oldstartO5
            O7start$code="O7(start)"
            O7end<-oldendO5
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O7O8s<-newO7s
            O7O8s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
            
          }  
          
          if(compareO5$time[j]<start & compareO5end$time[j]==stop ){#S2_This part handles O5s that start before each O4 starts and end as each O4 ends.
            oldendO5<-compareO5end[j,]
            oldstartO5<-compareO5[j,]
            
            O8start<-oldstartO5
            O8start$code="O8(start)"
            O8end<-startrow
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O7start<-startrow
            O7start$code="O7(start)"
            O7end<-stoprow
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O7O8s<-rbind(newO7s,newO8s)
            O7O8s$class<-paste(i,j)
            table<-rbind(table,newO7s) #adds new O7s and O8s
            
          }  
          
          if(compareO5$time[j]<start & compareO5end$time[j]<stop & compareO5end$time[j]>start ){#S1_This part handles all O5s that start before each O4 starts and end while O4 is on
            oldendO5<-compareO5end[j,]
            oldstartO5<-compareO5[j,]
            
            O8start<-oldstartO5
            O8start$code="O8(start)"
            O8end<-startrow
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O7start<-startrow
            O7start$code="O7(start)"
            O7end<-oldendO5
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O7O8s<-rbind(newO7s,newO8s)
            O7O8s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
          }  
          
          if(compareO5$time[j]>start & compareO5$time[j]<stop & compareO5end$time[j]==stop ){#E4_This part handles all O5s that start while O4 and end exactly as O4 ends
            oldendO5<-compareO5end[j,]
            oldstartO5<-compareO5[j,]
            
            Nope<-oldstartO5[-1,]
            newO8s<-Nope
            
            O7start<-oldstartO5
            O7start$code="O7(start)"
            O7end<-oldendO5
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O7O8s<-rbind(newO7s,newO8s)
            O7O8s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
            
          }  
          
          if(compareO5$time[j]==start & compareO5end$time[j]==stop ){#E3_This part handles all O5s that start and end exactly as O4 starts and ends
            oldendO5<-compareO5end[j,]
            oldstartO5<-compareO5[j,]
            
            Nope<-oldstartO5[-1,]
            newO8s<-Nope
            
            O7start<-oldstartO5
            O7start$code="O7(start)"
            O7end<-oldendO5
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O7O8s<-rbind(newO7s,newO8s)
            O7O8s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
            
          }  
          
          if(compareO5$time[j]==start & compareO5end$time[j]>stop ){#E2_This part handles O5s that start as O4 starts and end after O4 ends.
            oldendO5<-compareO5end[j,]
            oldstartO5<-compareO5[j,]
            
            O8start<-stoprow
            O8start$code="O8(start)"
            O8end<-oldendO5
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O7start<-oldstartO5
            O7start$code="O7(start)"
            O7end<-stoprow
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O7O8s<-rbind(newO7s,newO8s)
            O7O8s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
          } 
          
          if(compareO5$time[j]>start & compareO5$time[j]<stop & compareO5end$time[j]>stop ){#E1_This part handles O5s that start while O4 is on, then end after O4 is off
            oldendO5<-compareO5end[j,]
            oldstartO5<-compareO5[j,]
            
            O8start<-stoprow
            O8start$code="O8(start)"
            O8end<-oldendO5
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O7start<-oldstartO5
            O7start$code="O7(start)"
            O7end<-stoprow
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O7O8s<-rbind(newO7s,newO8s)
            O7O8s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
            
          }  
          
          if (compareO5$time[j]<start & compareO5end$time[j]>stop){ #SPAN_if start happens before bird goes OFF AND stop happens after bird comes back (OFF(end)), cut out middle and add new stop and start
            oldendO5<-compareO5end[j,]
            oldstartO5<-compareO5[j,]
            
            O8start<-rbind(oldstartO5,stoprow)
            O8start$code="O8(start)"
            O8end<-rbind(startrow,oldendO5)
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O7start<-startrow
            O7start$code="O7(start)"
            O7end<-stoprow
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O7O8s<-rbind(newO7s,newO8s)
            O7O8s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
            
          }
          
          if (compareO5$time[j]>start & compareO5end$time[j]<stop){ #MID_if start AND stop of a duration behavior occur while bird is OFF, delete
            oldendO5<-compareO5end[j,]
            oldstartO5<-compareO5[j,]
            
            O7start<-oldstartO5
            O7start$code="O7(start)"
            O7end<-oldendO5
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O7O8s<-newO7s
            O7O8s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
          }
          
          if ((compareO5$time[j]<start & compareO5end$time[j]<start) | (compareO5$time[j]>stop & compareO5end$time[j]>stop)){ #PRE OR POST_if start AND stop of a duration behavior occur while bird is not on O4
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            O7start<-rbind(oldstartO6,oldendO6)
            O7start$code=c("O7(start)","O7(end)")
            newO7s<-O7start
            
            table<-rbind(table,newO7s)
          }   
          
        }
      }
    }
    if(nrow(compareO6)==nrow(compareO6end) & nrow(compareO6)!=0 & nrow(compareO6end)!=0) { #ONLY DO THIS IF EQUAL STARTS/STOPS FOR O6 (and !=0)
      for(j in 1:length(compareO6[,1])){
        
        for(i in 1:O4bouts){
          start<-O4starts[i,1]
          stop<-O4ends[i,1]
          startrow<-O4starts[i,]
          stoprow<-O4ends[i,]
          
          
          
          if(compareO6$time[j]==start & compareO6end$time[j]<stop ){#S3_This part handles O6s that start exactly as O4 starts and end before each O4 ends.
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            
            
            
          }  
          if(compareO6$time[j]<start & compareO6end$time[j]==stop ){#S2_This part handles O6s that start before each O4 starts and end as each O4 ends.
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            
            O7start<-oldstartO6
            O7start$code="O7(start)"
            O7end<-startrow
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O8start<-startrow
            O8start$code="O8(start)"
            O8end<-stoprow
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O8O7s<-rbind(newO8s,newO7s)
            O8O7s$class<-paste(i,j)
            table<-rbind(table,newO7s) #adds new O8s and O7s
            
          }  
          
          if(compareO6$time[j]<start & compareO6end$time[j]<stop & compareO6end$time[j]>start ){#S1_This part handles all O6s that start before each O4 starts and end while O4 is on
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            
            O7start<-oldstartO6
            O7start$code="O7(start)"
            O7end<-startrow
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O8start<-startrow
            O8start$code="O8(start)"
            O8end<-oldendO6
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O8O7s<-rbind(newO8s,newO7s)
            O8O7s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
          }  
          
          if(compareO6$time[j]>start & compareO6$time[j]<stop & compareO6end$time[j]==stop ){#E4_This part handles all O6s that start while O4 and end exactly as O4 ends
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            
            
            
          }  
          
          if(compareO6$time[j]==start & compareO6end$time[j]==stop ){#E3_This part handles all O6s that start and end exactly as O4 starts and ends
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            
            
            
          }  
          
          if(compareO6$time[j]==start & compareO6end$time[j]>stop ){#E2_This part handles O6s that start as O4 starts and end after O4 ends.
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            
            O7start<-stoprow
            O7start$code="O7(start)"
            O7end<-oldendO6
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O8start<-oldstartO6
            O8start$code="O8(start)"
            O8end<-stoprow
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O8O7s<-rbind(newO8s,newO7s)
            O8O7s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
          } 
          
          if(compareO6$time[j]>start & compareO6$time[j]<stop & compareO6end$time[j]>stop ){#E1_This part handles O6s that start while O4 is on, then end after O4 is off
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            
            O7start<-stoprow
            O7start$code="O7(start)"
            O7end<-oldendO6
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O8start<-oldstartO6
            O8start$code="O8(start)"
            O8end<-stoprow
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O8O7s<-rbind(newO8s,newO7s)
            O8O7s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
            
          }  
          
          if (compareO6$time[j]<start & compareO6end$time[j]>stop){ #SPAN_if start happens before bird goes O4 AND stop happens after O4 goes OFF, cut out middle and add new stop and start
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            
            O7start<-rbind(oldstartO6,stoprow)
            O7start$code="O7(start)"
            O7end<-rbind(startrow,oldendO6)
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O8start<-startrow
            O8start$code="O8(start)"
            O8end<-stoprow
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O8O7s<-rbind(newO8s,newO7s)
            O8O7s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
            
          }
          
          if (compareO6$time[j]>start & compareO6end$time[j]<stop){ #MID_if start AND stop of a duration behavior occur while bird is OFF, delete
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            
            
          }
          
          if ((compareO6$time[j]<start & compareO6end$time[j]<start) | (compareO6$time[j]>stop & compareO6end$time[j]>stop)){ #PRE OR POST_if start AND stop of O6 occurs completely before or after an O4 bout
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            O7start<-rbind(oldstartO6,oldendO6)
            O7start$code=c("O7(start)","O7(end)")
            newO7s<-O7start
            
            table<-rbind(table,newO7s)
          }
        }
      }
    } 
    
  }	    
  ################################################
  #This if addresses all postural (O5, O6,O7,O8) issues (if there are no O4s)
  if(nrow(O4starts)==nrow(O4ends) & O4bouts==0){
    oldendO5<-compareO5end
    oldstartO5<-compareO5
    
    oldendO6<-compareO6end
    oldstartO6<-compareO6
    
    if(nrow(oldstartO5)>0){
      O8start<-oldstartO5
      O8start$code="O8(start)"
      O8end<-oldendO5
      O8end$code<-"O8(end)"
      newO8s<-rbind(O8start,O8end)
      newO8s$class<-5
    } else {
      newO8s<-oldendO5
      if(nrow(newO8s)>0){
        newO8s$class<-5
      }
    }
    
    if(nrow(oldstartO6)>0){
      O7start<-oldstartO6
      O7start$code="O7(start)"
      O7end<-compareO6end
      O7end$code<-"O7(end)"
      newO7s<-rbind(O7start,O7end)
      newO7s$class<-5
    } else {
      newO7s<-oldstartO6
      if(nrow(newO7s)>0){
        newO7s$class<-5
      }
      
    }
    
    O7O8s<-rbind(newO7s,newO8s)
    O7O8s$class<-5
    table<-rbind(table,O7O8s) 
  }
  
  midcheck<-table
  table<-unique(table) #ELIMINATES DUPLICATE ROWS
  table<-table[order(table$time),]	
  
  
  #ORDER (start/stop) FIXER FOR DURATION BEHAVIOR
  sevens<-table[grep("O7",table$code),] #O7 values
  if(length(sevens[,1])>0) {sevens$class<-"77"} #Standardize O7 class values
  sevens<-unique(sevens)
  
  
  sevens$keep<- ifelse((sevens$code==(data.table::shift(sevens$code,n=1L,type="lead"))) ,"no","yes")
  sevens<-subset(sevens,keep=="yes" | is.na(keep))[,-4]
  sevens<-subset(sevens,!(is.na(time)))
  
  table<-table[!(grepl("O7", table$code)),] #get rid of O7 valuse
  
  
  
  
  
  
  
  
  
  table<-rbind(table,sevens)
  
  
  
  
  table<-fixOFFs(table)	#Fixes OFF issues
  
  # 
  # 	AllO7<-table[which(table$code=="O7(start)" | table$code=="O7(end)" ),]
  # 	AllO7<-AllO7[order(AllO7$time),]	
  # 	
  # 	AllO8<-table[which(table$code=="O8(start)" | table$code=="O8(end)"),]
  # 	AllO8<-AllO8[order(AllO8$time),]	
  # 	
  # 	noO7O8<-table[-which(table$code=="O7(start)" | table$code=="O7(end)" | table$code=="O8(start)" | table$code=="O8(end)"),]
  # 	table<-rbind(noO7O8,AllO7,AllO8)
  #########################################################################	
  #If there are no O7s, then the bird was O8 the entire time
  TTolerance <- 0.0009
  
  
  if(nrow(table[which(table$code == "O7(start)"), ])!=0){ 
    firstO7<-min(table[which(table$code == "O7(start)"), ][,1],na.rm = TRUE) #Takes times of all O7 starts and finds first (minimum)
    lastO7<-max(table[which(table$code == "O7(end)"), ][,1],na.rm = TRUE) #Takes times of all O7 starts and finds last (maximum)
    # O7starttable<-rbind(table[which(table$code == "O7(end)"),],table[which(table$code == "O7(start)"),],table[which(table$code == "END"),])#binds all O7 starts, ends, and the Trial End
    O7starttable<-table[grep("O7", table$code), ]
    endtime<-table[which(table$code == "END"),][,1]
    
    O8starttable<-O7starttable
    # O8starttable<-O8starttable[-(which(O8starttable$code=="O7(start)" & O8starttable$time %in% OFFendTimes)),]
    # O8starttable<-O8starttable[-(which(O8starttable$code=="O7(end)" & O8starttable$time %in% OFFstartTimes)),]
    
    if (BirdWasOffatBeginning>0){
      if (isTRUE(all.equal(firstO7,TimeBirdFlewON,tolerance = 0.0009))){
        if (isTRUE(all.equal(lastO7,endtime,tolerance = 0.0009))) { #If O7 begins and ends a trial, substiutions should not be made on those obs
          O8starttable<-O8starttable[!(O8starttable$code == "O7(start)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, TimeBirdFlewON, tolerance = 0.0009)))), ]
          O8starttable<-O8starttable[!(O8starttable$code == "O7(end)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]
          O8starttable$code<-(gsub("O7(end)", "O8(start)", O8starttable$code, fixed = TRUE)) 
          O8starttable$code<-(gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE))#substituting O8starts for O7 ends, etc
          
          O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
          
          
          
        } else { #If O7 begins a trial, substitutions should not be made on that first obs, and O8(end) should be at the end
          O8starttable<-O8starttable[!(O8starttable$code == "O7(start)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, TimeBirdFlewON, tolerance = 0.0009)))), ]
          O8starttable$code<-(gsub("O7(end)", "O8(start)", O8starttable$code, fixed = TRUE)) 
          O8starttable$code<-(gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE))#substituting O8starts for O7 ends, etc
          O8starttable$code<-(gsub("END", "O8(end)", O8starttable$code, fixed = TRUE)) #substituting O8starts for O7 ends, etc
          
          O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
          
          O8ends<-O8starttable[1,]
          O8ends$time<-endtime
          O8ends$code<-"O8(end)"
          
          O8starttable<-rbind(O8starttable,O8ends)
          
        }
      } else { #IF O7 doesn't start the trial, O8 must
        if (isTRUE(all.equal(lastO7,endtime,tolerance = 0.0009))) {
          O8startsatzero<-O8starttable[1,]
          O8startsatzero$time<-TimeBirdFlewON
          O8startsatzero$code<-"O8(start)" 
          O8starttable<-rbind(O8startsatzero,O8starttable)#This (+3 lines above) adds an O8(start) at t=0
          
          O8starttable<-O8starttable[!(O8starttable$code == "O7(end)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]#Don't replace O7(end) if it's at the END
          O8starttable$code<-(gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE))
          O8starttable$code<-(gsub("O7(end)", "O8(start)", O8starttable$code, fixed = TRUE))
          
          O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
          
          
        } else { #When O7 neither starts the trial nor ends the trial, all subs can be made, END can be turned into O8(end), and O8(start) can be added at t=0
          O8startsatzero<-O8starttable[1,]
          O8startsatzero$time<-TimeBirdFlewON
          O8startsatzero$code<-"O8(start)" 
          O8starttable<-rbind(O8startsatzero,O8starttable)#This (+3 lines above) adds an O8(start) at t=0
          
          O8starttable$code<-(gsub("O7(end)", "O8(start)", gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE), fixed = TRUE)) #substituting O8starts for O7 ends, etc
          O8starttable$code<-(gsub("END", "O8(end)", O8starttable$code, fixed = TRUE)) #If 05 doesn't end the trial, O8 must
          
          O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
          
          O8ends<-O8starttable[1,]
          O8ends$time<-endtime
          O8ends$code<-"O8(end)"
          
          O8starttable<-rbind(O8starttable,O8ends)
          
          
        }
      }
    } else {
      if (isTRUE(all.equal(firstO7,0,tolerance = 0.0009))){
        if (isTRUE(all.equal(lastO7,endtime,tolerance = 0.0009))) { #If O7 begins and ends a trial, substiutions should not be made on those obs
          O8starttable<-O8starttable[!(O8starttable$code == "O7(start)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, 0, tolerance = 0.0009)))), ]
          O8starttable<-O8starttable[!(O8starttable$code == "O7(end)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]
          O8starttable$code<-(gsub("O7(end)", "O8(start)", O8starttable$code, fixed = TRUE)) 
          O8starttable$code<-(gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE))#substituting O8starts for O7 ends, etc
          
          O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
          
          
        } else { #If O7 begins a trial, substitutions should not be made on that first obs, and O8(end) should be at the end
          O8starttable<-O8starttable[!((O8starttable$code == "O7(end)" & abs(O8starttable$time - endtime) < TTolerance)), ]
          O8starttable$code<-(gsub("O7(end)", "O8(start)", O8starttable$code, fixed = TRUE)) 
          O8starttable$code<-(gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE))#substituting O8starts for O7 ends, etc
          O8starttable$code<-(gsub("END", "O8(end)", O8starttable$code, fixed = TRUE)) #substituting O8starts for O7 ends, etc
          
          O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
          
          O8ends<-O8starttable[1,]
          O8ends$time<-endtime
          O8ends$code<-"O8(end)"
          
          O8starttable<-rbind(O8starttable,O8ends)
          
        }
      } else { #IF O7 doesn't start the trial, O8 must
        if (isTRUE(all.equal(lastO7,endtime,tolerance = 0.0009))) {
          O8startsatzero<-O8starttable[1,]
          O8startsatzero$time<-0
          O8startsatzero$code<-"O8(start)" 
          O8starttable<-rbind(O8startsatzero,O8starttable)#This (+3 lines above) adds an O8(start) at t=0
          
          O8starttable<-O8starttable[!(O8starttable$code == "O7(end)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]#Don't replace O7(end) if it's at the END
          O8starttable$code<-(gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE))
          O8starttable$code<-(gsub("O7(end)", "O8(start)", O8starttable$code, fixed = TRUE))
          
          O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
          
          
        } else { #When O7 neither starts the trial nor ends the trial, all subs can be made, END can be turned into O8(end), and O8(start) can be added at t=0
          O8startsatzero<-O8starttable[1,]
          O8startsatzero$time<-0
          O8startsatzero$code<-"O8(start)" 
          O8starttable<-rbind(O8startsatzero,O8starttable)
          #This (+3 lines above) adds an O8(start) at t=0
          
          O8ends<-O8starttable[1,]
          O8ends$time<-endtime
          O8ends$code<-"O8(end)"
          
          O8starttable$code<-(gsub("O7(end)", "O8(start)", gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE), fixed = TRUE)) #substituting O8starts for O7 ends, etc
          O8starttable$code<-(gsub("END", "O8(end)", O8starttable$code, fixed = TRUE)) #If 05 doesn't end the trial, O8 must
          
          O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
          
          O8starttable<-rbind(O8starttable,O8ends)
          
          
        }
      }
    }
  } else {
    if(BirdWasOffatBeginning>0){
      O8starttable<-data.frame("time"=c(TimeBirdFlewON,TotalTime),"code"=c("O8(start)","O8(end)"),"class"=c(5,5))
    } else {
      O8starttable<-data.frame("time"=c(0,TotalTime),"code"=c("O8(start)","O8(end)"),"class"=c(5,5))
    }
  } 
  
  
  
  table5<-rbind(table,O8starttable)	
  
  table<-table5 #keeps table our working, modified dataframe
  table<-unique(table) #ELIMINATES DUPLICATE ROWS
  table<-table[order(table$time),]	
  
  
  #table<-fixOFFs(table)
  
  
  
  #############################################################	
  
  ##########REMOVES instances where there are multiple *starts* in a row
  # 	OrientationDurations<-c("O3","O4","O5","O6","O7","O8")	
  # 	
  # 	for (OD2 in OrientationDurations){  
  # 	  behav.observations <- table[grep(OD2, table$code), ]
  # 	  goo<-rle(behav.observations$code)
  # 	  
  # 	  firstone2<-behav.observations[(which(rep(goo$lengths>1,times=goo$lengths))),]
  # 	  firstone<-firstone2[1,]
  # 	  alternators<-behav.observations[-(which(rep(goo$lengths>1,times=goo$lengths))),]
  # 	  totalsequence<-rbind(firstone,alternators)
  # 	  totalsequence<-totalsequence[order(totalsequence$time),]    
  # 	  
  # 	
  # 	  kill<-firstone2[-1,]
  # 	  table<-table[!(table$time %in% kill$time & table$code %in% kill$code),]
  # 	  
  # 	}	
  ######################################################
  ####################################################
  #ELIMINATE FAULTILY CREATED DEPENDENT VARIABLES (e.g. O7(end) when 08(start) IF O8(start) = 0 )
  EndList<-c("BP1(end)","BP3(end)",
             "SS1(end)","SS2(end)",
             "O3(end)","O4(end)",
             "O5(end)","OFF(end)",
             "O6(end)","O7(end)","O8(end)",
             "OPMH.D(end)","OPMB.D(end)","OPMF.D(end)","OPMW.D(end)","OPMT.D(end)",
             "RM1(end)","RM2(end)","RM3(end)","RM4(end)",
             "MO2(end)","PU1(end)")
  
  for(xx in EndList){ #Nice little function that eliminates any "end" variables if they occur at t=0
    times2<-subset(table,code==xx)[,1]
    if(length(times2)>0){ 
      if(min(times2)==0){
        table<-table[-which(table$code == xx & table$time==0),]
      }
    }
  }
  
  for(xx in EndList){ #Nice little function that eliminates any "end" variables if they occur at the moment the bird flies on screen for the first time
    times<-subset(table,code==xx)[,1]
    if(length(times)>0){    
      if(min(times)==TimeBirdFlewON){
        table<-table[-which(table$code == xx & table$time==TimeBirdFlewON),]
      }
    }
  }	
  
  StartList2<-c("BP1(start)","BP3(start)",
                "SS1(start)","SS2(start)",
                "O3(start)","O4(start)",
                "O5(start)","OFF(end)",
                "O6(start)","O7(start)","O8(start)",
                "OPMH.D(start)","OPMB.D(start)","OPMF.D(start)","OPMW.D(start)","OPMT.D(start)",
                "RM1(start)","RM2(start)","RM3(start)","RM4(start)",
                "MO2(start)","PU1(start)")
  
  StartList<-c("BP1(start)","BP3(start)",
               "SS1(start)","SS2(start)",
               "O3(start)","O4(start)",
               "O5(start)",
               "O6(start)","O7(start)","O8(start)",
               "OPMH.D(start)","OPMB.D(start)","OPMF.D(start)","OPMW.D(start)","OPMT.D(start)",
               "RM1(start)","RM2(start)","RM3(start)","RM4(start)",
               "MO2(start)","PU1(start)")
  
  table<-fixOFFs(table)	#Fixes OFF issues
  
  
  for(xj in StartList){ #Nice little function that eliminates any "start" variables if they occur at the moment the bird flies off screen
    starttimes<-subset(table,code==xj)[,1]
    flyofftimes<-subset(table,code=="OFF(start)")[,1]
    if(length(starttimes)>0 & length(flyofftimes)){   
      table<-table[!(table$time %in% flyofftimes & table$code==xj),]        
    }
  }
  ############################################
  ###############################
  secondcheck<-table
  
  
  #############################################################  
  #OFF PERIODS---HOW MANY AND HOW MANY BEHAVIORS WERE MEASURED WHILE "OFF"
  
  happenedwhileoff<-data.frame()
  offstarts<-nrow(subset(table,code=="OFF(start)"))
  offstops<-nrow(subset(table,code=="OFF(end)"))
  off.ss.eq<-(offstarts==offstops)
  for (v in 1:nrow(subset(table,code=="OFF(start)"))){
    stopBegin<-subset(table,code=="OFF(start)")[v,1]
    stopEnd<-subset(table,code=="OFF(end)")[v,1]
    
    gap.wise.misses<-table[which(table$time < stopEnd & table$time > stopBegin),]
    happenedwhileoff<-rbind(happenedwhileoff,gap.wise.misses)
  }
  
  happenedwhileoff<-paste(happenedwhileoff$time,collapse="|")
  # 	
  #############################################################	
  table<-unique(table) #ELIMINATES DUPLICATE ROWS	
  
  
  
  
  
  
  #FOURTH DATA PROCESSING STEP (removes starts/stops of same behavior, at the same time)
  #These are induced if a duration-based orienatation spans an O4 switch, but if they do not change when O4 changes, they can be eliminated here	
  #############################################################	
  #############################################################	
  
  
  
  
  
  
  
  #############################################################################
  ###################################################################
  
  
  #DATA SUMMARIZING STAGE BEGINS HERE
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #Summary measurements extracted from BOPLog output	
  is.odd <- function(x) x %% 2 == 1 #defines is.odd function used to ensure even number of observations for state behaviors
  
  
  
  
  
  #creates empty dataframe that duration behavior summary info goes into
  summarydurationbehaviors<-data.frame(matrix(vector(),0,4,dimnames=list(c(),c("Behavior","Frequency","Duration","Bouts"))),stringsAsFactors=F)
  #Creates variable to be filled with duration behaviors with odd numbers
  odd.behaviors<-c()
  
  #Loops through all duration behaviors identified above and calculates duration, n, and mean time/bout
  for (m in DurationBehaviors){      
    behav.obs <- table[grep(m, table$code), ]
    behav.obs$time<-as.numeric(behav.obs$time)
    if(m=="SP"){
      
      SP1<- subset(behav.obs,code=="SP1")
      SP2<- subset(behav.obs,code=="SP2")   
      SP3<- subset(behav.obs,code=="SP3")  
      SP4<- subset(behav.obs,code=="SP4")
      
      if(nrow(SP1)!=0){
        SP1.SP2.duration<-0  
        for(i in 1:nrow(SP1)){
          SP1.SP2.duration<-SP1.SP2.duration+(SP2$time[i]-SP1$time[i])
        }
        SP1.SP2.freq<-nrow(SP1)
        SP1.SP2.bouts<-  SP1.SP2.duration/ SP1.SP2.freq
        SP1.SP2<-cbind("SP1.SP2",SP1.SP2.freq,SP1.SP2.duration,SP1.SP2.bouts)
      } else {
        SP1.SP2<-data.frame("SP1.SP2",0,0,0)
      }
      colnames(SP1.SP2)<-c("Behavior","Frequency","Duration","Bouts")
      
      if(nrow(SP2)!=0){
        SP2.SP3.duration<-0  
        for(i in 1:nrow(SP2)){
          SP2.SP3.duration<-SP2.SP3.duration+(SP3$time[i]-SP2$time[i])
        }
        SP2.SP3.freq<-nrow(SP2)
        SP2.SP3.bouts<-  SP2.SP3.duration/ SP2.SP3.freq
        SP2.SP3<-cbind("SP2.SP3",SP2.SP3.freq,SP2.SP3.duration,SP2.SP3.bouts)
      } else {
        SP2.SP3<-data.frame("SP2.SP3",0,0,0)
      }
      colnames(SP2.SP3)<-c("Behavior","Frequency","Duration","Bouts")
      
      if(nrow(SP3)!=0){
        SP3.SP4.duration<-0  
        for(i in 1:nrow(SP3)){
          SP3.SP4.duration<-SP3.SP4.duration+(SP4$time[i]-SP3$time[i])
        }
        SP3.SP4.freq<-nrow(SP3)
        SP3.SP4.bouts<-  SP3.SP4.duration/ SP3.SP4.freq
        SP3.SP4<-cbind("SP3.SP4",SP3.SP4.freq,SP3.SP4.duration,SP3.SP4.bouts)
      } else {
        SP3.SP4<-data.frame("SP3.SP4",0,0,0)
      }
      colnames(SP3.SP4)<-c("Behavior","Frequency","Duration","Bouts")
      
    } else {
      if(is.odd(nrow(behav.obs))) {
        behav.obs ##prints table of behaviors if n is odd, suggesting an unclosed duration. 
        odd.behaviors<-c(odd.behaviors,m)
        behav.freq<-"ODD"           #Assigns ODD if behav is mising
        behav.duration<-"ODD"
        behav.bouts<- "ODD" 
      }  else {
        if(nrow(behav.obs)==0) {#Assigns zeros if behav is mising,
          behav.freq<-0
          behav.duration<-0
          behav.bouts<-  0 
        } else {#Otherwise returns durations (behav.duration), frequency (behav.freq), and avg bout length (behav.bouts)
          starts<-behav.obs[grep("start", behav.obs$code), ]
          ends<-behav.obs[grep("end", behav.obs$code), ]
          duration<-0
          for(n in 1:nrow(starts)){
            duration<-duration+(ends$time[n]-starts$time[n])
          }# for every row in 'starts', which has all starts for a behavior, it calculates the difference between end and start and then adds this to 'duration'
          behav.freq<-nrow(starts)
          behav.duration<-duration
          behav.bouts<-  behav.duration/behav.freq
        }
      }
    }
    summary<-cbind(m,behav.freq,behav.duration,behav.bouts)
    colnames(summary)<-c("Behavior","Frequency","Duration","Bouts")
    summarydurationbehaviors<-rbind(summarydurationbehaviors,summary)
    
  }
  summarydurationbehaviors<-rbind(summarydurationbehaviors,SP1.SP2)
  summarydurationbehaviors<-rbind(summarydurationbehaviors,SP2.SP3)
  summarydurationbehaviors<-rbind(summarydurationbehaviors,SP3.SP4)
  
  rownames(summarydurationbehaviors)<-summarydurationbehaviors$Behavior
  summarydurationbehaviors<-summarydurationbehaviors[,-1]
  summarydurationbehaviors$Frequency<-(as.numeric(as.character(summarydurationbehaviors$Frequency)))
  summarydurationbehaviors$Duration<-(as.numeric(as.character(summarydurationbehaviors$Duration)))
  summarydurationbehaviors$Bouts<-(as.numeric(as.character(summarydurationbehaviors$Bouts)))
  summarydurationbehaviors$ids<-rownames(summarydurationbehaviors)
  
  rowwisedurations<-dcast(melt(summarydurationbehaviors, id.var="ids"), 1~ids+variable)
  rowwisedurations<-rowwisedurations[,-1]
  ######################
  #EVENT BEHAVIORS
  summaryeventehaviors<-data.frame(matrix(vector(),0,2,dimnames=list(c(),c("Behavior","Number"))),stringsAsFactors=F)
  for (m in EventBehaviors){       #Loop through each event behavior and counts its occurences
    eventtable<-table[table$code %in% EventBehaviors,]
    #use<-eventtable[jj %in% eventtable$code, ]
    behav.obs <- subset(eventtable,code==m)
    behav.obs$time<-as.numeric(behav.obs$time)
    number<-nrow(behav.obs)
    summary<-cbind(m,number)
    summaryeventehaviors<-rbind(summaryeventehaviors,summary)
  }
  rownames(summaryeventehaviors)<-summaryeventehaviors$m #applies rownames of behaviors
  summaryeventehaviors$number<-as.numeric(as.character(summaryeventehaviors$number)) #
  rowwiseevents<-as.data.frame(t(summaryeventehaviors)) #transpose function makes numbers lose their numeric classification
  rowwiseevents<-rowwiseevents[-1,]
  rowwiseevents<-as.data.frame(t(as.data.frame(sapply(rowwiseevents, function(x) as.numeric(as.character(x))))))  #this function, applied to every column, returns numeric values
  #############################################
  #COMBINED STATE AND EVENT BEHAVIOR SUMMARY DATA
  
  stateANDduration<-cbind(rowwisedurations,rowwiseevents)
  rownames(stateANDduration)<-NULL
  
  #################
  #Additional summary measures
  #####################################
  if (BirdWasOffatBeginning>0){
    accuratetime<-table
    accuratetime$time<-accuratetime$time-TimeBirdFlewON
    TotalTime<-accuratetime[which(accuratetime$code == "END"), ][,1]
    
    table$time<-table$time-TimeBirdFlewON+0.1
  }
  
  
  TotalTimeWatched<-TotalTime-stateANDduration$OFF_Duration
  PropTimeMoving<-stateANDduration$BP1_Duration/TotalTimeWatched
  
  Nfemales<-length(table[which(table$code == "Fpres"),][,1])
  Nfemales.disp<-length(table[which(table$code == "Fdisp"),][,1])
  Nmales<-length(table[which(table$code == "Mpres"),][,1])
  
  Time.upright <- TotalTimeWatched-stateANDduration$O3_Duration
  Time.inverted<- stateANDduration$O3_Duration
  
  
  O1.obs <- subset(table,code=="O1")
  O2.obs <- subset(table,code=="O2")
  timespent<-0
  femaleorientedbouts<-nrow(O1.obs)#Number of bouts oriented towards females
  if(nrow(O1.obs)!=0){ #IF there are is at least 1 obs of O1, then we can calculate time oriented towards female
    
    
    if(nrow(O2.obs)!=0){ #IF there are O1s and O2s, we can use O2s to calculate end of O1s
      
      
      for(j in 1:femaleorientedbouts){
        lookingstart<-as.numeric(as.character(O1.obs[j,"time"]))
        alltimesafterorientingtowardsfemale<-subset(O2.obs,time>lookingstart)
        lookingend<-0
        
        if(nrow(alltimesafterorientingtowardsfemale)>0){
          lookingend<-as.numeric(as.character(min(alltimesafterorientingtowardsfemale$time,na.rm=TRUE)))
        } else {	lookingend<-TotalTime
        }
        
        timespent<-timespent+(lookingend-lookingstart)
      }
    } else {#If there is an O1 and no subsequent O2s, then we use the end time (TotalTime) to calc duration spent oriented towards female
      timespent<-TotalTime-as.numeric(as.character(O1.obs[1,"time"]))
    }
  }
  
  Time.orientedtofemale <- timespent
  
  
  lengthodds<-length(odd.behaviors)
  
  
  odd.behaviors<-paste(odd.behaviors,collapse="-")
  
  behaviorsscored<-nrow(table)
  
  belongsto<-species.name
  
  unique.behaviors.scored<-as.data.frame(unique(table$code))
  colnames(unique.behaviors.scored)<-"uniques"
  n.unique.behavs<-length(unique.behaviors.scored[ grep("end", unique.behaviors.scored$uniques, invert = TRUE,ignore.case = TRUE) , ])
  
  
  #  printme<-data.frame(FileisNamed,belongsto,TotalTimeWatched,behaviorsscored,n.unique.behavs,Nmales,Nfemales,Nfemales.disp,Time.upright,Time.inverted,Time.orientedtofemale,femaleorientedbouts,stateANDduration,odd.behaviors,off.ss.eq,happenedwhileoff)
  
  origsourcefolder<-getwd()
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  
  
  compositebehavioralmetrictable<-table[table$code!="END" & table$code!="Mpres" & table$code!="Fpres"&
                                          # table$code!="OFF(start)"& table$code!="OFF(end)" &
                                          table$code!="Fdisp",]
  
  compositebehavioralmetrictable$time<-round(compositebehavioralmetrictable$time,digit=1)
  
  #adds 0.1 sec to OFF(starts) to avoid overlap with behaviors that end as birds goes OFF
  compositebehavioralmetrictable[compositebehavioralmetrictable$code=="OFF(start)","time"]<-compositebehavioralmetrictable[compositebehavioralmetrictable$code=="OFF(start)","time"]+0.1
  
  #subtracts 0.1 sec to OFF(end)s to avoid overlap with behaviors that start as bird comes back on
  compositebehavioralmetrictable[compositebehavioralmetrictable$code=="OFF(end)","time"]<-compositebehavioralmetrictable[compositebehavioralmetrictable$code=="OFF(end)","time"]-0.1
  
  newlength<-round(TotalTime/0.1)+2 #Complete # of 0.1 sec intervals for the video (adds 1, b/c we add 0.1 to the first behavior so it didn't occur at 0 = 0th column)
  #  allbehaviorstoconverttenthsecondintervals<-cbind(StartList2,EndList)
  
  
  AllBehaviors<-c("O5","O6","O3","O4",
                  "BP1","BP2","BP3",
                  "SS1","SS2",
                  "OFF", "O2",
                  
                  #"O1","O2",
                  "OPMH","OPMB","OPMF","OPMW","OPMT",
                  "OPAH","OPAC1","OPAC2","OPAC3","OPAW","OPATW","OPAMW","OPAMTT",
                  "OPABH","OPABB","OPABF","OPABW","OPABT",
                  "PC1","PC2","PC3","PC4",
                  #"RepWing","RepTail","RepHead","RepTorso",
                  "RM1","RM2","RM3","RM4",
                  "MO1","MO2","PU1",
                  "SP1","SP2","SP3","SP4")
  
  
  
  allbehaviorstolatercombine<-length(AllBehaviors)
  completebehaviormatrix<-matrix(data=NA,nrow=allbehaviorstolatercombine,ncol=newlength)
  row.names(completebehaviormatrix)<-AllBehaviors
  
  DurationBehaviors.nosp.no078<-c("O5","O6","O3","O4",
                                  "BP1","BP3","OFF","SS1","SS2",
                                  "OPMH.D","OPMB.D","OPMF.D","OPMW.D","OPMT.D",
                                  #"RepWing","RepTail","RepHead","RepTorso",
                                  #"RM1","RM2","RM3","RM4",
                                  "MO2","PU1")
  
  
  
  #Special class (0,1,2,3) determines which specials are treated as duration behaviors  
  if(specialclass==0){
    DurationBehaviors.toUse<-DurationBehaviors.nosp.no078
  }
  if(specialclass==1){
    DurationBehaviors.toUse<-c(DurationBehaviors.nosp.no078,"SP3")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP3"]<-"SP3(start)")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP4"]<-"SP3(end)")
  }
  if(specialclass==2){
    DurationBehaviors.toUse<-c(DurationBehaviors.nosp.no078,"SP1")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP1"]<-"SP1(start)")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP2"]<-"SP1(end)")
  }
  if(specialclass==3){
    DurationBehaviors.toUse<-c(DurationBehaviors.nosp.no078,c("SP1","SP3"))
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP1"]<-"SP1(start)")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP2"]<-"SP1(end)")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP3"]<-"SP3(start)")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP4"]<-"SP3(end)")
    
    
  }
  
  for (cb in 1:length(DurationBehaviors.toUse)){ #loop through all initially scored behaviors, but not Specials
    
    thisone<-DurationBehaviors.toUse[cb] #gives new name, "thisone" to the current behavior
    
    start.code<-paste(thisone,"(start)",sep="")
    end.code<-paste(thisone,"(end)",sep="")
    
    starts<-compositebehavioralmetrictable[compositebehavioralmetrictable$code==start.code, ] #df of starts
    ends<-compositebehavioralmetrictable[compositebehavioralmetrictable$code==end.code, ] #df of ends
    
    number.of.starts<-length(starts[,2])
    
    if(number.of.starts!=0){  #Only run the "fill-in loop" if the behavior was recorded at all    
      for (ooo in 1:number.of.starts)  { #This loop fills the matrix with 1s for every 1/10 sec where the behavior is occuring
        
        columnstofill<-10* (seq(starts[ooo,1],ends[ooo,1],by=0.1)   )
        
        if(grepl(".D",thisone)){
          thisone<-substr(thisone,1,nchar(thisone)-2)
        } 
        
        completebehaviormatrix[rownames(completebehaviormatrix)==thisone,columnstofill]<-1 #fills columns while behavior is active
        
        
      }
    }
  }  
  
  
  
  EventBehaviors.Fill<-c("BP2","O1","O2",
                         "OPMH","OPMB","OPMF","OPMW","OPMT","OPAH",
                         "OPAC1","OPAC2","OPAC3","OPAW","OPATW","OPAMW","OPAMTT",
                         "OPABH","OPABB","OPABF","OPABW","OPABT",
                         "PC1","PC2","PC3","PC4",
                         "MO1","SP1","SP2","SP3","SP4")
  
  event.length<-length(EventBehaviors)
  
  for (ebs in 1:event.length){
    
    behaviortopull<-EventBehaviors.Fill[ebs]
    behavs.to.fill<-compositebehavioralmetrictable[compositebehavioralmetrictable$code==behaviortopull, ]
    column.to.add<-10*behavs.to.fill$time
    column.to.add<-ifelse(column.to.add==0,1,column.to.add)
    
    completebehaviormatrix[rownames(completebehaviormatrix)==behaviortopull,column.to.add]<-1 #fills columns while behavior is active
    
  }
  
  
  
  compositebehaviorcategorizer<-function(dataframe,x){
    x2<-dataframe[,x,drop=FALSE]
    yyy<-as.data.frame(x2[complete.cases(x2),])
    behav<-paste(rownames(yyy),sep=".",collapse = ".")
    
  }  
  
  temporalbehaviorcategories<-matrix(data=NA,nrow=1,ncol=newlength)
  
  for(pdq in 1:newlength){
    temporalbehaviorcategories[1,pdq]<-compositebehaviorcategorizer(completebehaviormatrix,pdq)
  }
  
  
  
  as.data.frame(unique(temporalbehaviorcategories[1,]))
  
  
  no.OFFs<-as.matrix(temporalbehaviorcategories[,!grepl("OFF",temporalbehaviorcategories[1,])])#gets rid of behavior combos with OFF
  no.OFFs[no.OFFs == ""] <- NA
  no.OFFs<-na.omit(no.OFFs)
  
  
  
  ###########################################################################################################
  #Finds behaviors missing a body-axis (O5/O6) and adds axis from nearest (time-wise) axis
  deficient.behaviors<-no.OFFs[((!grepl("O5", no.OFFs[,1], fixed=TRUE)) & (!grepl("O6", no.OFFs[,1], fixed=TRUE)) ),]
  containing.behaviors<-no.OFFs[((grepl("O5", no.OFFs[,1], fixed=TRUE)) | (grepl("O6", no.OFFs[,1], fixed=TRUE)) ),]
  
  index.of.deficient<-which((!grepl("O5", no.OFFs[,1], fixed=TRUE)) & (!grepl("O6", no.OFFs[,1], fixed=TRUE)) )
  index.of.containing<-which((grepl("O5", no.OFFs[,1], fixed=TRUE)) | (grepl("O6", no.OFFs[,1], fixed=TRUE)) )
  
  ####
  #Little loop through all deficient (posture-absent) behaviors, adding posture from nearest time
  if(length(index.of.deficient)>0){
    for(defic in 1:length(index.of.deficient)){
      defbehav<-deficient.behaviors[defic]
      defind<-index.of.deficient[defic]
      
      replacementbehavior<-no.OFFs[which.min(abs(index.of.containing-defind)),]
      replacementposture<-ifelse(grepl("O5",replacementbehavior),"O5","O6")
      #replacementpostureDOT<-paste(replacementposture,".", sep ="")
      
      
      
      newbehaviorvalue<- paste(replacementposture,defbehav,sep=".")
      
      
      
      
      no.OFFs[defind,1]<-newbehaviorvalue
      
    }
  }  
  ####################################################
  #######################################################
  #no.OFFs<-as.matrix(gsub("O5.O6","O5",no.OFFs[,1]))#swaps out impossible O5/O6
  
  #need to deal with O5.O6s (maybe randomly choose one?)
  no.OFFs[((grepl("O5.O6", no.OFFs[,1], fixed=TRUE)) ),]
  simultaneousO5O6index<-grep("O5.O6", no.OFFs[,1], fixed=TRUE,value=FALSE)
  fixerO5O6er<-c("O5","O6")
  
  if(length(simultaneousO5O6index)>0){
    for(abc in 1:length(simultaneousO5O6index)){
      no.OFFs[simultaneousO5O6index[abc],1]<-as.matrix(gsub("O5.O6",fixerO5O6er[sample(1:2,1)],no.OFFs[simultaneousO5O6index[abc],1]))
    }
  } #randomly replaces any simultatneous O5.O6 with either O5 or O6
  
  
  
  #no.OFFs<-as.matrix(gsub("O5.O6",fixerO5O6er[sample(1:2,1)],no.OFFs[,1]))#swaps out impossible O5/O6
  
  
  #Counts occurences for each ID type
  summarytable<-table(no.OFFs)
  behaviorcounts<-table(no.OFFs)
  proportiontable<-prop.table(summarytable)
  uniquebehaviors<-length(summarytable)
  
  six<-no.OFFs
  #ALL BEHAVIORS and TRANSITIONS
  if(uniquebehaviors>1){
    
    ##This loop creates the denominator for the calculation of the inverse Simpson Index for behavioral diversity
    val2<-0
    for (u in 1:uniquebehaviors){
      proptosquare<-proportiontable[u]
      val <-proptosquare^2.0
      val2<-val+val2
    }			
    #inverse Simpson Index for behavioral diversity
    iSimpson<-1.0/as.numeric(val2)
    relativebehavioraldiversity<-iSimpson/(uniquebehaviors) #relative color diversity, accounting for # of color classes
    
    ##########################################################################################################
    #uses "createSequenceMatrix" function (from *markovchain*) to calculate transition matrix
    c.no.OFFs<-as.character(no.OFFs)
    
    
    TransitionMatrix.markov<-createSequenceMatrix(c.no.OFFs,toRowProbs = TRUE,sanitize=FALSE)
    TransitionMatrix.cum<-createSequenceMatrix(c.no.OFFs,sanitize=FALSE)
    
    #TransitionMatrix<-TransitionMatrix[-ncol(TransitionMatrix),-ncol(TransitionMatrix)]
    
    
    time.value<-sum(TransitionMatrix.cum)/10 #seconds watched
    
    
    noselfs<-TransitionMatrix.cum
    diag(noselfs)<-NA
    
    prop.table.excludeNAs<- function(x) {x/sum(x, na.rm=TRUE)}
    proportiontable.trans.Simpson<-prop.table.excludeNAs(noselfs) #get rid of diagonal/self-transitions for Simpson Analyses
    
    proportiontable.trans<-prop.table.excludeNAs(TransitionMatrix.cum) 
    
    
    if(uniquebehaviors>1){
      val3<-matrix(nrow=(((uniquebehaviors^2)-uniquebehaviors)/2),ncol=2)
      rownum<-1
      for (u in 1:uniquebehaviors){
        for(v in u:uniquebehaviors){#iteratively reduces columns analyzed in next loop to avoid double counting transitions 
          if (u!=v){
            val3[rownum,2]<-as.numeric(proportiontable.trans.Simpson[u,v]+proportiontable.trans.Simpson[v,u])
            val3[rownum,1]<-paste(row.names(proportiontable.trans.Simpson)[u],colnames(proportiontable.trans.Simpson)[v],sep="-")
            rownum<-rownum+1
          }
        }
      }
      val4<-as.data.frame(val3)
      val4[,2]<-as.numeric(as.character(val4[,2])) 
      
      simpsonT<-0
      
      for(z in 1:nrow(val4)){
        propsq<-val4[z,2]^2.0
        simpsonT<-simpsonT+propsq
      }
      
      
      iSimpsonT<-1.0/simpsonT #Calculates inverse Simpson index for horizontal transitions
      relativetransitiondiversity<-iSimpsonT/(uniquebehaviors*(uniquebehaviors-1)/2)  
    } else {
      
      iSimpsonT <- NA
      relativetransitiondiversity<- NA
    }
    
    
    PropSelfTransitioning<-sum(diag(TransitionMatrix.cum))/sum(TransitionMatrix.cum)
    
    network.TransitionMatrix.cum<-network(proportiontable.trans,directed=TRUE,loops=TRUE,hyper=FALSE,matrix.type="adjacency")
    
    
    
    showme<-graph_from_adjacency_matrix(proportiontable.trans,mode="directed",weighted=TRUE,diag=FALSE,add.colnames=TRUE)
    showme.undirected<-graph_from_adjacency_matrix(proportiontable.trans,mode="undirected",weighted=TRUE,diag=FALSE,add.colnames=TRUE)
    
    showme.looped<- graph_from_adjacency_matrix(proportiontable.trans,mode="directed",weighted=TRUE,diag=TRUE,add.colnames=TRUE)
    showme.looped.undirected<- graph_from_adjacency_matrix(proportiontable.trans,mode="undirected",weighted=TRUE,diag=TRUE,add.colnames=TRUE)
    
    
    
    E(showme)$width <- 1+E(showme)$weight/50
    V(showme)$nombres<-colnames(proportiontable.trans)
    
    E(showme.looped)$width <- 1+E(showme.looped)$weight/50
    V(showme.looped)$nombres<-colnames(proportiontable.trans)
    
    # Compute approx self-transition and use that to set node size:
    selfsizesscaled<-(diag(proportiontable.trans))*100+1
    V(showme)$selfsize<-sqrt(selfsizesscaled)
    V(showme)$size <- ((selfsizesscaled)^2)*9
    
    V(showme.looped)$selfsize<-log(selfsizesscaled)+1
    
    
    
    #####################################################
    # Network Clustering
    # l <- layout.fruchterman.reingold(showme.looped)
    # l <- layout_with_kk(showme.looped)
    # l <- layout_with_lgl(showme.looped)
    l <- layout_with_graphopt(showme.looped)
    l2 <-layout_with_graphopt(showme)
    #l <- layout_nicely(showme.looped)
    #behaviorplotlist
    
    
    #clp <- cluster_label_prop(showme.looped)
    clp <-cluster_label_prop(showme.looped.undirected) # input graph should be undirected to make sense.
    clp2<-cluster_label_prop(showme.undirected)
    # setwd("C:/Users/Rusty/Amazon Drive/BOP/Video_analysis")
    # source("PictoGrammerMega_Dec13ab.R", chdir = F)
    # #PictoGrammer<-function(proportiontable,hangingbird)
    #
    
    
    
    
    ################################################################################################################################
    #PictoGrammer(proportiontable.trans,hangingbird,uniquebehaviors)
    
    origwd<-dir
    setwd(dir<-"C:/Users/Rusty/Amazon Drive/BOP/Video_analysis/Pictograms")
    
    pictogram.names<-list.files(dir,full.names=FALSE,pattern=".png")
    
    proportiontable<-proportiontable
    hangingbird<-hangingbird
    
    
    selflooppoints<-read.table("SelfLoop.csv",sep=",",header=F)
    
    #Creates empty list, then loops through pictogram directory, pulling in each PG, adding it to the list, and naming it
    piclist<-list()
    for (pg in 1:56){
      img <- readPNG(pictogram.names[pg])  
      
      if(pg>23 && pg <42) {
        w <- matrix(rgb(img[,,1],img[,,2],img[,,3], img[,,4] * 0.75), nrow=dim(img)[1]) #0.5 is alpha
        rimg <- as.raster(w) # raster multilayer object
      } else {
        rimg <- as.raster(img) # raster multilayer object  
      }
      pgname<-gsub('.{4}$', '',pictogram.names[pg])
      piclist[[pg]]<-rimg
      names(piclist)[[pg]]<-pgname
    }
    
    #Loops through each unique behavior (from proportiontable) and creates a new, unique, behavioral pictogram
    behaviorplotlist<-list()
    par(bg="transparent")
    plot.new()
    goo=1
    #pdf("Behav30-trans-fix.pdf",width= 12, height= 12,family="NimbusRom")
    #op <- par(mfrow=c(6,6))# rows, columns
    #par(mar=c(1,1,1,1)) 
    
    
    colorweights<-E(showme.looped)$weight/min(E(showme.looped)$weight) 
    colorweights<-log(colorweights+0.01)+1
    
    Lab.palette<-colorRampPalette(c("blue",  "red"),space = "Lab")
    #  Lab.palette<-colorRampPalette(c("#ffeda0","#feb24c","#f03b20"),space = "Lab")
    
    
    E(showme.looped)$color<-Lab.palette(max(colorweights))[colorweights]
    
    col.tr <- grDevices::adjustcolor(E(showme.looped)$color, alpha=1) #adds translucency to arrows
    
    SELFCOLORS<-log(diag(proportiontable.trans)/min(E(showme.looped)$weight))+1 
    SELFCOLORS.colors<-SELFCOLORS
    #Little loop that adds same heatmap scaled colors to self-transitions, for use in behavior boxes
    for(ggg in 1:length(SELFCOLORS)){
      selfvalue<-SELFCOLORS[ggg]
      if(selfvalue>0){
        SELFCOLORS.colors[ggg]<-Lab.palette(max(colorweights))[selfvalue]
      } else {
        SELFCOLORS.colors[ggg]<-NA
      }
      
    }
    
    
    for (ub in 1:uniquebehaviors){
      par(bg="transparent")
      graphics::plot(goo, type ="l", xlab="", ylab="", xlim=c(0, 100), ylim=c(0, 100),axes=FALSE)
      
      #53
      if( (length(grep("SS1", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`SS1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #54
      if( (length(grep("SS2", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`SS2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      ##########################################################  
      #2
      if( (length(grep("BP1", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`BP1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #1
      if( (length(grep("BP1.BP2", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`BP1-BP2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #3
      if( (length(grep("BP3", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`BP3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #4
      if( (length(grep("O1", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`O1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #5
      if( (length(grep("O2", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      ###############################################
      
      #If hangingbird  
      if(hangingbird==1){  
        #6
        if( (length(grep("O5.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1){
          rasterImage(piclist$`O5-O3-O4-MO1or2a`, xleft = 1, ybottom=1, xright=100,ytop=90)  
        }
        #7
        if( (length(grep("O5.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
          rasterImage(piclist$`O5-O3-O4a`, xleft = 1, ybottom=1, xright=100,ytop=90)  
        }
        
        #10
        if( (length(grep("O6.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1){
          rasterImage(piclist$`O6-O3-O4-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
        }
        #11
        if( (length(grep("O6.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
          rasterImage(piclist$`O6-O3-O4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
        }  
      }
      
      #If non-hangingbird
      if(hangingbird==0){  
        #8
        if( (length(grep("O5.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1){
          rasterImage(piclist$`O5-O3-O4-MO1or2b`, xleft = 1, ybottom=1, xright=100,ytop=90)  
        }
        #9
        if( (length(grep("O5.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
          rasterImage(piclist$`O5-O3-O4b`, xleft = 1, ybottom=1, xright=100,ytop=90)  
        }
        
        #10.a
        if( (length(grep("O6.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1){
          rasterImage(piclist$`O6-O3-O4-MO1or2a`, xleft = 1, ybottom=1, xright=100,ytop=90)  
        }
        #11.a
        if( (length(grep("O6.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
          rasterImage(piclist$`O6-O3-O4a`, xleft = 1, ybottom=1, xright=100,ytop=90)  
        }  
      }
      
      #12
      if( (length(grep("O5.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1 && (length(grep("O3.O4.O5", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O5-O4-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #13
      if( (length(grep("O5.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0 && (length(grep("O3.O4.O5", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O5-O4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }  
      #14
      if( (length(grep("O6.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1 && (length(grep("O3.O4.O6", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O6-O4-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #15
      if( (length(grep("O6.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0 && (length(grep("O3.O4.O6", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O6-O4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      
      #18
      if( (length(grep("O5", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1 && (length(grep("O4", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O5-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #19
      if( (length(grep("O5", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0 && (length(grep("O4", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O5`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #20
      if( (length(grep("O6", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1 && (length(grep("O4", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O6-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #21
      if( (length(grep("O6", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0 && (length(grep("O4", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O6`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      ###############################################  
      #22
      if( (length(grep("OPABB", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPABB`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }  
      #23
      if( (length(grep("OPABF", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPABF`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      } 
      #24
      if( (length(grep("OPABH", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPABH`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #25
      if( (length(grep("OPABT", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPABT`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #26
      if( (length(grep("OPABW", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPABW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #27
      if( (length(grep("OPAC1", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPAC1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #28
      if( (length(grep("OPAC2", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPAC2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #29
      if( (length(grep("OPAC3", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPAC3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #30
      if( (length(grep("OPAH", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPAH`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #31
      if( (length(grep("OPAMTT", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPAMTT`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #32
      if( (length(grep("OPAMW", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPAMW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #33
      if( (length(grep("OPATW", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPATW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #34
      if( (length(grep("OPAW", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPAW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #35
      if( (length(grep("OPMB", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPMB`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #36
      if( (length(grep("OPMF", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPMF`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #37
      if( (length(grep("OPMH", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPMH`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #38
      if( (length(grep("OPMT", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPMT`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #39
      if( (length(grep("OPMW", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPMW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      ######################################
      #40
      if( (length(grep("PC1", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`PC1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #41
      if( (length(grep("PC2", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`PC2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #42
      if( (length(grep("PC3", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`PC3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #43
      if( (length(grep("PC4", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`PC4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #44
      if( (length(grep("PU1", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`PU1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      ######################################  
      #45
      if( (length(grep("RM1", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`RM1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #46
      if( (length(grep("RM2", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`RM2`, xleft = 1, ybottom=1, xright=100,ytop=90)    
      }
      #47
      if( (length(grep("RM3", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`RM3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #48
      if( (length(grep("RM4", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`RM4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      ######################################  
      #49
      if( (length(grep("SP1", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`SP1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #50
      if( (length(grep("SP2", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`SP2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }  
      #51
      if( (length(grep("SP3", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`SP3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      } 
      #52
      if( (length(grep("SP4", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`SP4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }  
      
      ############  
      text(50,88,ub,cex=1.5) #ADD BEHAVIOR ID NUMBER
      ########################################### 
      ############  
      # box(which="plot",lty="solid",lwd=30,col=SELFCOLORS.colors[ub]) #ADD HEATMAPPED COLORED BOX AROUND PLOT, CORRESPONDING TO LOOPS
      arrow.head.color<-ifelse(is.na(SELFCOLORS.colors[ub]),NA,SELFCOLORS.colors[ub])
      
      #   diagram::curvedarrow(from=c(45,80), to=c(55,80), lwd = 7, lty = 1, 
      #                        # lcol = "black", arr.col = "black", 
      #                        lcol = SELFCOLORS.colors[ub], arr.col = NA,
      #                        arr.pos = 1, curve = -2, dr = 0.01, arr.type="triangle",endhead=TRUE,
      #                        segment = c(0, 1))
      #   
      # polygon(x=c(58,52,55),y=c(86,86,75),col=arrow.head.color,border=arrow.head.color)
      if(!(is.na(arrow.head.color))) {
        polygon(x=30+0.4*selflooppoints[,1],y=79+0.25*selflooppoints[,2],col=arrow.head.color)
      }
      ########################################### 
      setwd("C:/Users/Rusty/Amazon Drive/BOP/Video_analysis/Pictograms/temp")
      pngplot <- recordPlot() 
      xyz<-names(proportiontable)[ub]
      individualbehaviorname<-paste(xyz,"png",sep=".")
      
      png(file = individualbehaviorname, bg = "transparent")
      replayPlot(pngplot)
      dev.off()
      
      completepicture<-readPNG(individualbehaviorname)
      rimg <- as.raster(completepicture) # raster multilayer object
      
      newpgname<-gsub('.{4}$', '',xyz)
      behaviorplotlist[[ub]]<-rimg
      names(behaviorplotlist)[[ub]]<-xyz
      
      setwd(dir)
    }
    
    #dev.off()
    
    
    
    #  behaviorplotlist is a named list corresponding to all unique behaviors from a given trial
    
    
    
    setwd(origwd)
    
    
    
    
    
    
    
    
    
    # NETWORK SUMMARY VARIABLES
    
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    #6.1 Density
    
    #The proportion of present edges from all possible edges in the network.
    
    #densityscore<-ecount(showme.looped)/(vcount(showme.looped)*(vcount(showme.looped)-1)) #for a directed network
    densityscore<-ecount(showme.looped)/(vcount(showme.looped)^2)#for a directed network WITH self-loops
    
    #Mean degree
    deg <- degree(showme.looped,loops=FALSE)
    mean.degree<-mean(deg)
    
    #Mean path length
    mean_path_length<-mean_distance(showme.looped, directed=T)
    
    
    #network diameter is the longest geodesic distance (length of the shortest path between two nodes) in the network.
    network.diameter<-diameter(showme.looped, directed=T,weights=NA)
    
    
    #average clustering coefficient
    avg.cluster.coeff<-transitivity(showme.looped)
    
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
    
    printme<-data.frame(FileisNamed,belongsto,TotalTimeWatched,behaviorsscored,n.unique.behavs,
                        uniquebehaviors,densityscore,mean.degree,mean_path_length,network.diameter,avg.cluster.coeff,
                        iSimpson,relativebehavioraldiversity,iSimpsonT,relativetransitiondiversity,
                        PropSelfTransitioning,Smallworldness,
                        Nmales,Nfemales,Nfemales.disp,Time.upright,Time.inverted,Time.orientedtofemale,
                        femaleorientedbouts,stateANDduration,odd.behaviors,off.ss.eq,happenedwhileoff)
    
    
    
    
    
    
    #CONNECTIVITY measures reflect degree to which network differs from complete network
    #Edge density: % of edges compared to maximum
    #Average degree: Avg. number of links
    #Average path length: avg of shortest pathes between reachable nodes
    #Network diameter: longest of shortest paths
    
    #CENTRALITY measures quantify heterogeneity in network structure
    #Average clustering coefficient:
    #Components
    
    
    
    
    
    
    
    
    
    
    
    #V(showme.looped)$size<-80
    netname<-paste(name2,"Network_cust.pdf",sep='_')
    
    setwd(location)
    
    pdf(netname,width= dimensions[1], height= dimensions[2],family="NimbusRom")
    op <- par(mfrow=c(1,1))# rows, columns
    par(mar=c(0,0,0,0)) 
    # plot(clp, showme.looped,edge.curved=.1,edge.arrow.size = 0.45,vertex.size=18,layout=l,rescale=FALSE,
    #      xlim=range(l[,1]),ylim=range(l[,2]))
    
    
    # plot(clp, showme.looped,edge.curved=.1,edge.arrow.size = 0.45,vertex.size=18,layout=l,rescale=TRUE,
    #      xlim=c(-1,1),ylim=c(-1,1))
    
    V(showme.looped)$raster <- behaviorplotlist
    V(showme)$raster <- behaviorplotlist
    
    
    
    
    # print(  plot(clp, showme.looped,edge.curved=.1,layout=l,rescale=TRUE,vertex.shape="raster",
    #        vertex.label=NA, 
    #        #vertex.size=15*(V(showme.looped)$selfsize), # resizes vertices based on self-self transitions
    #        vertex.size=25,
    #        edge.color=col.tr,
    #        edge.arrow.size=0.35,
    #        
    #        edge.width=3,
    #        edge.arrow.width=1.6,
    #        edge.loop.angle=1.5,
    #        
    #        xlim=c(-1,1),ylim=c(-1,1)))
    
    print(  plot(clp, showme,edge.curved=.1,layout=l,rescale=TRUE,vertex.shape="raster",
                 vertex.label=NA, 
                 #vertex.size=15*(V(showme.looped)$selfsize), # resizes vertices based on self-self transitions
                 vertex.size=25,
                 edge.color=col.tr,
                 edge.arrow.size=0.35,
                 
                 edge.width=3,
                 edge.arrow.width=1.6,
                 edge.loop.angle=1.5,
                 
                 xlim=c(-2,2),ylim=c(-1,1)))
    
    print(text(0.9,-0.75,paste("Timescored(s) = ",time.value,sep=""),col="black",cex=0.75,adj=(0)))
    print(text(0.9,-0.8,paste("behaviors = ",uniquebehaviors,sep=""),col="black",cex=0.75,adj=(0)))  
    print(text(0.9,-0.85,paste("mean deg = ",round(mean.degree,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
    print(text(0.9,-0.9,paste("density score = ",round(densityscore,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
    print(text(0.9,-0.95,paste("mean path L = ",round(mean_path_length,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
    print(text(0.9,-1,paste("network.diameter = ",network.diameter,sep=""),col="black",cex=0.75,adj=(0)))  
    print(text(0.9,-1.05,paste("Smallworldness = ",Smallworldness,sep=""),col="black",cex=0.75,adj=(0)))
    
    
    #tkplot(clp, showme.looped)
    #####  par("usr")
    
    # library(plotrix)
    # 
    # scalelist<-list()
    # gettinsmaller<-seq(from=1,to=0.52, by=-0.02)
    # for (vc in 1:25){
    # 
    #   
    #   gettingsmaller<-gettinsmaller[vc]
    #   pointlocations<-cbind(rescale(l[,1],gettingsmaller*c(-1,1)),rescale(l[,2],gettingsmaller*c(-1,1)))
    #   par(new=TRUE)
    #   plot(pointlocations,xlab="", ylab="",xlim=c(-1,1),ylim=c(-1,1),pch=vc)
    # 
    # }
    # dev.off()
    
    
    
    # for(bbb in 1:uniquebehaviors){
    #   rasterImage(behaviorplotlist[[bbb]], xleft = pointlocations[bbb,1], xright=pointlocations[bbb,1]+0.15,
    #               ybottom=pointlocations[bbb,2], ytop=pointlocations[bbb,2]+0.15)  
    #   par(new=TRUE)
    #   plot(pointlocations,xlab="", ylab="",xlim=c(-1,1),ylim=c(-1,1))
    # }
    # 
    # pointlocations<-cbind(rescale(l[,1],0.92*c(-1,1)),rescale(l[,2],0.92*c(-1,1)))
    # pointlocations<-pointlocations
    # 
    # for(bbb in 1:uniquebehaviors){
    #   rasterImage(behaviorplotlist[[bbb]], xleft = pointlocations[bbb,1], xright=pointlocations[bbb,1]+0.15,
    #               ybottom=pointlocations[bbb,2], ytop=pointlocations[bbb,2]+0.15)  
    #   par(new=TRUE)
    #   
    # }
    
    dev.off()
    
    #Another network program
    #geph
    
  }
  
  #FILTERED BEHAVIORS (minus behaviors that only happened once)
  no.OFFs.filtered1<-no.OFFs[which(duplicated(no.OFFs)),,drop=FALSE]
  
  summarytable.filtered1<-table(no.OFFs.filtered1)
  behaviorcounts.filtered1<-table(no.OFFs.filtered1)
  proportiontable.filtered1<-prop.table(summarytable.filtered1)
  uniquebehaviors.filtered1<-length(summarytable.filtered1)  
  
  seven<-no.OFFs.filtered1
  
  if(uniquebehaviors.filtered1>1){  
    ##This loop creates the denominator for the calculation of the inverse Simpson Index for behavioral diversity
    val2.filtered1<-0
    for (u in 1:uniquebehaviors.filtered1){
      proptosquare<-proportiontable[u]
      val <-proptosquare^2.0
      val2.filtered1<-val+val2.filtered1
    }			
    #inverse Simpson Index for behavioral diversity
    iSimpson.filtered1<-1.0/as.numeric(val2.filtered1)
    relativebehavioraldiversity.filtered1<-iSimpson.filtered1/(uniquebehaviors.filtered1) #relative color diversity, accounting for # of color classes
    
    ##########################################################################################################
    #uses "createSequenceMatrix" function (from *markovchain*) to calculate transition matrix
    c.no.OFFs<-as.character(no.OFFs.filtered1)
    
    
    TransitionMatrix.markov<-createSequenceMatrix(c.no.OFFs,toRowProbs = TRUE,sanitize=FALSE)
    TransitionMatrix.cum.filtered1<-createSequenceMatrix(c.no.OFFs,sanitize=FALSE)
    
    #TransitionMatrix<-TransitionMatrix[-ncol(TransitionMatrix),-ncol(TransitionMatrix)]
    
    time.value<-sum(TransitionMatrix.cum.filtered1)/10 #seconds watched
    
    
    
    noselfs<-TransitionMatrix.cum.filtered1
    diag(noselfs)<-NA
    
    prop.table.excludeNAs<- function(x) {x/sum(x, na.rm=TRUE)}
    proportiontable.trans.Simpson<-prop.table.excludeNAs(noselfs) #get rid of diagonal/self-transitions for Simpson Analyses
    
    proportiontable.trans<-prop.table.excludeNAs(TransitionMatrix.cum.filtered1) 
    
    
    if(uniquebehaviors.filtered1>1){
      val3<-matrix(nrow=(((uniquebehaviors.filtered1^2)-uniquebehaviors.filtered1)/2),ncol=2)
      rownum<-1
      for (u in 1:uniquebehaviors.filtered1){
        for(v in u:uniquebehaviors.filtered1){#iteratively reduces columns analyzed in next loop to avoid double counting transitions 
          if (u!=v){
            val3[rownum,2]<-as.numeric(proportiontable.trans.Simpson[u,v]+proportiontable.trans.Simpson[v,u])
            val3[rownum,1]<-paste(row.names(proportiontable.trans.Simpson)[u],colnames(proportiontable.trans.Simpson)[v],sep="-")
            rownum<-rownum+1
          }
        }
      }
      val4<-as.data.frame(val3)
      val4[,2]<-as.numeric(as.character(val4[,2])) 
      
      simpsonT<-0
      
      for(z in 1:nrow(val4)){
        propsq<-val4[z,2]^2.0
        simpsonT<-simpsonT+propsq
      }
      
      
      iSimpsonT.filtered1<-1.0/simpsonT #Calculates inverse Simpson index for horizontal transitions
      relativetransitiondiversity.filtered1<-iSimpsonT.filtered1/(uniquebehaviors.filtered1*(uniquebehaviors.filtered1-1)/2)  
    } else {
      
      iSimpsonT.filtered1 <- NA
      relativetransitiondiversity.filtered1<- NA
    }
    
    
    PropSelfTransitioning.filtered1<-sum(diag(TransitionMatrix.cum.filtered1))/sum(TransitionMatrix.cum.filtered1)
    
    network.TransitionMatrix.cum.filtered1<-network(proportiontable.trans,directed=TRUE,loops=TRUE,hyper=FALSE,matrix.type="adjacency")
    
    
    
    showme<-graph_from_adjacency_matrix(proportiontable.trans,mode="directed",weighted=TRUE,diag=FALSE,add.colnames=TRUE)
    showme.undirected<-graph_from_adjacency_matrix(proportiontable.trans,mode="undirected",weighted=TRUE,diag=FALSE,add.colnames=TRUE)
    
    showme.looped<- graph_from_adjacency_matrix(proportiontable.trans,mode="directed",weighted=TRUE,diag=TRUE,add.colnames=TRUE)
    showme.looped.undirected<- graph_from_adjacency_matrix(proportiontable.trans,mode="undirected",weighted=TRUE,diag=TRUE,add.colnames=TRUE)
    
    
    
    E(showme)$width <- 1+E(showme)$weight/50
    V(showme)$nombres<-colnames(proportiontable.trans)
    
    E(showme.looped)$width <- 1+E(showme.looped)$weight/50
    V(showme.looped)$nombres<-colnames(proportiontable.trans)
    
    # Compute approx self-transition and use that to set node size:
    selfsizesscaled<-(diag(proportiontable.trans))*100+1
    V(showme)$selfsize<-sqrt(selfsizesscaled)
    V(showme)$size <- ((selfsizesscaled)^2)*9
    
    V(showme.looped)$selfsize<-log(selfsizesscaled)+1
    
    
    
    #####################################################
    # Network Clustering
    # l <- layout.fruchterman.reingold(showme.looped)
    # l <- layout_with_kk(showme.looped)
    # l <- layout_with_lgl(showme.looped)
    l <- layout_with_graphopt(showme.looped)
    l2 <-layout_with_graphopt(showme)
    #l <- layout_nicely(showme.looped)
    #behaviorplotlist
    
    
    #clp <- cluster_label_prop(showme.looped)
    clp <-cluster_label_prop(showme.looped.undirected) # input graph should be undirected to make sense.
    clp2<-cluster_label_prop(showme.undirected)
    # setwd("C:/Users/Rusty/Amazon Drive/BOP/Video_analysis")
    # source("PictoGrammerMega_Dec13ab.R", chdir = F)
    # #PictoGrammer<-function(proportiontable,hangingbird)
    #
    
    
    
    
    ################################################################################################################################
    #PictoGrammer(proportiontable.trans,hangingbird,uniquebehaviors)
    
    origwd<-dir
    setwd(dir<-"C:/Users/Rusty/Amazon Drive/BOP/Video_analysis/Pictograms")
    
    pictogram.names<-list.files(dir,full.names=FALSE,pattern=".png")
    
    proportiontable<-proportiontable
    hangingbird<-hangingbird
    
    
    selflooppoints<-read.table("SelfLoop.csv",sep=",",header=F)
    
    #Creates empty list, then loops through pictogram directory, pulling in each PG, adding it to the list, and naming it
    piclist<-list()
    for (pg in 1:56){
      img <- readPNG(pictogram.names[pg])  
      
      if(pg>23 && pg <42) {
        w <- matrix(rgb(img[,,1],img[,,2],img[,,3], img[,,4] * 0.75), nrow=dim(img)[1]) #0.5 is alpha
        rimg <- as.raster(w) # raster multilayer object
      } else {
        rimg <- as.raster(img) # raster multilayer object  
      }
      pgname<-gsub('.{4}$', '',pictogram.names[pg])
      piclist[[pg]]<-rimg
      names(piclist)[[pg]]<-pgname
    }
    
    #Loops through each unique behavior (from proportiontable) and creates a new, unique, behavioral pictogram
    behaviorplotlist<-list()
    par(bg="transparent")
    plot.new()
    goo=1
    #pdf("Behav30-trans-fix.pdf",width= 12, height= 12,family="NimbusRom")
    #op <- par(mfrow=c(6,6))# rows, columns
    #par(mar=c(1,1,1,1)) 
    
    
    colorweights<-E(showme.looped)$weight/min(E(showme.looped)$weight) 
    colorweights<-log(colorweights+0.01)+1
    
    Lab.palette<-colorRampPalette(c("blue",  "red"),space = "Lab")
    #  Lab.palette<-colorRampPalette(c("#ffeda0","#feb24c","#f03b20"),space = "Lab")
    
    
    E(showme.looped)$color<-Lab.palette(max(colorweights))[colorweights]
    
    col.tr <- grDevices::adjustcolor(E(showme.looped)$color, alpha=1) #adds translucency to arrows
    
    SELFCOLORS<-log(diag(proportiontable.trans)/min(E(showme.looped)$weight))+1 
    SELFCOLORS.colors<-SELFCOLORS
    #Little loop that adds same heatmap scaled colors to self-transitions, for use in behavior boxes
    for(ggg in 1:length(SELFCOLORS)){
      selfvalue<-SELFCOLORS[ggg]
      if(selfvalue>0){
        SELFCOLORS.colors[ggg]<-Lab.palette(max(colorweights))[selfvalue]
      } else {
        SELFCOLORS.colors[ggg]<-NA
      }
      
    }
    
    
    for (ub in 1:uniquebehaviors.filtered1){
      par(bg="transparent")
      graphics::plot(goo, type ="l", xlab="", ylab="", xlim=c(0, 100), ylim=c(0, 100),axes=FALSE)
      
      #53
      if( (length(grep("SS1", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`SS1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #54
      if( (length(grep("SS2", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`SS2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      ##########################################################  
      #2
      if( (length(grep("BP1", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`BP1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #1
      if( (length(grep("BP1.BP2", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`BP1-BP2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #3
      if( (length(grep("BP3", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`BP3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #4
      if( (length(grep("O1", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`O1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #5
      if( (length(grep("O2", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      ###############################################
      
      #If hangingbird  
      if(hangingbird==1){  
        #6
        if( (length(grep("O5.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1){
          rasterImage(piclist$`O5-O3-O4-MO1or2a`, xleft = 1, ybottom=1, xright=100,ytop=90)  
        }
        #7
        if( (length(grep("O5.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
          rasterImage(piclist$`O5-O3-O4a`, xleft = 1, ybottom=1, xright=100,ytop=90)  
        }
        
        #10
        if( (length(grep("O6.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1){
          rasterImage(piclist$`O6-O3-O4-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
        }
        #11
        if( (length(grep("O6.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
          rasterImage(piclist$`O6-O3-O4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
        }  
      }
      
      #If non-hangingbird
      if(hangingbird==0){  
        #8
        if( (length(grep("O5.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1){
          rasterImage(piclist$`O5-O3-O4-MO1or2b`, xleft = 1, ybottom=1, xright=100,ytop=90)  
        }
        #9
        if( (length(grep("O5.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
          rasterImage(piclist$`O5-O3-O4b`, xleft = 1, ybottom=1, xright=100,ytop=90)  
        }
        
        #10.a
        if( (length(grep("O6.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1){
          rasterImage(piclist$`O6-O3-O4-MO1or2a`, xleft = 1, ybottom=1, xright=100,ytop=90)  
        }
        #11.a
        if( (length(grep("O6.O3.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0){
          rasterImage(piclist$`O6-O3-O4a`, xleft = 1, ybottom=1, xright=100,ytop=90)  
        }  
      }
      
      #12
      if( (length(grep("O5.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1 && (length(grep("O3.O4.O5", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O5-O4-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #13
      if( (length(grep("O5.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0 && (length(grep("O3.O4.O5", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O5-O4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }  
      #14
      if( (length(grep("O6.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1 && (length(grep("O3.O4.O6", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O6-O4-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #15
      if( (length(grep("O6.O4", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0 && (length(grep("O3.O4.O6", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O6-O4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      
      #18
      if( (length(grep("O5", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1 && (length(grep("O4", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O5-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #19
      if( (length(grep("O5", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0 && (length(grep("O4", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O5`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #20
      if( (length(grep("O6", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==1 && (length(grep("O4", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O6-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #21
      if( (length(grep("O6", names(proportiontable)[ub])))==1 && (length(grep("MO", names(proportiontable)[ub])))==0 && (length(grep("O4", names(proportiontable)[ub])))==0){
        rasterImage(piclist$`O6`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      ###############################################  
      #22
      if( (length(grep("OPABB", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPABB`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }  
      #23
      if( (length(grep("OPABF", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPABF`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      } 
      #24
      if( (length(grep("OPABH", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPABH`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #25
      if( (length(grep("OPABT", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPABT`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #26
      if( (length(grep("OPABW", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPABW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #27
      if( (length(grep("OPAC1", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPAC1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #28
      if( (length(grep("OPAC2", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPAC2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #29
      if( (length(grep("OPAC3", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPAC3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #30
      if( (length(grep("OPAH", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPAH`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #31
      if( (length(grep("OPAMTT", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPAMTT`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #32
      if( (length(grep("OPAMW", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPAMW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #33
      if( (length(grep("OPATW", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPATW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #34
      if( (length(grep("OPAW", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPAW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #35
      if( (length(grep("OPMB", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPMB`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #36
      if( (length(grep("OPMF", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPMF`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #37
      if( (length(grep("OPMH", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPMH`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #38
      if( (length(grep("OPMT", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPMT`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #39
      if( (length(grep("OPMW", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`OPMW`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      ######################################
      #40
      if( (length(grep("PC1", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`PC1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #41
      if( (length(grep("PC2", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`PC2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #42
      if( (length(grep("PC3", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`PC3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #43
      if( (length(grep("PC4", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`PC4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #44
      if( (length(grep("PU1", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`PU1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      ######################################  
      #45
      if( (length(grep("RM1", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`RM1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #46
      if( (length(grep("RM2", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`RM2`, xleft = 1, ybottom=1, xright=100,ytop=90)    
      }
      #47
      if( (length(grep("RM3", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`RM3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #48
      if( (length(grep("RM4", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`RM4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      ######################################  
      #49
      if( (length(grep("SP1", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`SP1`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }
      #50
      if( (length(grep("SP2", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`SP2`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }  
      #51
      if( (length(grep("SP3", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`SP3`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      } 
      #52
      if( (length(grep("SP4", names(proportiontable)[ub])))==1){
        rasterImage(piclist$`SP4`, xleft = 1, ybottom=1, xright=100,ytop=90)  
      }  
      
      ############  
      text(50,88,ub,cex=1.5) #ADD BEHAVIOR ID NUMBER
      ########################################### 
      ############  
      # box(which="plot",lty="solid",lwd=30,col=SELFCOLORS.colors[ub]) #ADD HEATMAPPED COLORED BOX AROUND PLOT, CORRESPONDING TO LOOPS
      arrow.head.color<-ifelse(is.na(SELFCOLORS.colors[ub]),NA,SELFCOLORS.colors[ub])
      
      #   diagram::curvedarrow(from=c(45,80), to=c(55,80), lwd = 7, lty = 1, 
      #                        # lcol = "black", arr.col = "black", 
      #                        lcol = SELFCOLORS.colors[ub], arr.col = NA,
      #                        arr.pos = 1, curve = -2, dr = 0.01, arr.type="triangle",endhead=TRUE,
      #                        segment = c(0, 1))
      #   
      # polygon(x=c(58,52,55),y=c(86,86,75),col=arrow.head.color,border=arrow.head.color)
      if(!(is.na(arrow.head.color))) {
        polygon(x=30+0.4*selflooppoints[,1],y=79+0.25*selflooppoints[,2],col=arrow.head.color)
      }
      ########################################### 
      setwd("C:/Users/Rusty/Amazon Drive/BOP/Video_analysis/Pictograms/temp")
      pngplot <- recordPlot() 
      xyz<-names(proportiontable)[ub]
      individualbehaviorname<-paste(xyz,"png",sep=".")
      
      png(file = individualbehaviorname, bg = "transparent")
      replayPlot(pngplot)
      dev.off()
      
      completepicture<-readPNG(individualbehaviorname)
      rimg <- as.raster(completepicture) # raster multilayer object
      
      newpgname<-gsub('.{4}$', '',xyz)
      behaviorplotlist[[ub]]<-rimg
      names(behaviorplotlist)[[ub]]<-xyz
      
      setwd(dir)
    }
    
    #dev.off()
    
    
    
    #  behaviorplotlist is a named list corresponding to all unique behaviors from a given trial
    
    
    
    setwd(origwd)
    
    
    
    
    
    
    
    
    
    # NETWORK SUMMARY VARIABLES
    
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    #6.1 Density
    
    #The proportion of present edges from all possible edges in the network.
    
    densityscore.filtered1<-ecount(showme.looped)/(vcount(showme.looped)*(vcount(showme.looped)-1)) #for a directed network
    
    #Mean degree
    deg <- degree(showme.looped,loops=FALSE)
    mean.degree.filtered1<-mean(deg)
    
    #Mean path length
    mean_path_length.filtered1<-mean_distance(showme.looped, directed=T)
    
    
    #network diameter is the longest geodesic distance (length of the shortest path between two nodes) in the network.
    network.diameter.filtered1<-diameter(showme.looped, directed=T,weights=NA)
    
    
    #average clustering coefficient
    avg.cluster.coeff.filtered1<-transitivity(showme.looped)
    
    ############################################
    #SMALL WORLDNESS --- video-wise, filtered
    #number of nodes/vertices in graph
    vertices<- uniquebehaviors.filtered1
    #number of edges in G(n,m) graph
    edges<- sum(TransitionMatrix.cum.filtered1!=0)
    
    rando.network.filtered1<-sample_gnm(n=vertices, m=edges, directed = TRUE, loops = TRUE)
    
    Trobserved<-avg.cluster.coeff.filtered1
    mean.Trrandom<-transitivity(rando.network.filtered1)
    
    SPobserved<-mean_path_length.filtered1
    mean.SPrandom<-mean_distance(rando.network.filtered1, directed=T)  
    
    Smallworldness.filtered1<- (Trobserved/mean.Trrandom)/(SPobserved/mean.SPrandom)
    ############################################
    
    
    printme<-data.frame(FileisNamed,belongsto,TotalTimeWatched,behaviorsscored,n.unique.behavs,
                        uniquebehaviors,densityscore,mean.degree,mean_path_length,network.diameter,avg.cluster.coeff,
                        iSimpson,relativebehavioraldiversity,iSimpsonT,relativetransitiondiversity,
                        PropSelfTransitioning,Smallworldness,
                        
                        uniquebehaviors.filtered1,densityscore.filtered1,mean.degree.filtered1,mean_path_length.filtered1,network.diameter.filtered1,avg.cluster.coeff.filtered1,
                        iSimpson.filtered1,relativebehavioraldiversity.filtered1,iSimpsonT.filtered1,relativetransitiondiversity.filtered1,
                        PropSelfTransitioning.filtered1,Smallworldness.filtered1,
                        Nmales,Nfemales,Nfemales.disp,Time.upright,Time.inverted,Time.orientedtofemale,
                        femaleorientedbouts,stateANDduration,odd.behaviors,off.ss.eq,happenedwhileoff)
    
    
    
    
    
    
    #CONNECTIVITY measures reflect degree to which network differs from complete network
    #Edge density: % of edges compared to maximum
    #Average degree: Avg. number of links
    #Average path length: avg of shortest pathes between reachable nodes
    #Network diameter: longest of shortest paths
    
    #CENTRALITY measures quantify heterogeneity in network structure
    #Average clustering coefficient:
    #Components
    
    
    
    
    
    
    
    
    
    
    
    #V(showme.looped)$size<-80
    netname<-paste(name2,"Network_filtered1_cust.pdf",sep='_')
    
    setwd(location)
    pdf(netname,width= dimensions[1], height= dimensions[2],family="NimbusRom")
    op <- par(mfrow=c(1,1))# rows, columns
    par(mar=c(0,0,0,0)) 
    # plot(clp, showme.looped,edge.curved=.1,edge.arrow.size = 0.45,vertex.size=18,layout=l,rescale=FALSE,
    #      xlim=range(l[,1]),ylim=range(l[,2]))
    
    
    # plot(clp, showme.looped,edge.curved=.1,edge.arrow.size = 0.45,vertex.size=18,layout=l,rescale=TRUE,
    #      xlim=c(-1,1),ylim=c(-1,1))
    
    V(showme.looped)$raster <- behaviorplotlist
    V(showme)$raster <- behaviorplotlist
    
    
    
    
    # print(  plot(clp, showme.looped,edge.curved=.1,layout=l,rescale=TRUE,vertex.shape="raster",
    #        vertex.label=NA, 
    #        #vertex.size=15*(V(showme.looped)$selfsize), # resizes vertices based on self-self transitions
    #        vertex.size=25,
    #        edge.color=col.tr,
    #        edge.arrow.size=0.35,
    #        
    #        edge.width=3,
    #        edge.arrow.width=1.6,
    #        edge.loop.angle=1.5,
    #        
    #        xlim=c(-1,1),ylim=c(-1,1)))
    
    print(  plot(clp, showme,edge.curved=.1,layout=l,rescale=TRUE,vertex.shape="raster",
                 vertex.label=NA, 
                 #vertex.size=15*(V(showme.looped)$selfsize), # resizes vertices based on self-self transitions
                 vertex.size=25,
                 edge.color=col.tr,
                 edge.arrow.size=0.35,
                 
                 edge.width=3,
                 edge.arrow.width=1.6,
                 edge.loop.angle=1.5,
                 
                 xlim=c(-1,1),ylim=c(-1,1)))
    
    print(text(0.9,-0.75,paste("Timescored(s) = ",time.value,sep=""),col="black",cex=0.75,adj=(0)))
    print(text(0.9,-0.8,paste("behaviors = ",uniquebehaviors.filtered1,sep=""),col="black",cex=0.75,adj=(0)))  
    print(text(0.9,-0.85,paste("mean deg = ",round(mean.degree.filtered1,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
    print(text(0.9,-0.9,paste("density score = ",round(densityscore.filtered1,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
    print(text(0.9,-0.95,paste("mean path L = ",round(mean_path_length.filtered1,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
    print(text(0.9,-1,paste("network.diameter = ",network.diameter.filtered1,sep=""),col="black",cex=0.75,adj=(0)))  
    print(text(0.9,-1.05,paste("Smallworldness = ",Smallworldness.filtered1,sep=""),col="black",cex=0.75,adj=(0)))
    
    #tkplot(clp, showme.looped)
    #####  par("usr")
    
    # library(plotrix)
    # 
    # scalelist<-list()
    # gettinsmaller<-seq(from=1,to=0.52, by=-0.02)
    # for (vc in 1:25){
    # 
    #   
    #   gettingsmaller<-gettinsmaller[vc]
    #   pointlocations<-cbind(rescale(l[,1],gettingsmaller*c(-1,1)),rescale(l[,2],gettingsmaller*c(-1,1)))
    #   par(new=TRUE)
    #   plot(pointlocations,xlab="", ylab="",xlim=c(-1,1),ylim=c(-1,1),pch=vc)
    # 
    # }
    # dev.off()
    
    
    
    # for(bbb in 1:uniquebehaviors){
    #   rasterImage(behaviorplotlist[[bbb]], xleft = pointlocations[bbb,1], xright=pointlocations[bbb,1]+0.15,
    #               ybottom=pointlocations[bbb,2], ytop=pointlocations[bbb,2]+0.15)  
    #   par(new=TRUE)
    #   plot(pointlocations,xlab="", ylab="",xlim=c(-1,1),ylim=c(-1,1))
    # }
    # 
    # pointlocations<-cbind(rescale(l[,1],0.92*c(-1,1)),rescale(l[,2],0.92*c(-1,1)))
    # pointlocations<-pointlocations
    # 
    # for(bbb in 1:uniquebehaviors){
    #   rasterImage(behaviorplotlist[[bbb]], xleft = pointlocations[bbb,1], xright=pointlocations[bbb,1]+0.15,
    #               ybottom=pointlocations[bbb,2], ytop=pointlocations[bbb,2]+0.15)  
    #   par(new=TRUE)
    #   
    # }
    
    dev.off()
    
    #Another network program
    #geph
  } else {
    c.no.OFFs<-as.character(no.OFFs.filtered1)
    TransitionMatrix.cum.filtered1<-  createSequenceMatrix(c.no.OFFs,sanitize=FALSE)
    
    densityscore.filtered1<-mean.degree.filtered1<-mean_path_length.filtered1<-network.diameter.filtered1<-avg.cluster.coeff.filtered1<-iSimpson.filtered1<-relativebehavioraldiversity.filtered1<-iSimpsonT.filtered1<-relativetransitiondiversity.filtered1<-PropSelfTransitioning.filtered1<-Smallworldness.filtered1<-NA
    
    printme<-data.frame(FileisNamed,belongsto,TotalTimeWatched,behaviorsscored,n.unique.behavs,
                        uniquebehaviors,densityscore,mean.degree,mean_path_length,network.diameter,avg.cluster.coeff,
                        iSimpson,relativebehavioraldiversity,iSimpsonT,relativetransitiondiversity,
                        PropSelfTransitioning,Smallworldness,
                        
                        uniquebehaviors.filtered1,densityscore.filtered1,mean.degree.filtered1,mean_path_length.filtered1,network.diameter.filtered1,avg.cluster.coeff.filtered1,
                        iSimpson.filtered1,relativebehavioraldiversity.filtered1,iSimpsonT.filtered1,relativetransitiondiversity.filtered1,
                        PropSelfTransitioning.filtered1,Smallworldness.filtered1,
                        Nmales,Nfemales,Nfemales.disp,Time.upright,Time.inverted,Time.orientedtofemale,
                        femaleorientedbouts,stateANDduration,odd.behaviors,off.ss.eq,happenedwhileoff)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  ######################################################################################################################################
  #####################################################################################################################################
  
  ###PLOTTING PORTION OF SCRIPT
  #####################################################################################################################################
  EventBehaviors<-c("BP2","O1","O2","OPMH","OPMB","OPMF","OPMW","OPMT","OPAH","OPAC1","OPAC2","OPAC3","OPAW","OPATW","OPAMW","OPAMTT",
                    "OPABH","OPABB","OPABF","OPABW","OPABT","PC1","PC2","PC3","PC4","MO1","SP1","SP2","SP3","SP4")
  
  
  DurationBehaviors2<-c("BP1","BP3","SS1","SS2","O3","O4","O5","OFF","O6","O7","O8","OPMH.D","OPMB.D","OPMF.D","OPMW.D","OPMT.D","RepWing","RepTail","RepHead","RepTorso",
                        "RM1","RM2","RM3","RM4","MO2","PU1")
  SPDurations.start<-c("SP1","SP2","SP3")
  SPDurations.stop<-c("SP2","SP3","SP4")
  
  forplot<-matrix(nrow=(TotalTime*10+3),ncol=length(c(DurationBehaviors2,EventBehaviors))) #Creates empty matrix with 1 column/behavior and 10 rows/sec watched (0.1s resolution) (+2 accounts for 0 and buffer at end)
  colnames(forplot)<-c(DurationBehaviors2,EventBehaviors)
  
  fulltime<-as.character(seq(0,TotalTime+0.2,0.1))
  rownames(forplot)<-fulltime
  forplot<-data.frame(forplot)
  
  q.time<-as.character(round(table$time/0.1)*0.1)
  table<-cbind(table,q.time)
  
  #eventdata<-table[which(table$code %in% EventBehaviors),]
  
  for (jj in EventBehaviors){
    eventtable<-table[which(table$code %in% EventBehaviors),]
    #behav.obs <- eventtable[grep(m, eventtable$code), ]
    use<-subset(eventtable,code==jj)
    use<-unique(use)
    rowlist<-match(use$q.time,rownames(forplot))
    forplot[rowlist,jj]<-1 #Adds 1 to every 0.1 sec time window that matches for each behavior at each time
    
  }
  
  
  for (kk in DurationBehaviors2){
    use<-table[grep(kk, table$code), ]
    use<-unique(use)
    bouts<-nrow(use)/2
    if(bouts>0.6){#Only if there are >0 bouts, do the following (adding in 1s during the time periods when each duration behavior is active)
      for (ll in 1:bouts){
        start<-(ll-1)*2+1
        end<-start+1
        subuse<-use[start:end,]
        listotimes<-as.numeric(as.character(subuse$q.time))
        timeon<-seq(listotimes[1],listotimes[2],0.1)
        timeon2<-as.character(timeon)
        
        matchingrows<-match(timeon2,rownames(forplot)) #Can return non-matches for true matches based on floating point values
        forplot[matchingrows,kk]<-1 #Adds 1 to every 0.1 sec time window that matches for each behavior at each time
        
      }
    }
    
  }
  
  
  for (mm in 1:length(SPDurations.start)){
    use.starts<-table[grep(SPDurations.start[mm], table$code), ]
    use.ends<-table[grep(SPDurations.stop[mm], table$code), ]
  }
  
  
  DistinctColors<-c("#ae7b4a","#6264de","#48cd60","#9744bf","#50a428","#da5ac9","#91bf3c","#be77e6","#5ebd63","#a62d8a","#4dcc93","#d84498","#368939","#7753b5","#b4b43b","#587ee4",
                    "#e1aa34","#695faf","#e8792d","#4da9d9","#d04228","#4ed6c7","#db3973","#319a6b","#d3384c","#3bbac6","#e56d5c","#57b597","#e482d5","#69882b","#a756a7","#82bc75",
                    "#dd6a93","#2d7350","#b197e3","#a68222","#7399dc","#c5792b","#4167a5","#d5ae69","#7f5e9c","#aab56e","#984167","#5c894d","#e293c1","#5a6624","#a86597","#8b8239",
                    "#a74851","#e7986b","#7f5724","#dd8383","#a84e28")
  
  
  notimes<-forplot
  
  offperiod<-notimes[,c("OFF")]
  offperiod[is.na(offperiod)]<-0 #REPLACE NAs with ZERO
  
  notimes[,c("OFF")]<-NULL #removes fulltime and  OFF columns
  #remove columns with NO data
  columnstoremove<-c()
  for(rr in 1:ncol(notimes)){
    if(sum(notimes[,rr],na.rm = TRUE)==0){
      columnstoremove<-c(columnstoremove,rr)
    }
  }
  cleaned<-notimes[,-c(columnstoremove)]
  cleaned[is.na(cleaned)]<-0 #REPLACE NAs with ZERO
  
  str(cleaned)
  nnn.pred <- nrow(cleaned)
  behavs<-ncol(cleaned)
  date.seq <- seq(from=0, to=1, length=nnn.pred)
  rowtot <- apply(cleaned, 1, sum, na.rm=T )
  coltot <- apply(cleaned, 2, sum, na.rm=T )
  
  
  modifiedtime<-file.info(filename)[5]
  
  
  shadowtext <- function(x, y=NULL, labels, col='white', bg='black', 
                         theta= seq(0, 2*pi, length.out=50), r=0.1, ... ) {
    
    xy <- xy.coords(x,y)
    xo <- r*strwidth('A')
    yo <- r*strheight('A')
    
    # draw background text with small shift in x and y in background colour
    for (i in theta) {
      text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
    }
    # draw actual text in exact xy position in foreground colour
    text(xy$x, xy$y, labels, col=col, ... )
  }
  
  
  
  
  height<-behavs/50
  
  
  
  
  #pdf("C",width= 8, height= 12,family="NimbusRom")
  prefix<-paste("Unequal",odd.behaviors,sep="_")
  newunequalname<-paste(prefix,name,sep="_")
  newunequalname2<-paste(location,newunequalname,sep="/")
  
  newname<-ifelse(lengthodds<1,fullname,newunequalname2)
  
  
  
  pdf(newname,width= 8, height= 12,family="NimbusRom")
  # Cummulative Value Plot
  par(mar=c(1,0,0,0)) 
  op <- par(mfrow=c(behavs+1,1))# rows, columns (one more than behavs in 'cleaned', to include OFF periods)
  
  #Plot OFF periods at top
  yyy.old <- offperiod
  yyy <- c( yyy.old, rev(rep(0, length(yyy.old))) )
  
  behaviorname<-"OffScreen"
  
  print(plot( seq(0, 1, nnn.pred), seq(0, 1, nnn.pred),type="n",xlim = c(-.05,1.0),ylim = c(0,1.0),xlab="", ylab="",cex.main = 1.0,main = "",axes = F,col.main = "black" ))
  # Below 
  print(axis(1,at = c(0.00,1),labels=c("","")))
  
  print(axis(2, at = c(0, 0.25, 0.5, 0.75, 1.0), labels = c("0", "", "0.5","", "1")))	
  
  print(box(lwd=2))
  xxx <- c( date.seq, rev(date.seq))
  print(polygon(xxx, yyy, col="black", border=NA))
  OFFxxs<-xxx
  OFFyys<-yyy
  #text(0,0.5,behaviorname,,col="black",cex=1.75,font=4)
  #text(0,0.5,behaviorname,bg="black",col=DistinctColors[iii.nlcd],cex=1.5,font=4)
  print(shadowtext(-0.09,height,behaviorname,col="white",cex=2,adj=(0)))
  print(text(-0.09,.95,modifiedtime[1,1],col="black",cex=0.75,adj=(0)))
  print(text(-0.09,.65,name,col="black",cex=0.75,adj=(0)))
  
  summarytextsize<-15/behavs
  print(shadowtext(1,0.95,paste("Duration",TotalTime,sep=" "),col="red",cex=summarytextsize,adj=(1)))
  print(shadowtext(1,0.75,paste("Unequal",odd.behaviors,sep="="),col="red",cex=summarytextsize,adj=(1)))
  print(shadowtext(1,0.55,paste("SSequal",off.ss.eq,sep="="),col="red",cex=summarytextsize,adj=(1)))
  print(shadowtext(1,0.35,paste("WhileOff",happenedwhileoff,sep="="),col="red",cex=summarytextsize,adj=(1)))
  print(shadowtext(1,0.15,paste("Nmales",Nmales,sep="="),col="red",cex=summarytextsize,adj=(1)))
  
  
  print(axis(1,at = c(0.5),labels=""))
  # print(text(0.489,0.10,round(TotalTime/2,digits=1)))
  
  
  for (iii.nlcd in 1:behavs){
    yyy.old <- cleaned[,iii.nlcd]
    yyy <- c( yyy.old, rev(rep(0, length(yyy.old))) )
    
    behaviorname<-colnames(cleaned)[iii.nlcd]
    
    print(plot( seq(0, 1, nnn.pred), seq(0, 1, nnn.pred),type="n",xlim = c(-.05,1.0),ylim = c(0,1.0),xlab="", ylab="",cex.main = 1.0,main = "",axes = F,col.main = "black" ))
    # Below 
    print(axis(1,at = c(0.00,1),labels=c("","")))
    print(axis(2, at = c(0, 0.25, 0.5, 0.75, 1.0), labels = c("0", "", "0.5","", "1")))	
    
    print(box(lwd=2))
    xxx <- c( date.seq, rev(date.seq))
    
    print(polygon(xxx, yyy, col=DistinctColors[iii.nlcd], border=NA))
    
    print(polygon(OFFxxs, OFFyys, col=rgb(.2,.2,.2,alpha=.2), border=NA))
    #text(0,0.5,behaviorname,,col="black",cex=1.75,font=4)
    #text(0,0.5,behaviorname,bg="black",col=DistinctColors[iii.nlcd],cex=1.5,font=4)
    print(shadowtext(-0.09,height,behaviorname,col=DistinctColors[iii.nlcd],cex=2,adj=(0)))
    print(axis(1,at = c(0.5),labels=""))
    #print(text(0.48,0.10,round(TotalTime/2,digits=1)))
    
    
  } 
  
  dev.off()
  
  
  
  ######################################################################################################################################  
  ######################################################################################################################################
  
  
  setwd(origsourcefolder)
  
  sevencomponents<-list(printme,TransitionMatrix.cum,TransitionMatrix.cum.filtered1,behaviorcounts,behaviorcounts.filtered1,six,seven)
  
  return(sevencomponents)
}
#testone and testtwo are the names of two matrices you want to combine
CombineTwoTransitionMatrices<-function(testone,testtwo){
  #tests if first matrix (testone) is of length zero---if it is, we just use testtwo
  if(length(row.names(testone))>0){
    
  shared.variables<-rownames(testtwo)[which(rownames(testtwo) %in% rownames(testone))]
  

  
  
  #tests if there are any shared variables
    if(length(shared.variables)==0){
      not.in.testone<-rownames(testtwo)
      not.in.testtwo<-rownames(testone)
      behav.structure<-c(not.in.testone,not.in.testtwo)
    } else {
      not.in.testone<-rownames(testtwo)[-which(rownames(testtwo) %in% rownames(testone))]
      not.in.testtwo<-rownames(testone)[-which(rownames(testone) %in% rownames(testtwo))]
      behav.structure<-c(shared.variables,not.in.testone,not.in.testtwo)  
    }
  
  structured.mat<-matrix(data=NA,ncol=1,nrow=length(behav.structure))
  row.names(structured.mat)<-behav.structure
  
  for(ww in 1:length(behav.structure)){
    good<-behav.structure[ww]
    origcolumn.num<-which(colnames(testone)==good)#returns only exact match
    seconmatcolumn.num<-which(colnames(testtwo)==good)
    #is this behavior ("good") present in the original (testone)?
    if(length(origcolumn.num)>0){
      #is this behavior ("good") present in the second matrix (testtwo)?
      if(length(seconmatcolumn.num)>0){
        aa<-testone[,origcolumn.num,drop=FALSE]
        bb<-testtwo[,good,drop=FALSE]
        
        a1<-aa[shared.variables,,drop=FALSE]+bb[shared.variables,,drop=FALSE]
        b1<-bb[not.in.testone,,drop=FALSE]
        c1<-aa[not.in.testtwo,,drop=FALSE]
        
        newcol<-rbind(a1,b1,c1)

        
        
        #After the first column (ww>1), merge "newcol" to newmatrix
        if(ww>1){
          newmatrix<-transform(merge(newmatrix,newcol,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
        } else {
          newmatrix<-transform(merge(structured.mat[,-1],newcol,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
        }
      } else {
        newcol<-testone[,good,drop=FALSE]
        
        #After the first column (ww>1), merge "newcol" to newmatrix
        if(ww>1){
          newmatrix<-transform(merge(newmatrix,newcol,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
        } else {
      
          newmatrix<-transform(merge(structured.mat[,-1],newcol,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
          
        }
      }
    } else { #to get to here, we must have a behavior that is not in testone
      
      newcol<-testtwo[,good,drop=FALSE]
      
      #After the first column (ww>1), merge "newcol" to newmatrix
      if(ww>1){
        newmatrix<-transform(merge(newmatrix,newcol,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
      } else {
        
        structured.mat<-matrix(data=NA,ncol=1,nrow=length(behav.structure))
        row.names(structured.mat)<-behav.structure
        newmatrix<-transform(merge(structured.mat[,-1],newcol,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
        
      }
    }
  
  }
  
  
  #re-organize and make symmetrical
  newmatrix<-newmatrix[,rownames(newmatrix)]
  newmatrix[is.na(newmatrix)] <- 0
  newmatrix<-as.matrix(newmatrix)
    
  } else {
  newmatrix<-testtwo
}  
  return(newmatrix)
  
}


#include a (nrow(MegaMatrix),b(colnames(MegaMatrix),c(length of behavcounts),d(behavcounts)))
megabehaviorsummarizer<-function(a,b,c,d){
  Megauniquebehaviors<-a
  Totalbehaviors<-matrix(data=0,ncol=a,nrow=1)
  colnames(Totalbehaviors)<-b
  
  for (tt in 1:c){
    check<-d[[tt]]
    Totalbehaviors[1,match(rownames(check),colnames(Totalbehaviors))]<-check+Totalbehaviors[1,match(rownames(check),colnames(Totalbehaviors))]
  }
  
  
  Megaproportiontable<-prop.table(Totalbehaviors)
  ##This loop creates the denominator for the calculation of the inverse Simpson Index for behavioral diversity
  val2<-0
  for (u in 1:Megauniquebehaviors){
    proptosquare<-Megaproportiontable[u]
    val <-proptosquare^2.0
    val2<-val+val2
  }			
  #inverse Simpson Index for behavioral diversity
  Mega.iSimpson<-1.0/as.numeric(val2)
  Mega.relativebehavioraldiversity<-Mega.iSimpson/(Megauniquebehaviors) #relative color diversity, accounting for # of color classes
  behavsumm<-list(Mega.iSimpson,Mega.relativebehavioraldiversity)
  names(behavsumm)<-c("Mega.iSimpson","Mega.relativebehavioraldiversity")
  
  return(behavsumm)
} 

#uses the MegaMatrix to return iSimpson transitional values
megatransitionsummarizer<-function(MegaMatrix){
  
  noselfs<-MegaMatrix
  diag(noselfs)<-NA
  
  prop.table.excludeNAs<- function(x) {x/sum(x, na.rm=TRUE)}
  proportiontable.trans<-prop.table.excludeNAs(noselfs)
  
  uniquebehaviors<-nrow(MegaMatrix)
  
  if(uniquebehaviors>1){
    val3<-matrix(nrow=(((uniquebehaviors^2)-uniquebehaviors)/2),ncol=2)
    rownum<-1
    for (u in 1:uniquebehaviors){
      for(v in u:uniquebehaviors){#iteratively reduces columns analyzed in next loop to avoid double counting transitions 
        if (u!=v){
          val3[rownum,2]<-as.numeric(proportiontable.trans[u,v]+proportiontable.trans[v,u])
          val3[rownum,1]<-paste(row.names(proportiontable.trans)[u],colnames(proportiontable.trans)[v],sep="-")
          rownum<-rownum+1
        }
      }
    }
    val4<-as.data.frame(val3)
    val4[,2]<-as.numeric(as.character(val4[,2])) 
    
    simpsonT<-0
    
    for(z in 1:nrow(val4)){
      propsq<-val4[z,2]^2.0
      simpsonT<-simpsonT+propsq
    }
    
    
    iSimpsonT<-1.0/simpsonT #Calculates inverse Simpson index for horizontal transitions
    relativetransitiondiversity<-iSimpsonT/(uniquebehaviors*(uniquebehaviors-1)/2)  
  } else {
    
    iSimpsonT <- NA
    relativetransitiondiversity<- NA
  }
  
  PropSelfTransitioning<-sum(diag(MegaMatrix))/sum(MegaMatrix)
  
  transsumm<-list(iSimpsonT,relativetransitiondiversity,PropSelfTransitioning)
  names(transsumm)<-c("Mega.iSimpson.T","Mega.relativeTransitiondiversity","PropSelfTransitioning")
  
  return(transsumm)
  
  
}



#uses MegaMatrix to calculate species-wide network + summaries
MegaNetwork<-function(MegaMatrix,hangingbird,species.name,location,filtered){
  
  
  proportiontable<-prop.table(MegaMatrix)
  
  uniquebehaviors<-nrow(MegaMatrix)
  
  time.value<-sum(MegaMatrix)/10 #seconds watched
  
  network.TransitionMatrix.cum<-network(proportiontable,directed=TRUE,loops=TRUE,hyper=FALSE,matrix.type="adjacency")
  
  
  
  showme<-graph_from_adjacency_matrix(proportiontable,mode="directed",weighted=TRUE,diag=FALSE,add.colnames=TRUE)
  showme.undirected<-graph_from_adjacency_matrix(proportiontable,mode="undirected",weighted=TRUE,diag=FALSE,add.colnames=TRUE)
  
  showme.looped<- graph_from_adjacency_matrix(proportiontable,mode="directed",weighted=TRUE,diag=TRUE,add.colnames=TRUE)
  showme.looped.undirected<- graph_from_adjacency_matrix(proportiontable,mode="undirected",weighted=TRUE,diag=TRUE,add.colnames=TRUE)
  
  
  
  E(showme)$width <- 1+E(showme)$weight/50
  V(showme)$nombres<-colnames(proportiontable)
  
  E(showme.looped)$width <- 1+E(showme.looped)$weight/50
  V(showme.looped)$nombres<-colnames(proportiontable)
  
  # Compute approx self-transition and use that to set node size:
  selfsizesscaled<-(diag(proportiontable))*100+1
  V(showme)$selfsize<-sqrt(selfsizesscaled)
  V(showme)$size <- ((selfsizesscaled)^2)*9
  
  V(showme.looped)$selfsize<-log(selfsizesscaled)+1
  
  
  
  #####################################################
  # Network Clustering
  # l <- layout.fruchterman.reingold(showme.looped)
  # l <- layout_with_kk(showme.looped)
  # l <- layout_with_lgl(showme.looped)
  l <- layout_with_graphopt(showme.looped)
  l2 <-layout_with_graphopt(showme)
  #l <- layout_nicely(showme.looped)
  #behaviorplotlist
  
  
  #clp <- cluster_label_prop(showme.looped)
  clp <-cluster_label_prop(showme.looped.undirected) # input graph should be undirected to make sense.
  clp2<-cluster_label_prop(showme.undirected)
  # setwd("C:/Users/Rusty/Amazon Drive/BOP/Video_analysis")
  # source("PictoGrammerMega_Dec13ab.R", chdir = F)
  # #PictoGrammer<-function(proportiontable,hangingbird)
  #
  
  
  
  
  ################################################################################################################################
  #PictoGrammer(proportiontable,hangingbird,uniquebehaviors)
  
  origwd<-dir
  setwd(dir<-"C:/Users/Rusty/Amazon Drive/BOP/Video_analysis/Pictograms")
  
  pictogram.names<-list.files(dir,full.names=FALSE,pattern=".png")
  
  proportiontable<-proportiontable
  hangingbird<-hangingbird
  
  selflooppoints<-read.table("SelfLoop.csv",sep=",",header=F)
  
  #Creates empty list, then loops through pictogram directory, pulling in each PG, adding it to the list, and naming it
  piclist<-list()
  for (pg in 1:56){
    img <- readPNG(pictogram.names[pg])  
    
    if(pg>23 && pg <42) {
      w <- matrix(rgb(img[,,1],img[,,2],img[,,3], img[,,4] * 0.75), nrow=dim(img)[1]) #0.5 is alpha
      rimg <- as.raster(w) # raster multilayer object
    } else {
      rimg <- as.raster(img) # raster multilayer object  
    }
    pgname<-gsub('.{4}$', '',pictogram.names[pg])
    piclist[[pg]]<-rimg
    names(piclist)[[pg]]<-pgname
  }
  
  #Loops through each unique behavior (from proportiontable) and creates a new, unique, behavioral pictogram
  behaviorplotlist<-list()
  par(bg="transparent")
  plot.new()
  goo=1
  #pdf("Behav30-trans-fix.pdf",width= 12, height= 12,family="NimbusRom")
  #op <- par(mfrow=c(6,6))# rows, columns
  #par(mar=c(1,1,1,1)) 
  
  
  colorweights<-E(showme.looped)$weight/min(E(showme.looped)$weight) 
  colorweights<-log(colorweights+0.01)+1
  
  Lab.palette<-colorRampPalette(c("blue",  "red"),space = "Lab")
  #  Lab.palette<-colorRampPalette(c("#ffeda0","#feb24c","#f03b20"),space = "Lab")
  
  
  E(showme.looped)$color<-Lab.palette(max(colorweights))[colorweights]
  
  col.tr <- grDevices::adjustcolor(E(showme.looped)$color, alpha=1) #adds translucency to arrows
  
  SELFCOLORS<-log(diag(proportiontable)/min(E(showme.looped)$weight))+1 
  SELFCOLORS.colors<-SELFCOLORS
  #Little loop that adds same heatmap scaled colors to self-transitions, for use in behavior boxes
  for(ggg in 1:length(SELFCOLORS)){
    selfvalue<-SELFCOLORS[ggg]
    if(selfvalue>0){
      SELFCOLORS.colors[ggg]<-Lab.palette(max(colorweights))[selfvalue]
    } else {
      SELFCOLORS.colors[ggg]<-NA
    }
    
  }
  
  
  for (ub in 1:uniquebehaviors){
    par(bg="transparent")
    graphics::plot(goo, type ="l", xlab="", ylab="", xlim=c(0, 100), ylim=c(0, 100),axes=FALSE)
    
    #53
    if( (length(grep("SS1", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`SS1`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #54
    if( (length(grep("SS2", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`SS2`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    ##########################################################  
    #2
    if( (length(grep("BP1", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`BP1`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #1
    if( (length(grep("BP1.BP2", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`BP1-BP2`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #3
    if( (length(grep("BP3", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`BP3`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #4
    if( (length(grep("O1", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`O1`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #5
    if( (length(grep("O2", colnames(proportiontable)[ub])))==1 && (length(grep("MO", colnames(proportiontable)[ub])))==0){
      rasterImage(piclist$`O2`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    ###############################################
    
    #If hangingbird  
    if(hangingbird==1){  
      #6
      if( (length(grep("O5.O3.O4", colnames(proportiontable)[ub])))==1 && (length(grep("MO", colnames(proportiontable)[ub])))==1){
        rasterImage(piclist$`O5-O3-O4-MO1or2a`, xleft = 1, ybottom=1, xright=100,ytop=100)  
      }
      #7
      if( (length(grep("O5.O3.O4", colnames(proportiontable)[ub])))==1 && (length(grep("MO", colnames(proportiontable)[ub])))==0){
        rasterImage(piclist$`O5-O3-O4a`, xleft = 1, ybottom=1, xright=100,ytop=100)  
      }
      
      #10
      if( (length(grep("O6.O3.O4", colnames(proportiontable)[ub])))==1 && (length(grep("MO", colnames(proportiontable)[ub])))==1){
        rasterImage(piclist$`O6-O3-O4-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=100)  
      }
      #11
      if( (length(grep("O6.O3.O4", colnames(proportiontable)[ub])))==1 && (length(grep("MO", colnames(proportiontable)[ub])))==0){
        rasterImage(piclist$`O6-O3-O4`, xleft = 1, ybottom=1, xright=100,ytop=100)  
      }  
    }
    
    #If non-hangingbird
    if(hangingbird==0){  
      #8
      if( (length(grep("O5.O3.O4", colnames(proportiontable)[ub])))==1 && (length(grep("MO", colnames(proportiontable)[ub])))==1){
        rasterImage(piclist$`O5-O3-O4-MO1or2b`, xleft = 1, ybottom=1, xright=100,ytop=100)  
      }
      #9
      if( (length(grep("O5.O3.O4", colnames(proportiontable)[ub])))==1 && (length(grep("MO", colnames(proportiontable)[ub])))==0){
        rasterImage(piclist$`O5-O3-O4b`, xleft = 1, ybottom=1, xright=100,ytop=100)  
      }
      
      #10.a
      if( (length(grep("O6.O3.O4", colnames(proportiontable)[ub])))==1 && (length(grep("MO", colnames(proportiontable)[ub])))==1){
        rasterImage(piclist$`O6-O3-O4-MO1or2a`, xleft = 1, ybottom=1, xright=100,ytop=100)  
      }
      #11.a
      if( (length(grep("O6.O3.O4", colnames(proportiontable)[ub])))==1 && (length(grep("MO", colnames(proportiontable)[ub])))==0){
        rasterImage(piclist$`O6-O3-O4a`, xleft = 1, ybottom=1, xright=100,ytop=100)  
      }  
    }
    
    #12
    if( (length(grep("O5.O4", colnames(proportiontable)[ub])))==1 && (length(grep("MO", colnames(proportiontable)[ub])))==1 && (length(grep("O3.O4.O5", colnames(proportiontable)[ub])))==0){
      rasterImage(piclist$`O5-O4-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #13
    if( (length(grep("O5.O4", colnames(proportiontable)[ub])))==1 && (length(grep("MO", colnames(proportiontable)[ub])))==0 && (length(grep("O3.O4.O5", colnames(proportiontable)[ub])))==0){
      rasterImage(piclist$`O5-O4`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }  
    #14
    if( (length(grep("O6.O4", colnames(proportiontable)[ub])))==1 && (length(grep("MO", colnames(proportiontable)[ub])))==1 && (length(grep("O3.O4.O6", colnames(proportiontable)[ub])))==0){
      rasterImage(piclist$`O6-O4-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #15
    if( (length(grep("O6.O4", colnames(proportiontable)[ub])))==1 && (length(grep("MO", colnames(proportiontable)[ub])))==0 && (length(grep("O3.O4.O6", colnames(proportiontable)[ub])))==0){
      rasterImage(piclist$`O6-O4`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    
    #18
    if( (length(grep("O5", colnames(proportiontable)[ub])))==1 && (length(grep("MO", colnames(proportiontable)[ub])))==1 && (length(grep("O4", colnames(proportiontable)[ub])))==0){
      rasterImage(piclist$`O5-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #19
    if( (length(grep("O5", colnames(proportiontable)[ub])))==1 && (length(grep("MO", colnames(proportiontable)[ub])))==0 && (length(grep("O4", colnames(proportiontable)[ub])))==0){
      rasterImage(piclist$`O5`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #20
    if( (length(grep("O6", colnames(proportiontable)[ub])))==1 && (length(grep("MO", colnames(proportiontable)[ub])))==1 && (length(grep("O4", colnames(proportiontable)[ub])))==0){
      rasterImage(piclist$`O6-MO1or2`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #21
    if( (length(grep("O6", colnames(proportiontable)[ub])))==1 && (length(grep("MO", colnames(proportiontable)[ub])))==0 && (length(grep("O4", colnames(proportiontable)[ub])))==0){
      rasterImage(piclist$`O6`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    ###############################################  
    #22
    if( (length(grep("OPABB", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPABB`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }  
    #23
    if( (length(grep("OPABF", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPABF`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    } 
    #24
    if( (length(grep("OPABH", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPABH`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #25
    if( (length(grep("OPABT", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPABT`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #26
    if( (length(grep("OPABW", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPABW`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #27
    if( (length(grep("OPAC1", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPAC1`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #28
    if( (length(grep("OPAC2", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPAC2`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #29
    if( (length(grep("OPAC3", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPAC3`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #30
    if( (length(grep("OPAH", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPAH`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #31
    if( (length(grep("OPAMTT", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPAMTT`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #32
    if( (length(grep("OPAMW", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPAMW`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #33
    if( (length(grep("OPATW", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPATW`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #34
    if( (length(grep("OPAW", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPAW`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #35
    if( (length(grep("OPMB", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPMB`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #36
    if( (length(grep("OPMF", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPMF`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #37
    if( (length(grep("OPMH", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPMH`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #38
    if( (length(grep("OPMT", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPMT`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #39
    if( (length(grep("OPMW", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`OPMW`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    ######################################
    #40
    if( (length(grep("PC1", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`PC1`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #41
    if( (length(grep("PC2", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`PC2`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #42
    if( (length(grep("PC3", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`PC3`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #43
    if( (length(grep("PC4", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`PC4`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #44
    if( (length(grep("PU1", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`PU1`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    ######################################  
    #45
    if( (length(grep("RM1", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`RM1`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #46
    if( (length(grep("RM2", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`RM2`, xleft = 1, ybottom=1, xright=100,ytop=100)    
    }
    #47
    if( (length(grep("RM3", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`RM3`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #48
    if( (length(grep("RM4", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`RM4`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    ######################################  
    #49
    if( (length(grep("SP1", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`SP1`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }
    #50
    if( (length(grep("SP2", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`SP2`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }  
    #51
    if( (length(grep("SP3", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`SP3`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    } 
    #52
    if( (length(grep("SP4", colnames(proportiontable)[ub])))==1){
      rasterImage(piclist$`SP4`, xleft = 1, ybottom=1, xright=100,ytop=100)  
    }  
    
    ############  
    text(50,95,ub,cex=3) #ADD BEHAVIOR ID NUMBER
    ########################################### 
    ############  
  # box(which="plot",lty="solid",lwd=20,col=SELFCOLORS.colors[ub]) #ADD HEATMAPPED COLORED BOX AROUND PLOT, CORRESPONDING TO LOOPS
  #  graphics::plot(goo, type ="l", xlab="", ylab="", xlim=c(0, 100), ylim=c(0, 100),axes=FALSE)
    #add heatmapped "self" loop arrows
    arrow.head.color<-ifelse(is.na(SELFCOLORS.colors[ub]),NA,SELFCOLORS.colors[ub])
    
    #   diagram::curvedarrow(from=c(45,80), to=c(55,80), lwd = 7, lty = 1, 
    #                        # lcol = "black", arr.col = "black", 
    #                        lcol = SELFCOLORS.colors[ub], arr.col = NA,
    #                        arr.pos = 1, curve = -2, dr = 0.01, arr.type="triangle",endhead=TRUE,
    #                        segment = c(0, 1))
    #   
    # polygon(x=c(58,52,55),y=c(86,86,75),col=arrow.head.color,border=arrow.head.color)
    if(!(is.na(arrow.head.color))) {
      polygon(x=30+0.4*selflooppoints[,1],y=79+0.25*selflooppoints[,2],col=arrow.head.color)
    }
    
    

    ########################################### 
    setwd("C:/Users/Rusty/Amazon Drive/BOP/Video_analysis/Pictograms/temp")
    pngplot <- recordPlot() 
    xyz<-colnames(proportiontable)[ub]
    individualbehaviorname<-paste(xyz,"png",sep=".")
    
    png(file = individualbehaviorname, bg = "transparent")
    replayPlot(pngplot)
    dev.off()
    
    completepicture<-readPNG(individualbehaviorname)
    rimg <- as.raster(completepicture) # raster multilayer object
    
    newpgname<-gsub('.{4}$', '',xyz)
    behaviorplotlist[[ub]]<-rimg
    names(behaviorplotlist)[[ub]]<-xyz
    
    setwd(dir)
  }
  
  #dev.off()
  
  
  
  #  behaviorplotlist is a named list corresponding to all unique behaviors from a given trial
  
  
  
  setwd(origwd)
  
  
  
  
  
  
  
  
  
  # NETWORK SUMMARY VARIABLES
  
  #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
  #6.1 Density
  
  #The proportion of present edges from all possible edges in the network.
  
#densityscore<-ecount(showme.looped)/(vcount(showme.looped)*(vcount(showme.looped)-1)) #for a directed network
  densityscore<-ecount(showme.looped)/(vcount(showme.looped)^2)#for a directed network WITH self-loops
  
  #Mean degree
  deg <- degree(showme.looped,loops=FALSE)
  mean.degree<-mean(deg)
  
  #Mean path length
  mean_path_length<-mean_distance(showme.looped, directed=T)
  
  
  #network diameter is the longest geodesic distance (length of the shortest path between two nodes) in the network.
  network.diameter<-diameter(showme.looped, directed=T,weights=NA)
  
  
  #average clustering coefficient
  avg.cluster.coeff<-transitivity(showme.looped,)
  
  ############################################
  #SMALL WORLDNESS ---- MEGA NETWORKS
  #number of nodes/vertices in graph
  vertices<- nrow(MegaMatrix)
  #number of edges in G(n,m) graph
  edges<- sum(MegaMatrix!=0)
  
  rando.network<-sample_gnm(n=vertices, m=edges, directed = TRUE, loops = TRUE)
  
  Trobserved<-avg.cluster.coeff
  mean.Trrandom<-transitivity(rando.network)
  
  SPobserved<-mean_path_length
  mean.SPrandom<-mean_distance(rando.network, directed=T)  
  
  Smallworldness<- (Trobserved/mean.Trrandom)/(SPobserved/mean.SPrandom)
  ############################################
  
  
  
  printme<-data.frame(uniquebehaviors,densityscore,mean.degree,mean_path_length,network.diameter,avg.cluster.coeff,Smallworldness)
  
  if(filtered==0){
  netname<-paste(species.name,"SpeciesNetwork.pdf",sep="-")
  } else {
    netname<-paste(species.name,"SpeciesNetwork_filtered.pdf",sep="-")
  }
  
  
  setwd(location)
  pdf(netname,width= 24, height= 24,family="NimbusRom")
  op <- par(mfrow=c(1,1))# rows, columns
  par(mar=c(0,0,0,0)) 
  # plot(clp, showme.looped,edge.curved=.1,edge.arrow.size = 0.45,vertex.size=18,layout=l,rescale=FALSE,
  #      xlim=range(l[,1]),ylim=range(l[,2]))
  
  
  # plot(clp, showme.looped,edge.curved=.1,edge.arrow.size = 0.45,vertex.size=18,layout=l,rescale=TRUE,
  #      xlim=c(-1,1),ylim=c(-1,1))
  
  V(showme.looped)$raster <- behaviorplotlist
  V(showme)$raster <- behaviorplotlist
  
  
  print(  plot(clp, showme,edge.curved=.1,layout=l,rescale=TRUE,vertex.shape="raster",
               vertex.label=NA, 
               #vertex.size=15*(V(showme.looped)$selfsize), # resizes vertices based on self-self transitions
               vertex.size=25,
               edge.color=col.tr,
               edge.arrow.size=0.35,
               
               edge.width=3,
               edge.arrow.width=1.6,
               edge.loop.angle=1.5,
               
               xlim=c(-1,1),ylim=c(-1,1)))
  
  print(text(0.9,-0.75,paste("Timescored(s) = ",time.value,sep=""),col="black",cex=0.75,adj=(0)))
  print(text(0.9,-0.8,paste("behaviors = ",uniquebehaviors,sep=""),col="black",cex=0.75,adj=(0)))  
  print(text(0.9,-0.85,paste("mean deg = ",round(mean.degree,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
  print(text(0.9,-0.9,paste("density score = ",round(densityscore,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
  print(text(0.9,-0.95,paste("mean path L = ",round(mean_path_length,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
  print(text(0.9,-1,paste("network.diameter = ",network.diameter,sep=""),col="black",cex=0.75,adj=(0)))  
  print(text(0.9,-1.05,paste("Smallworldness = ",Smallworldness,sep=""),col="black",cex=0.75,adj=(0)))
  
  dev.off()
  
  
  
  return(printme)
  
}

#Uses "no.OFFs" sequence of behaviors to generate summary stats of behavioral complexity, for use in sliding window analysis
networksummarizer<-function(z){
  summarytable<-table(z)
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
  Shannon<-diversity(dftable,index = "shannon")
  propShannon<-Shannon/log(specnumber(dftable))
  redundancy<-1-propShannon
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
    
    val4.s[,2]<-ifelse(val4.s[,2]==0,NA,val4.s[,2])
    val4.s<-na.omit(val4.s)
    
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
    uniquetransitions<-nrow(val4.s)+uniquebehaviors
    
  } else {
    Shannon.T<-NA
    trueShannon.T<-NA
    propShannon.T<-NA
    redundancy.T<-NA
    uniquetransitions<-1
  }
  
  
  #########################################################
  
  PropSelfTransitioning<-sum(diag(TransitionMatrix.cum))/sum(TransitionMatrix.cum)
  
  network.TransitionMatrix.cum<-network(proportiontable.trans,directed=TRUE,loops=TRUE,hyper=FALSE,matrix.type="adjacency")
  
  # NETWORK SUMMARY VARIABLES
  showme.looped<- graph_from_adjacency_matrix(proportiontable.trans,mode="directed",weighted=TRUE,diag=TRUE,add.colnames=TRUE)
  undir.looped<- graph_from_adjacency_matrix(proportiontable.trans,mode="undirected",weighted=TRUE,diag=TRUE,add.colnames=TRUE)
  
  ########################################################
  
  # Klein DJ, Randić M. Resistance distance. Journal of Mathematical Chemistry 1993; 12: 81–95.
  # Ellens W, Kooij RE. 2013. Graph measures and network robustness. arXiv:1311.5064v1
  
  # Klein and Randić [30] found that the effective graph resistance of a connected network can be written as a function of all 
  # non-zero Laplacian eigenvalues of the network.
  
  # Yang et al. 2016. The Rationality of Four Metrics of Network Robustness: A Viewpoint of Robust Growth of 
  # Generalized Meshes. PLOS ONE. https://doi.org/10.1371/journal.pone.0161077
  
  lap.mat<-laplacian_matrix(undir.looped, norm=FALSE, sparse=FALSE)
  ev.lm<-eigen(lap.mat)
  ev.lm<-unlist(ev.lm[1])
  sorted.Laplacian.eigenvalues<-sort(ev.lm)
  value.ER<-0
  for(lam in 2:length(sorted.Laplacian.eigenvalues)){
    value.ER<-value.ER+1/sorted.Laplacian.eigenvalues[lam]
  }
  
  EffectiveGraphResistance<-value.ER*length(sorted.Laplacian.eigenvalues)
  
  if(is.na(EffectiveGraphResistance))
    EffectiveGraphResistance<-NA
  
  
  
  
  ########################################################################
  #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
  #6.1 Density
  
  #The proportion of present edges from all possible edges in the network.
  #Edge density: % of edges compared to maximum
  #densityscore<-ecount(showme.looped)/(vcount(showme.looped)*(vcount(showme.looped)-1)) #for a directed network
  densityscore<-ecount(showme.looped)/(vcount(showme.looped)^2)#for a directed network WITH self-loops
  
  
  #Mean degree
  deg <- degree(showme.looped,loops=FALSE)
  mean.degree<-mean(deg)
  
  #Mean path length
  mean_path_length<-mean_distance(showme.looped, directed=T)
  
  #network diameter is the longest geodesic distance (length of the shortest path between two nodes) in the network.
  network.diameter<-diameter(showme.looped, directed=T,weights=NA)
  
  #average clustering coefficient
  avg.cluster.coeff<-transitivity(showme.looped)
  
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
  #Modularity
  wtc<-cluster_walktrap(showme.looped)
  
  #for some reason the pseudo-warning ("Modularity is implemented for undirected graphs only.") is
  # damn-near unsuppressable---double team on it (invisible+capture.output!)
  options(warn=-1)
  invisible(capture.output(EB.boo<-cluster_edge_betweenness(showme.looped),type="message"))
  
  fg.comm<-cluster_fast_greedy(undir.looped)
  
  infomapped<-infomap.community(undir.looped)
  
  clp <-cluster_label_prop(undir.looped) # input graph should be undirected to make sense.
  options(warn=0)
  ##################
  wt.modularity.score<-modularity(wtc)
  wt.Nmodules<-length(unique(wtc$membership))
  
  EB.modularity.score<-modularity(EB.boo)
  EB.Nmodules<-length(unique(EB.boo$membership))
  
  fg.modularity.score<-modularity(fg.comm)
  fg.Nmodules<-length(unique(fg.comm$membership))
  
  im.modularity.score<-modularity(infomapped)
  im.Nmodules<-length(unique(infomapped$membership))
  
  clp.modularity.score<-modularity(clp)
  clp.Nmodules<-length(unique(clp$membership))
  
  
  ###############################
  printme<-data.frame(uniquebehaviors,uniquetransitions,
                      Shannon,trueShannon,propShannon,redundancy,
                      Shannon.T,trueShannon.T, propShannon.T,redundancy.T,
                      
                      densityscore,mean.degree,mean_path_length,network.diameter,avg.cluster.coeff,
                      PropSelfTransitioning,Smallworldness,
                      wt.modularity.score,wt.Nmodules,
                      EB.modularity.score,EB.Nmodules,
                      fg.modularity.score,fg.Nmodules,
                      im.modularity.score,im.Nmodules,
                      clp.modularity.score,clp.Nmodules,
                      EffectiveGraphResistance)
  
  
  #CONNECTIVITY measures reflect degree to which network differs from complete network
  #Edge density: % of edges compared to maximum
  #Average degree: Avg. number of links
  #Average path length: avg of shortest pathes between reachable nodes
  #Network diameter: longest of shortest paths
  
  #CENTRALITY measures quantify heterogeneity in network structure
  #Average clustering coefficient:
  #Components
  
  
  
  return(printme)
  
}





#Uses CurveDataTotal, and a complexity variable you define, to find maximal values for each window size for each video
CurveAccumulation<-function(CurveDataTotal){

  plotlist<-list()  
  plotlisted<-list()
  DFCurveDataTotal<-data.frame(CurveDataTotal)
  
  complexityvariables<-colnames(DFCurveDataTotal)
  complexityvariables<-complexityvariables[-c(1:2)]
  
  SixHundredList<-list()
  
  for (q in 1:length(complexityvariables)) {
    investigate<-complexityvariables[q]

    for(abc in c(na.omit(unique(DFCurveDataTotal[,1])))){ 
      
      currentvideo<-DFCurveDataTotal[which(DFCurveDataTotal$videoindex==abc),]
      
      
      for(efg in c(na.omit(unique(currentvideo$windowsizetouse)))){
        currentwindow<-currentvideo[which(currentvideo$windowsizetouse==efg),]
        

          currentwindowbehav <- currentwindow[!is.infinite(currentwindow[,investigate]),investigate]
          currentwindowbehav<-na.omit(currentwindowbehav)
          valuestoplot<-matrix(c(efg,max(currentwindowbehav)),nrow=1)
          colnames(valuestoplot)<-c("Window",investigate)
          
        if(efg==100){
          baseline<-matrix(c(0,0),nrow=1)
          videoplot<-rbind(baseline,valuestoplot)
        } else {
          videoplot<-rbind(videoplot,valuestoplot)
        }
      }
      
      videoplot<-na.omit(videoplot)
      
      
      noINF<-DFCurveDataTotal[!is.infinite((DFCurveDataTotal[,investigate])),]
      maxuse<-max(na.omit(noINF[,investigate]))
      
      
      videoplot.df<-data.frame(videoplot)

      jittered.Xs<-c(0,jitter(videoplot.df[-1,1],amount=3))
      
      #Only plot if there are at least 2 time windows to use (100 & 200)
      if(length(c(na.omit(unique(currentvideo$windowsizetouse))))>2){
        formula.plot <- as.formula(paste (investigate, "~ SSasymp(Window, Asym, R0, lrc)",sep=""))
        
  
        out <- tryCatch(nls( formula.plot , data = videoplot.df),error = function(e) { cat('In error handler\n'); print(e); e })#function checks if nls model returns error
              #If any of the values are "error", returns TRUE (and we skip)
        if(any(class(out) == "error")){
          #maxuse<-max(na.omit(DFCurveDataTotal[,investigate]))
          
          use.Xs<-jittered.Xs
          use.Ys<-videoplot.df[,2]
          line.Xs<-videoplot.df[,"Window"]
          line.Ys<-videoplot.df[,investigate]
          relevantcolor<-"dodgerblue4"
          linetype<-2
          
#           plot(jittered.Xs,videoplot.df[,2],xlim=c(0,1000),ylim=c(0,maxuse),col="dodgerblue4",xlab="Window",ylab=investigate)
#           lines(c(0,na.omit(unique(currentvideo$windowsizetouse))), videoplot.df[,investigate], lty = 2, col = "dodgerblue4")
        }  else {
        
            fm<-nls( formula.plot , data = videoplot.df)
            predictedYs<-predict(fm)
            
            use.Xs<-jittered.Xs
            use.Ys<-videoplot.df[,2]
            line.Xs<-videoplot.df[,1]
            line.Ys<-predictedYs
            relevantcolor<-"green4"
            linetype<-1
            
#             plot(jittered.Xs,videoplot.df[,2],xlim=c(0,1000),ylim=c(0,maxuse),xlab="Window",ylab=investigate,col="green4")
#             lines(videoplot.df[,1],predictedYs,lty=1,lwd=2,col="green4")
        }
        
        
        
        #         xCurve<-seq(0,1000,1)
        #         LO<-loess(predictedYs~jittered.Xs)
        #         predict(LO,xCurve)
        #         lines(xCurve,predict(LO,xCurve), lty = 1, col = "red",lwd=2)
      }
      if(length(c(na.omit(unique(currentvideo$windowsizetouse))))==2){
        #formula.plot <- as.formula(paste (investigate, "~ SSasymp(window, Asym, R0, lrc)",sep=""))
        #fm <- nls( formula.plot , data = videoplot.df)
        
        
        use.Xs<-videoplot.df[,1]
        use.Ys<-videoplot.df[,2]
        line.Xs<-c(0,na.omit(unique(currentvideo$windowsizetouse)))
        line.Ys<-videoplot.df[,investigate]
        relevantcolor<-"blue"        
        linetype<-2
        
#         plot(videoplot.df,xlim=c(0,1000),ylim=c(0,maxuse),col="blue",xlab="Window",ylab=investigate)
#         lines(c(0,na.omit(unique(currentvideo$windowsizetouse))), videoplot.df[,investigate], lty = 2, col = "blue")
      }  
      
      
      if(abc!=1){
      par(new=TRUE)
      xlab2=NA
      ylab2=NA
          
      } else {
      xlab2="window"
      ylab2=investigate
      }
      
      plotted<-{plot(use.Xs,use.Ys,xlim=c(0,1000),ylim=c(0,maxuse),xlab=xlab2,ylab=ylab2,col=relevantcolor,cex.axis = 1.5, cex.lab = 4); lines(line.Xs,line.Ys,lty=linetype,lwd=2,col=relevantcolor,cex.axis = 1.5, cex.lab = 2)}
      
      if(abc==1){
        Speciesplot<-videoplot
      } else {
        Speciesplot<-rbind(Speciesplot,videoplot)
      }
      
      
      
    }
  
  plotlist[[q]]<- plotted
  
  Speciesplot.df<-data.frame(Speciesplot)
  sixhunds<-Speciesplot.df[which(Speciesplot.df$Window==600 & Speciesplot.df[,2]!=-Inf),]
  
  SixHundredList[[q]]<-sixhunds
  
  }

plotsandsixhund<-list(plotlist,SixHundredList)  
  
return(plotsandsixhund)  
  
}


#Takes a .csv and splits it up into *splits* based on any "OFF" periods, then generates BehavioralSequences for analysis
SequenceCreator<-function(filename,location,species.name,specialclass,hb){
  
  #function for pulling right side of text string, for naming purposes    
  right = function (string, char){
    substr(string,nchar(string)-(char-1),nchar(string))
  }
  hangingbird<-hb
  ######Function for eliminating event behaviors during "OFF" and fixing duration behaviors
  #DURATION BEHAVIORS
  DurationBehaviors<-c("BP1","BP3","SS1","SS2","O3","O4","O5","OFF","O6","O7","O8","OPMH.D","OPMB.D","OPMF.D","OPMW.D","OPMT.D","RepWing","RepTail","RepHead","RepTorso",
                       "RM1","RM2","RM3","RM4","MO2","PU1","SP")
  DurationBehaviors.nosp<-c("BP1","BP3","SS1","SS2","O3","O4","O5","OFF","O6","O7","O8","OPMH.D","OPMB.D","OPMF.D","OPMW.D","OPMT.D","RepWing","RepTail","RepHead","RepTorso",
                            "RM1","RM2","RM3","RM4","MO2","PU1")
  DurationBehaviors.nosp.noOFF<-c("BP1","BP3","SS1","SS2","O3","O4","O5","O6","O7","O8","OPMH.D","OPMB.D","OPMF.D","OPMW.D","OPMT.D","RepWing","RepTail","RepHead","RepTorso",
                                  "RM1","RM2","RM3","RM4","MO2","PU1")
  
  EventBehaviors<-c("BP2","O1","O2","OPMH","OPMB","OPMF","OPMW","OPMT","OPAH","OPAC1","OPAC2","OPAC3","OPAW","OPATW","OPAMW","OPAMTT",
                    "OPABH","OPABB","OPABF","OPABW","OPABT","PC1","PC2","PC3","PC4","MO1","SP1","SP2","SP3","SP4")
  
  fixOFFs<-function(tabletofix){  
    offstarts<-nrow(subset(tabletofix,code=="OFF(start)"))
    offstops<-nrow(subset(tabletofix,code=="OFF(end)"))	
    tabletofix<-tabletofix[order(tabletofix$time),]
    if(offstarts>0){#If there are "OFFS", this if statement eliminates any behaviors measured during these OFFs and adds "STOPS"
      #At the beginning of the OFF period and adds "STARTS" at the end of the OFF period. Also deletes behaviors that happen while OFF.
      offstartings<-subset(tabletofix,code=="OFF(start)")[,1]
      offstoppins<-subset(tabletofix,code=="OFF(end)")[,1]
      
      for (vvv in 1:offstarts){#offstarts
        tabletofix<-tabletofix[order(tabletofix$time),]
        #stopBegin<-subset(tabletofix,code=="OFF(start)")[vvv,1]
        #stopEnd<-subset(tabletofix,code=="OFF(end)")[vvv,1]
        stopBegin<-offstartings[vvv]
        stopEnd<-offstoppins[vvv]
        
        events<- tabletofix[tabletofix$code %in% EventBehaviors,]
        gap.wise.misses<-events[which(events$time < stopEnd & events$time > stopBegin),]
        
        if(nrow(gap.wise.misses)>0){
          
          tabletofix<-tabletofix[-(which((tabletofix$time%in% gap.wise.misses$time) & (tabletofix$code%in% gap.wise.misses$code))),] #Removes rows of eventbehaviors which happen while birds is OFF screen
          #[-(which(tabletofix$time == oldend$time & tabletofix$code==oldend$code)),]
        }
        #testlist<-"O4"
        for (db in DurationBehaviors.nosp){  ##testlist #DurationBehaviors.nosp
          behav.obs <- tabletofix[grep(db, tabletofix$code), ]
          bouts<-nrow(behav.obs)/2
          if(bouts>0.6){
            starts<-behav.obs[grep("start", behav.obs$code), ]
            ends<-behav.obs[grep("end", behav.obs$code), ]
            
            oldend<-starts[1,][-1,]
            newend<-starts[1,][-1,]
            oldstart<-starts[1,][-1,]
            newstart<-starts[1,][-1,]
            
            if(nrow(starts)==nrow(ends)){ #only do this if starts/ends are equal for a given behavior
              for(nnn in 1:nrow(starts)){
                if (starts[nnn,1]==stopBegin & ends[nnn,1]<stopEnd & ends[nnn,2]!="OFF(end)"){ #S3_ If start occurs as OFF starts, and end occurs while OFF is active, delete
                  starttimetoeliminate<-starts[nnn,1]
                  stoptimetoeliminate<-ends[nnn,1]
                  startremove<-behav.obs[which(behav.obs$time == starttimetoeliminate),] 
                  endremove<-behav.obs[which(behav.obs$time == stoptimetoeliminate),] 
                  
                  gap.wise.misses2<-rbind(startremove,endremove)
                  
                  tabletofix<-tabletofix[-((match(gap.wise.misses2$time,tabletofix$time))),]
                }
                
                if(starts[nnn,1]<stopBegin & ends[nnn,1]==stopEnd & ends[nnn,2]!="OFF(end)" ){#S2_This part handles all duration behaviors that start before each (nnn) OFF starts and end as each OFF ends.
                  oldend<-ends[nnn,]
                  
                  
                  newend<-ends[1,]
                  newend$time<-stopBegin
                  
                  
                  tabletofix<-rbind(tabletofix,newend) #adds new end at start of OFF period
                  tabletofix<-tabletofix[-(which(tabletofix$time == oldend$time & tabletofix$code==oldend$code)),] #gets rid of old end behavior
                }  
                
                if(starts[nnn,1]<stopBegin & ends[nnn,1]<stopEnd & ends[nnn,1]>stopBegin & ends[nnn,2]!="OFF(end)"){#S1_This part handles all duration behaviors that start before each (nnn) OFF starts and end while OFF is on.
                  oldend<-ends[nnn,]
                  
                  
                  newend<-ends[1,]
                  newend$time<-stopBegin
                  
                  
                  tabletofix<-rbind(tabletofix,newend) #adds new end at start of OFF period
                  tabletofix<-tabletofix[-(which(tabletofix$time == oldend$time & tabletofix$code==oldend$code)),] #gets rid of old end behavior
                }  
                
                if(starts[nnn,1]>stopBegin & starts[nnn,1]<stopEnd & ends[nnn,1]==stopEnd & ends[nnn,2]!="OFF(end)" ){#E4_This part handles all duration behaviors that start while OFF and ends with ON
                  
                  oldstart<-starts[nnn,]
                  
                  newstart<-starts[1,]
                  newstart$time<-stopEnd
                  
                  tabletofix<-rbind(tabletofix,newstart) #
                  tabletofix<-tabletofix[-(which(tabletofix$time == oldstart$time & tabletofix$code==oldstart$code)),] #gets rid of old start behavior
                  
                }  
                
                
                if(starts[nnn,1]==stopBegin & ends[nnn,1]==stopEnd & ends[nnn,2]!="OFF(end)" ){#E3_This part handles all duration behaviors that start and end exactly as the OFF starts and ends
                  oldend<-ends[nnn,]
                  oldstart<-starts[nnn,]
                  
                  
                  tabletofix<-tabletofix[-(which(tabletofix$time == oldend$time & tabletofix$code==oldend$code)),] #gets rid of old end behavior
                  tabletofix<-tabletofix[-(which(tabletofix$time == oldstart$time & tabletofix$code==oldstart$code)),] #gets rid of old start behavior
                  
                }  
                
                if(starts[nnn,1]==stopBegin & ends[nnn,1]>stopEnd & ends[nnn,2]!="OFF(end)" ){#E2_This part handles all duration behaviors that start as OFF starts and end after each OFF ends.
                  oldstart<-starts[nnn,]
                  
                  
                  newstart<-starts[1,]
                  newstart$time<-stopEnd
                  
                  
                  tabletofix<-rbind(tabletofix,newstart) #adds new end at start of OFF period
                  tabletofix<-tabletofix[-(which(tabletofix$time == oldstart$time & tabletofix$code==oldstart$code)),] #gets rid of old end behavior
                } 
                
                if(starts[nnn,1]>stopBegin & starts[nnn,1]<stopEnd & ends[nnn,1]>stopEnd & ends[nnn,2]!="OFF(end)"){#E1_This part handles all duration behaviors that start while OFF is on, then after OFF is off
                  oldstart<-starts[nnn,]
                  
                  
                  newstart<-starts[1,]
                  newstart$time<-stopEnd
                  
                  
                  tabletofix<-rbind(tabletofix,newstart) #adds new end at start of OFF period
                  tabletofix<-tabletofix[-(which(tabletofix$time == oldstart$time & tabletofix$code==oldstart$code)),] #gets rid of old end behavior
                }  
                
                
                if (starts[nnn,1]<stopBegin & ends[nnn,1]>stopEnd & ends[nnn,2]!="OFF(end)"){ #SPAN_if start happens before bird goes OFF AND stop happens after bird comes back (OFF(end)), cut out middle and add new stop and start
                  starttimetoadd<-starts[nnn,]
                  starttimetoadd$time<-stopEnd
                  
                  stoptimetoadd<-ends[nnn,]
                  stoptimetoadd$time<-stopBegin
                  
                  
                  addstartstop<-rbind(starttimetoadd,stoptimetoadd)
                  
                  tabletofix<-rbind(tabletofix,addstartstop)
                }
                
                
                if (starts[nnn,1]>stopBegin & ends[nnn,1]<stopEnd & ends[nnn,2]!="OFF(end)"){ #MID_if start AND stop of a duration behavior occur while bird is OFF, delete
                  starttimetoeliminate<-starts[nnn,1]
                  stoptimetoeliminate<-ends[nnn,1]
                  startremove<-behav.obs[which(behav.obs$time == starttimetoeliminate),] 
                  endremove<-behav.obs[which(behav.obs$time == stoptimetoeliminate),] 
                  
                  gap.wise.misses2<-rbind(startremove,endremove)
                  
                  tabletofix<-tabletofix[-((match(gap.wise.misses2$time,tabletofix$time))),]
                }
                
                
              }# for every row in 'starts', which has all starts for a behavior, it calculates the difference between end and start and then adds this to 'duration'
            }    
          }
        }
        
      }
      
    } 
    
    #REMOVES behaviors that have a start and stop at the same time
    
    
    for (db2 in DurationBehaviors.nosp){  
      behav.observations <- tabletofix[grep(db2, tabletofix$code), ]
      startings<-behav.observations[grep("start",behav.observations$code),]
      endings<-behav.observations[grep("end",behav.observations$code),]
      
      startstoremove<-startings[(startings$time %in% endings$time),]
      endstoremove<-endings[(endings$time %in% startings$time),]
      kill<-rbind(startstoremove,endstoremove)
      tabletofix<-tabletofix[!(tabletofix$time %in% kill$time & tabletofix$code %in% kill$code),]
      
    }	
    
    tabletofix[order(tabletofix$time),]->table
    
    return(table)
  }  
 
  # setting wd, and performing initial naming tasks related to file output 
  ############################
  dir<-getwd()
  saveithere<-paste(dir,"/checks/",sep="")
  FileisNamed<-right(filename,46)
  
  name<-paste(FileisNamed,"pdf",sep=".")
  name<-gsub("/", ".", name, fixed = TRUE)
  name<-gsub(".csv", "", name, fixed = TRUE)
  name<-gsub("000000_000", "", name, fixed = TRUE)
  name<-gsub("_2000-01-01T","",name, fixed = TRUE)
  name<-gsub("ut_oneatatime.","",name, fixed = TRUE)
  name<-gsub("eatatime.","",name, fixed = TRUE)
  name<-gsub("sis.FixedCSVs.","",name, fixed = TRUE)
  
  
  name2<-gsub("/", ".", FileisNamed, fixed = TRUE)
  name2<-gsub(".csv", "", name2, fixed = TRUE)
  name2<-gsub("000000_000", "", name2, fixed = TRUE)
  name2<-gsub("_2000-01-01T","",name2, fixed = TRUE)
  name2<-gsub("ut_oneatatime.","",name2, fixed = TRUE)
  name2<-gsub("eatatime.","",name2, fixed = TRUE)
  name2<-gsub("sis.FixedCSVs.","",name2, fixed = TRUE)
  
  
  fullname<-paste(location,name,sep="/")
  ###########  
  #read in table
  read.table(filename,sep=",",header=TRUE)->table
  
  #Order selections according to begin time
  table[order(table$time),]->table
  
  
  
  #FIRST STEP PREPROCESSING
  
  #Early files used "Mpres" instead of "M pres"...these were incorrect and are removed in this step		
  rowswithMpresnospace<-which(table$code == "Mpres")
  if (length(rowswithMpresnospace)!=0){
    table<-table[-((which(table$code == "Mpres"))),] #Removes 'Mpres'
  }
  
  ####Remove spaces from 'code' identifiers
  table$code<-gsub(" ", "", table$code, fixed = TRUE)
  
  #Order selections according to begin time	
  table[order(table$time),]->table	
  
  #ELIMINATES DUPLICATE ROWS	
  table<-unique(table) 
  ###########
  #DATA CHECKING LOCATION, BEHAVIOR BY BEHAVIOR
  behaviortoinvestigate<-"OPMF"
  info<-table[grep(behaviortoinvestigate, table$code), ]
  rle(table[grep(behaviortoinvestigate, table$code), ][,2])
  ######################
  
  #Make repetitively pressed buttons (e.g. OPMW, OPAW, OPAC3) into duration of RM behaviors (e.g. RM1)
  
  flappinglist<-c("OPAW","OPAC3","OPAMW")	
  wingsflapping<-subset(table,code %in% flappinglist)	
  wingsflapping<-subset(wingsflapping,!duplicated(time))
  if(nrow(wingsflapping)>0){
    wingsflapping$timetonext<-c(diff(wingsflapping$time),1000)
    wingsflapping$keep<- ifelse((wingsflapping$code==(shift(wingsflapping$code,1L))) | (wingsflapping$code!=(shift(wingsflapping$code,1L)) & wingsflapping$timetonext>0.2),"yes","no")
    wingsflapping<-subset(wingsflapping,keep=="yes")
    wingsflapping<-wingsflapping[,c(1:3)]
  }
  
  
  taillist<-c("OPMT")
  tailflapping<-subset(table,code %in% taillist)	
  
  headlist<-c("OPAH")
  headflapping<-subset(table,code %in% headlist)
  
  torsolist<-c("OPAC1","OPAC2","OPATW","OPAMTT")
  torsoflapping<-subset(table, code %in% torsolist)
  torsoflapping<-subset(torsoflapping,!duplicated(time))
  if(nrow(torsoflapping)>0){
    torsoflapping$timetonext<-c(diff(torsoflapping$time),1000)
    torsoflapping$keep<- ifelse((torsoflapping$code==(shift(torsoflapping$code,1L))) | (torsoflapping$code!=(shift(torsoflapping$code,1L)) & torsoflapping$timetonext>0.2),"yes","no")
    torsoflapping<-subset(torsoflapping,keep=="yes")
    torsoflapping<-torsoflapping[,c(1:3)]
  }
  
  
  
  repeptiveornamentaccentuations<-list(wingsflapping,tailflapping,headflapping,torsoflapping)	
  
  durationofrepetbehavs<-table[1,]	
  durationofrepetbehavs<-durationofrepetbehavs[-1,]
  RepBehavLabels<-c("RepWing","RepTail","RepHead","RepTorso")
  
  for (RMs in 1:length(repeptiveornamentaccentuations)){	
    #flappers<-tailflapping
    
    newlabel<-RepBehavLabels[RMs]
    flappers<-repeptiveornamentaccentuations[[RMs]]
    flappers<-flappers[order(flappers$time),]
    
    if(nrow(flappers)>1){
      attach(flappers)
      
      timetonext<-c(diff(time),"FALSE")	
      TFseq<-timetonext<1 
      NewVec<-vector(length=length(TFseq))
      
      for (gg in 1:length(TFseq)) {
        
        xnow<-TFseq[gg]
        pre<-ifelse(gg>1,TFseq[(gg-1)],NA)
        nex<-ifelse(gg==length(TFseq),NA,TFseq[(gg+1)])
        
        if(is.na(pre)){  
          if(xnow==TRUE & is.na(pre)){
            newvalue<-"start"
          } else {
            if (xnow==FALSE & is.na(pre)){
              newvalue<-"off"
            }
          }
        } else {
          if (is.na(nex)){
            if (xnow==FALSE & (pre)==FALSE){
              newvalue<-"off"
            }
            if (xnow==FALSE & (pre)==TRUE){
              newvalue<-"end"
            }
          } else {
            if(xnow==TRUE & (nex)==TRUE & (pre)==FALSE){
              newvalue<-"start"
            } 
            if (xnow==TRUE & (nex)==FALSE & (pre)==TRUE){
              newvalue<-"on"
            }
            if(xnow==TRUE & (nex)==FALSE & (pre)==FALSE){
              newvalue<-"start"
            } 
            if (xnow==TRUE & (nex)==TRUE & (pre)==TRUE){
              newvalue<-"on"
            }
            if (xnow==FALSE & (nex)==TRUE & (pre)==TRUE){
              newvalue<-"end"
            }
            if (xnow==FALSE & (nex)==TRUE & (pre)==FALSE){
              newvalue<-"off"
            }
            if (xnow==FALSE & (nex)==FALSE & (pre)==FALSE){
              newvalue<-"off"
            }
            if (xnow==FALSE & (nex)==FALSE & (pre)==TRUE){
              newvalue<-"end"
            }
          }
        }
        NewVec[gg]<-newvalue
      }
      detach(flappers)
      
      flappers$times<-NewVec
      startsofrep<-subset(flappers,times=="start")
      if(nrow(startsofrep)>0){startsofrep$code<-paste(newlabel,"(start)",sep="")}
      endsofrep<-subset(flappers,times=="end")
      if(nrow(endsofrep)>0){endsofrep$code<-paste(newlabel,"(end)",sep="")}
      
      thesebehaviors<-rbind(startsofrep,endsofrep)
      thesebehaviors<-thesebehaviors[,-4]
      
      
      durationofrepetbehavs<-rbind(durationofrepetbehavs,thesebehaviors)
      durationofrepetbehavs<-durationofrepetbehavs[order(durationofrepetbehavs$time),]	
    }
  }	
  
  table<-rbind(table,durationofrepetbehavs)	
  
  
  #SECOND STEP PROCESSING (handling OFF periods)	
  
  
  #############################################################	
  #if last OFF (end) is the same as the END, this cuts off analysis at that point
  
  if(nrow(subset(table,code=="OFF(end)"))>0){
    if(isTRUE(all.equal(max((subset(table,code=="OFF(end)")[,1]),na.rm=TRUE),subset(table,code=="END")[,1],tolerance = 0.0009))){	  
      OriginalEndingTime<-	max((subset(table,code=="OFF(end)")[,1]),na.rm=TRUE)
      TimeBirdFlewOff<-tail(subset(table,code=="OFF(start)")[,1],1)
      FinishEnd<-data.frame("time"=TimeBirdFlewOff,"code"="END","class"=0)
      beforeflewaway<-subset(table,time<=TimeBirdFlewOff & code !="OFF(start)")
      previousOFFs<-subset(table,time<TimeBirdFlewOff & code =="OFF(start)")
      
      AtEnd<-subset(table,time==OriginalEndingTime & code !="END" & code != "OFF(end)")
      if(nrow(AtEnd)>0){
        AtEnd$time<-TimeBirdFlewOff
      }
      
      
      tabletable<-rbind(beforeflewaway,previousOFFs,AtEnd,FinishEnd)
      table<-tabletable[order(tabletable$time),]
    }
  }	
  
  #Finds END time and uses to calculate TotalTime
  TotalTime<-table[which(table$code == "END"), ][,1]	
  
  ######if first OFF (start) is at t=0, this cuts off analysis during this first 'missing' bit
  BirdWasOffatBeginning<-0
  TimeBirdFlewON<-0
  if(nrow(subset(table,code=="OFF(start)"))>0){
    if(isTRUE(all.equal(min((subset(table,code=="OFF(start)")[,1]),na.rm=TRUE),0,tolerance = 0.0009))){	  
      BirdWasOffatBeginning<-1
      OriginalStartingTime<-	0
      TimeBirdFlewON<-subset(table,code=="OFF(end)")[1,1]
      #FinishEnd<-data.frame("time"=TimeBirdFlewOff,"code"="END","class"=0)
      afterflewON<-subset(table,time>=TimeBirdFlewON & code !="OFF(end)")
      subsequentOFFs<-subset(table,time>TimeBirdFlewON & code =="OFF(end)")
      
      AtStart<-subset(table,time==0 &  code != "OFF(start)")
      if(nrow(AtStart)>0){
        AtStart$time<-TimeBirdFlewON
      }
      
      
      tabletabletable<-rbind(afterflewON,subsequentOFFs,AtStart)
      table<-tabletabletable[order(tabletabletable$time),]
      
      
    }
  }	
  
  
  OFFtimes<-table[grep("OFF", table$code), ][,1]
  OFFstartTimes<-subset(table,code=="OFF(start)")[,1]  
  OFFendTimes<-subset(table,code=="OFF(end)")[,1]  
  
  #############################################################	
  
  #PLOT A	
  
  
  table<-table[order(table$time),]	
  table<-fixOFFs(table)	
  
  #Adds O2 for every changed in O3
  O2additions<-rbind(table[which(table$code == "O3(start)"),],table[which(table$code == "O3(end)"),])
  O2additions$code<-(gsub("O3(end)", "O2", gsub("O3(start)", "O2", O2additions$code, fixed = TRUE), fixed = TRUE)) #substituting
  table<-rbind(table,O2additions)	
  
  
  #THIRD STEP PROCESSING (handling body orientation issues, e.g. O5s and O6s, O4 contingencies)	
  #############################################################			
  #############################################################	
  #############################################################	
  #Finds END time and uses to calculate TotalTime
  TotalTime<-table[which(table$code == "END"), ][,1]
  
  
  #If there are no 05s, then the bird was upright the entire time
  if(nrow(table[which(table$code == "O5(start)"), ])!=0){ 
    firstO5<-min(table[which(table$code == "O5(start)"), ][,1],na.rm = TRUE) #Takes times of all O5 starts and finds first (minimum)
    lastO5<-max(table[which(table$code == "O5(end)"), ][,1],na.rm = TRUE) #Takes times of all O5 starts and finds last (maximum)
    # O5starttable<-rbind(table[which(table$code == "O5(end)"),],table[which(table$code == "O5(start)"),],table[which(table$code == "END"),])#binds all O5 starts, ends, and the Trial End
    O5starttable<-table[grep("O5", table$code), ]
    endtime<-table[which(table$code == "END"),][,1]
    
    O6starttable<-O5starttable
    # O6starttable<-O6starttable[-(which(O6starttable$code=="O5(start)" & O6starttable$time %in% OFFendTimes)),]
    # O6starttable<-O6starttable[-(which(O6starttable$code=="O5(end)" & O6starttable$time %in% OFFstartTimes)),]
    
    if (BirdWasOffatBeginning>0){
      if (isTRUE(all.equal(firstO5,TimeBirdFlewON,tolerance = 0.0009))){
        if (isTRUE(all.equal(lastO5,endtime,tolerance = 0.0009))) { #If O5 begins and ends a trial, substiutions should not be made on those obs
          O6starttable<-O6starttable[!(O6starttable$code == "O5(start)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, TimeBirdFlewON, tolerance = 0.0009)))), ]
          O6starttable<-O6starttable[!(O6starttable$code == "O5(end)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]
          O6starttable$code<-(gsub("O5(end)", "O6(start)", O6starttable$code, fixed = TRUE)) 
          O6starttable$code<-(gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE))#substituting O6starts for O5 ends, etc
          
          O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
          
          
          
        } else { #If O5 begins a trial, substitutions should not be made on that first obs, and O6(end) should be at the end
          O6starttable<-O6starttable[!(O6starttable$code == "O5(start)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, TimeBirdFlewON, tolerance = 0.0009)))), ]
          O6starttable$code<-(gsub("O5(end)", "O6(start)", O6starttable$code, fixed = TRUE)) 
          O6starttable$code<-(gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE))#substituting O6starts for O5 ends, etc
          O6starttable$code<-(gsub("END", "O6(end)", O6starttable$code, fixed = TRUE)) #substituting O6starts for O5 ends, etc
          
          O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
          
          O6ends<-O6starttable[1,]
          O6ends$time<-endtime
          O6ends$code<-"O6(end)"
          
          O6starttable<-rbind(O6starttable,O6ends)
          
        }
      } else { #IF O5 doesn't start the trial, O6 must
        if (isTRUE(all.equal(lastO5,endtime,tolerance = 0.0009))) {
          O6startsatzero<-O6starttable[1,]
          O6startsatzero$time<-TimeBirdFlewON
          O6startsatzero$code<-"O6(start)" 
          O6starttable<-rbind(O6startsatzero,O6starttable)#This (+3 lines above) adds an O6(start) at t=0
          
          O6starttable<-O6starttable[!(O6starttable$code == "O5(end)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]#Don't replace O5(end) if it's at the END
          O6starttable$code<-(gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE))
          O6starttable$code<-(gsub("O5(end)", "O6(start)", O6starttable$code, fixed = TRUE))
          
          O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
          
          
        } else { #When O5 neither starts the trial nor ends the trial, all subs can be made, END can be turned into O6(end), and O6(start) can be added at t=0
          O6startsatzero<-O6starttable[1,]
          O6startsatzero$time<-TimeBirdFlewON
          O6startsatzero$code<-"O6(start)" 
          O6starttable<-rbind(O6startsatzero,O6starttable)#This (+3 lines above) adds an O6(start) at t=0
          
          O6starttable$code<-(gsub("O5(end)", "O6(start)", gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE), fixed = TRUE)) #substituting O6starts for O5 ends, etc
          O6starttable$code<-(gsub("END", "O6(end)", O6starttable$code, fixed = TRUE)) #If 05 doesn't end the trial, O6 must
          
          O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
          
          O6ends<-O6starttable[1,]
          O6ends$time<-endtime
          O6ends$code<-"O6(end)"
          
          O6starttable<-rbind(O6starttable,O6ends)
          
          
        }
      }
    } else {
      if (isTRUE(all.equal(firstO5,0,tolerance = 0.0009))){
        if (isTRUE(all.equal(lastO5,endtime,tolerance = 0.0009))) { #If O5 begins and ends a trial, substiutions should not be made on those obs
          O6starttable<-O6starttable[!(O6starttable$code == "O5(start)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, 0, tolerance = 0.0009)))), ]
          O6starttable<-O6starttable[!(O6starttable$code == "O5(end)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]
          O6starttable$code<-(gsub("O5(end)", "O6(start)", O6starttable$code, fixed = TRUE)) 
          O6starttable$code<-(gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE))#substituting O6starts for O5 ends, etc
          
          O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
          
          
        } else { #If O5 begins a trial, substitutions should not be made on that first obs, and O6(end) should be at the end
          O6starttable<-O6starttable[!(O6starttable$code == "O5(start)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, 0, tolerance = 0.0009)))), ]
          O6starttable$code<-(gsub("O5(end)", "O6(start)", O6starttable$code, fixed = TRUE)) 
          O6starttable$code<-(gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE))#substituting O6starts for O5 ends, etc
          O6starttable$code<-(gsub("END", "O6(end)", O6starttable$code, fixed = TRUE)) #substituting O6starts for O5 ends, etc
          
          O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
          
          O6ends<-O6starttable[1,]
          O6ends$time<-endtime
          O6ends$code<-"O6(end)"
          
          O6starttable<-rbind(O6starttable,O6ends)
          
        }
      } else { #IF O5 doesn't start the trial, O6 must
        if (isTRUE(all.equal(lastO5,endtime,tolerance = 0.0009))) {
          O6startsatzero<-O6starttable[1,]
          O6startsatzero$time<-0
          O6startsatzero$code<-"O6(start)" 
          O6starttable<-rbind(O6startsatzero,O6starttable)#This (+3 lines above) adds an O6(start) at t=0
          
          O6starttable<-O6starttable[!(O6starttable$code == "O5(end)" & sapply(O6starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]#Don't replace O5(end) if it's at the END
          O6starttable$code<-(gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE))
          O6starttable$code<-(gsub("O5(end)", "O6(start)", O6starttable$code, fixed = TRUE))
          
          O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
          
          
        } else { #When O5 neither starts the trial nor ends the trial, all subs can be made, END can be turned into O6(end), and O6(start) can be added at t=0
          O6startsatzero<-O6starttable[1,]
          O6startsatzero$time<-0
          O6startsatzero$code<-"O6(start)" 
          O6starttable<-rbind(O6startsatzero,O6starttable)
          #This (+3 lines above) adds an O6(start) at t=0
          
          O6ends<-O6starttable[1,]
          O6ends$time<-endtime
          O6ends$code<-"O6(end)"
          
          O6starttable$code<-(gsub("O5(end)", "O6(start)", gsub("O5(start)", "O6(end)", O6starttable$code, fixed = TRUE), fixed = TRUE)) #substituting O6starts for O5 ends, etc
          O6starttable$code<-(gsub("END", "O6(end)", O6starttable$code, fixed = TRUE)) #If 05 doesn't end the trial, O6 must
          
          O6starttable$class<-(gsub("0", "5", O6starttable$class, fixed = TRUE))
          
          O6starttable<-rbind(O6starttable,O6ends)
          
          
        }
      }
    }
  } else {
    if(BirdWasOffatBeginning>0){
      O6starttable<-data.frame("time"=c(TimeBirdFlewON,TotalTime),"code"=c("O6(start)","O6(end)"),"class"=c(5,5))
    } else {
      O6starttable<-data.frame("time"=c(0,TotalTime),"code"=c("O6(start)","O6(end)"),"class"=c(5,5))
    }
  } 
  
  table2<-rbind(table,O6starttable)
  
  
  table<-table2 #keeps table our working, modified dataframe
  table<-unique(table) #ELIMINATES DUPLICATE ROWS
  table<-table[order(table$time),]	
  table<-fixOFFs(table)	#Fixes OFF issues
  
  #PLOT D
  #PLOTB<-table
  
  
  #Order selections according to begin time
  table<-table[order(table$time),]
  
  O4active<-rbind(table[which(table$code == "O4(start)"),],table[which(table$code == "O4(end)"),])
  O4active<-O4active[order(O4active$time),]#Order selections according to begin time
  O4starts<-table[which(table$code == "O4(start)"),]
  O4ends<-table[which(table$code == "O4(end)"),]
  O4bouts<-length(table[which(table$code == "O4(start)"),][,1])
  
  
  compareO5<-table[which(table$code == "O5(start)"),]
  compareO5end<-table[which(table$code == "O5(end)"),]
  compareO6<-table[which(table$code == "O6(start)"),]
  compareO6end<-table[which(table$code == "O6(end)"),]
  
  #This if addresses all postural (O5, O6,O7,O8) issues (if there are at least one or more O4s)
  if(nrow(O4starts)==nrow(O4ends) & O4bouts>0){   #IF BRANCH IS NOT VERTICAL THE ENTIRE TIME, WE NEED THIS CONTINGENCY (ELSE FOR CASES WHEN WE NEVER HAVE ANY O4 PERIODS/BOUTS)
    
    #
    if(nrow(compareO5)==nrow(compareO5end) & nrow(compareO5)!=0 & nrow(compareO5end)!=0) { #ONLY DO THIS IF EQUAL STARTS/STOPS FOR O5 (and !=0)
      for(j in 1:length(compareO5[,1])){
        
        for(i in 1:O4bouts){
          start<-O4starts[i,1]
          stop<-O4ends[i,1]
          startrow<-O4starts[i,]
          stoprow<-O4ends[i,]
          
          if(compareO5$time[j]==start & compareO5end$time[j]<stop ){#S3_This part handles O5s that start exactly as O4 starts and end before each O4 ends.
            oldendO5<-compareO5end[j,]
            oldstartO5<-compareO5[j,]
            
            O7start<-oldstartO5
            O7start$code="O7(start)"
            O7end<-oldendO5
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O7O8s<-newO7s
            O7O8s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
            
          }  
          
          if(compareO5$time[j]<start & compareO5end$time[j]==stop ){#S2_This part handles O5s that start before each O4 starts and end as each O4 ends.
            oldendO5<-compareO5end[j,]
            oldstartO5<-compareO5[j,]
            
            O8start<-oldstartO5
            O8start$code="O8(start)"
            O8end<-startrow
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O7start<-startrow
            O7start$code="O7(start)"
            O7end<-stoprow
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O7O8s<-rbind(newO7s,newO8s)
            O7O8s$class<-paste(i,j)
            table<-rbind(table,newO7s) #adds new O7s and O8s
            
          }  
          
          if(compareO5$time[j]<start & compareO5end$time[j]<stop & compareO5end$time[j]>start ){#S1_This part handles all O5s that start before each O4 starts and end while O4 is on
            oldendO5<-compareO5end[j,]
            oldstartO5<-compareO5[j,]
            
            O8start<-oldstartO5
            O8start$code="O8(start)"
            O8end<-startrow
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O7start<-startrow
            O7start$code="O7(start)"
            O7end<-oldendO5
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O7O8s<-rbind(newO7s,newO8s)
            O7O8s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
          }  
          
          if(compareO5$time[j]>start & compareO5$time[j]<stop & compareO5end$time[j]==stop ){#E4_This part handles all O5s that start while O4 and end exactly as O4 ends
            oldendO5<-compareO5end[j,]
            oldstartO5<-compareO5[j,]
            
            Nope<-oldstartO5[-1,]
            newO8s<-Nope
            
            O7start<-oldstartO5
            O7start$code="O7(start)"
            O7end<-oldendO5
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O7O8s<-rbind(newO7s,newO8s)
            O7O8s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
            
          }  
          
          if(compareO5$time[j]==start & compareO5end$time[j]==stop ){#E3_This part handles all O5s that start and end exactly as O4 starts and ends
            oldendO5<-compareO5end[j,]
            oldstartO5<-compareO5[j,]
            
            Nope<-oldstartO5[-1,]
            newO8s<-Nope
            
            O7start<-oldstartO5
            O7start$code="O7(start)"
            O7end<-oldendO5
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O7O8s<-rbind(newO7s,newO8s)
            O7O8s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
            
          }  
          
          if(compareO5$time[j]==start & compareO5end$time[j]>stop ){#E2_This part handles O5s that start as O4 starts and end after O4 ends.
            oldendO5<-compareO5end[j,]
            oldstartO5<-compareO5[j,]
            
            O8start<-stoprow
            O8start$code="O8(start)"
            O8end<-oldendO5
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O7start<-oldstartO5
            O7start$code="O7(start)"
            O7end<-stoprow
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O7O8s<-rbind(newO7s,newO8s)
            O7O8s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
          } 
          
          if(compareO5$time[j]>start & compareO5$time[j]<stop & compareO5end$time[j]>stop ){#E1_This part handles O5s that start while O4 is on, then end after O4 is off
            oldendO5<-compareO5end[j,]
            oldstartO5<-compareO5[j,]
            
            O8start<-stoprow
            O8start$code="O8(start)"
            O8end<-oldendO5
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O7start<-oldstartO5
            O7start$code="O7(start)"
            O7end<-stoprow
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O7O8s<-rbind(newO7s,newO8s)
            O7O8s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
            
          }  
          
          if (compareO5$time[j]<start & compareO5end$time[j]>stop){ #SPAN_if start happens before bird goes OFF AND stop happens after bird comes back (OFF(end)), cut out middle and add new stop and start
            oldendO5<-compareO5end[j,]
            oldstartO5<-compareO5[j,]
            
            O8start<-rbind(oldstartO5,stoprow)
            O8start$code="O8(start)"
            O8end<-rbind(startrow,oldendO5)
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O7start<-startrow
            O7start$code="O7(start)"
            O7end<-stoprow
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O7O8s<-rbind(newO7s,newO8s)
            O7O8s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
            
          }
          
          if (compareO5$time[j]>start & compareO5end$time[j]<stop){ #MID_if start AND stop of a duration behavior occur while bird is OFF, delete
            oldendO5<-compareO5end[j,]
            oldstartO5<-compareO5[j,]
            
            O7start<-oldstartO5
            O7start$code="O7(start)"
            O7end<-oldendO5
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O7O8s<-newO7s
            O7O8s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
          }
          
          if ((compareO5$time[j]<start & compareO5end$time[j]<start) | (compareO5$time[j]>stop & compareO5end$time[j]>stop)){ #PRE OR POST_if start AND stop of a duration behavior occur while bird is not on O4
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            O7start<-rbind(oldstartO6,oldendO6)
            O7start$code=c("O7(start)","O7(end)")
            newO7s<-O7start
            
            table<-rbind(table,newO7s)
          }   
          
        }
      }
    }
    if(nrow(compareO6)==nrow(compareO6end) & nrow(compareO6)!=0 & nrow(compareO6end)!=0) { #ONLY DO THIS IF EQUAL STARTS/STOPS FOR O6 (and !=0)
      for(j in 1:length(compareO6[,1])){
        
        for(i in 1:O4bouts){
          start<-O4starts[i,1]
          stop<-O4ends[i,1]
          startrow<-O4starts[i,]
          stoprow<-O4ends[i,]
          
          
          
          if(compareO6$time[j]==start & compareO6end$time[j]<stop ){#S3_This part handles O6s that start exactly as O4 starts and end before each O4 ends.
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            
            
            
          }  
          if(compareO6$time[j]<start & compareO6end$time[j]==stop ){#S2_This part handles O6s that start before each O4 starts and end as each O4 ends.
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            
            O7start<-oldstartO6
            O7start$code="O7(start)"
            O7end<-startrow
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O8start<-startrow
            O8start$code="O8(start)"
            O8end<-stoprow
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O8O7s<-rbind(newO8s,newO7s)
            O8O7s$class<-paste(i,j)
            table<-rbind(table,newO7s) #adds new O8s and O7s
            
          }  
          
          if(compareO6$time[j]<start & compareO6end$time[j]<stop & compareO6end$time[j]>start ){#S1_This part handles all O6s that start before each O4 starts and end while O4 is on
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            
            O7start<-oldstartO6
            O7start$code="O7(start)"
            O7end<-startrow
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O8start<-startrow
            O8start$code="O8(start)"
            O8end<-oldendO6
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O8O7s<-rbind(newO8s,newO7s)
            O8O7s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
          }  
          
          if(compareO6$time[j]>start & compareO6$time[j]<stop & compareO6end$time[j]==stop ){#E4_This part handles all O6s that start while O4 and end exactly as O4 ends
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            
            
            
          }  
          
          if(compareO6$time[j]==start & compareO6end$time[j]==stop ){#E3_This part handles all O6s that start and end exactly as O4 starts and ends
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            
            
            
          }  
          
          if(compareO6$time[j]==start & compareO6end$time[j]>stop ){#E2_This part handles O6s that start as O4 starts and end after O4 ends.
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            
            O7start<-stoprow
            O7start$code="O7(start)"
            O7end<-oldendO6
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O8start<-oldstartO6
            O8start$code="O8(start)"
            O8end<-stoprow
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O8O7s<-rbind(newO8s,newO7s)
            O8O7s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
          } 
          
          if(compareO6$time[j]>start & compareO6$time[j]<stop & compareO6end$time[j]>stop ){#E1_This part handles O6s that start while O4 is on, then end after O4 is off
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            
            O7start<-stoprow
            O7start$code="O7(start)"
            O7end<-oldendO6
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O8start<-oldstartO6
            O8start$code="O8(start)"
            O8end<-stoprow
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O8O7s<-rbind(newO8s,newO7s)
            O8O7s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
            
          }  
          
          if (compareO6$time[j]<start & compareO6end$time[j]>stop){ #SPAN_if start happens before bird goes O4 AND stop happens after O4 goes OFF, cut out middle and add new stop and start
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            
            O7start<-rbind(oldstartO6,stoprow)
            O7start$code="O7(start)"
            O7end<-rbind(startrow,oldendO6)
            O7end$code<-"O7(end)"
            newO7s<-rbind(O7start,O7end)
            newO7s$class<-paste(i,j)
            
            O8start<-startrow
            O8start$code="O8(start)"
            O8end<-stoprow
            O8end$code<-"O8(end)"
            newO8s<-rbind(O8start,O8end)
            
            O8O7s<-rbind(newO8s,newO7s)
            O8O7s$class<-paste(i,j)
            table<-rbind(table,newO7s) 
            
          }
          
          if (compareO6$time[j]>start & compareO6end$time[j]<stop){ #MID_if start AND stop of a duration behavior occur while bird is OFF, delete
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            
            
          }
          
          if ((compareO6$time[j]<start & compareO6end$time[j]<start) | (compareO6$time[j]>stop & compareO6end$time[j]>stop)){ #PRE OR POST_if start AND stop of O6 occurs completely before or after an O4 bout
            oldendO6<-compareO6end[j,]
            oldstartO6<-compareO6[j,]
            O7start<-rbind(oldstartO6,oldendO6)
            O7start$code=c("O7(start)","O7(end)")
            newO7s<-O7start
            
            table<-rbind(table,newO7s)
          }
        }
      }
    } 
    
  }	    
  ################################################
  #This if addresses all postural (O5, O6,O7,O8) issues (if there are no O4s)
  if(nrow(O4starts)==nrow(O4ends) & O4bouts==0){
    oldendO5<-compareO5end
    oldstartO5<-compareO5
    
    oldendO6<-compareO6end
    oldstartO6<-compareO6
    
    if(nrow(oldstartO5)>0){
      O8start<-oldstartO5
      O8start$code="O8(start)"
      O8end<-oldendO5
      O8end$code<-"O8(end)"
      newO8s<-rbind(O8start,O8end)
      newO8s$class<-5
    } else {
      newO8s<-oldendO5
      if(nrow(newO8s)>0){
        newO8s$class<-5
      }
    }
    
    if(nrow(oldstartO6)>0){
      O7start<-oldstartO6
      O7start$code="O7(start)"
      O7end<-compareO6end
      O7end$code<-"O7(end)"
      newO7s<-rbind(O7start,O7end)
      newO7s$class<-5
    } else {
      newO7s<-oldstartO6
      if(nrow(newO7s)>0){
        newO7s$class<-5
      }
      
    }
    
    O7O8s<-rbind(newO7s,newO8s)
    O7O8s$class<-5
    table<-rbind(table,O7O8s) 
  }
  
  midcheck<-table
  table<-unique(table) #ELIMINATES DUPLICATE ROWS
  table<-table[order(table$time),]	
  
  
  #ORDER (start/stop) FIXER FOR DURATION BEHAVIOR
  sevens<-table[grep("O7",table$code),] #O7 values
  if(length(sevens[,1])>0) {sevens$class<-"77"} #Standardize O7 class values
  sevens<-unique(sevens)
  
  
  sevens$keep<- ifelse((sevens$code==(data.table::shift(sevens$code,n=1L,type="lead"))) ,"no","yes")
  sevens<-subset(sevens,keep=="yes" | is.na(keep))[,-4]
  sevens<-subset(sevens,!(is.na(time)))
  
  table<-table[!(grepl("O7", table$code)),] #get rid of O7 valuse
  
  
  
  
  
  
  
  
  
  table<-rbind(table,sevens)
  
  
  
  
  table<-fixOFFs(table)	#Fixes OFF issues
  
  # 
  # 	AllO7<-table[which(table$code=="O7(start)" | table$code=="O7(end)" ),]
  # 	AllO7<-AllO7[order(AllO7$time),]	
  # 	
  # 	AllO8<-table[which(table$code=="O8(start)" | table$code=="O8(end)"),]
  # 	AllO8<-AllO8[order(AllO8$time),]	
  # 	
  # 	noO7O8<-table[-which(table$code=="O7(start)" | table$code=="O7(end)" | table$code=="O8(start)" | table$code=="O8(end)"),]
  # 	table<-rbind(noO7O8,AllO7,AllO8)
  #########################################################################	
  #If there are no O7s, then the bird was O8 the entire time
  TTolerance <- 0.0009
  
  
  if(nrow(table[which(table$code == "O7(start)"), ])!=0){ 
    firstO7<-min(table[which(table$code == "O7(start)"), ][,1],na.rm = TRUE) #Takes times of all O7 starts and finds first (minimum)
    lastO7<-max(table[which(table$code == "O7(end)"), ][,1],na.rm = TRUE) #Takes times of all O7 starts and finds last (maximum)
    # O7starttable<-rbind(table[which(table$code == "O7(end)"),],table[which(table$code == "O7(start)"),],table[which(table$code == "END"),])#binds all O7 starts, ends, and the Trial End
    O7starttable<-table[grep("O7", table$code), ]
    endtime<-table[which(table$code == "END"),][,1]
    
    O8starttable<-O7starttable
    # O8starttable<-O8starttable[-(which(O8starttable$code=="O7(start)" & O8starttable$time %in% OFFendTimes)),]
    # O8starttable<-O8starttable[-(which(O8starttable$code=="O7(end)" & O8starttable$time %in% OFFstartTimes)),]
    
    if (BirdWasOffatBeginning>0){
      if (isTRUE(all.equal(firstO7,TimeBirdFlewON,tolerance = 0.0009))){
        if (isTRUE(all.equal(lastO7,endtime,tolerance = 0.0009))) { #If O7 begins and ends a trial, substiutions should not be made on those obs
          O8starttable<-O8starttable[!(O8starttable$code == "O7(start)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, TimeBirdFlewON, tolerance = 0.0009)))), ]
          O8starttable<-O8starttable[!(O8starttable$code == "O7(end)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]
          O8starttable$code<-(gsub("O7(end)", "O8(start)", O8starttable$code, fixed = TRUE)) 
          O8starttable$code<-(gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE))#substituting O8starts for O7 ends, etc
          
          O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
          
          
          
        } else { #If O7 begins a trial, substitutions should not be made on that first obs, and O8(end) should be at the end
          O8starttable<-O8starttable[!(O8starttable$code == "O7(start)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, TimeBirdFlewON, tolerance = 0.0009)))), ]
          O8starttable$code<-(gsub("O7(end)", "O8(start)", O8starttable$code, fixed = TRUE)) 
          O8starttable$code<-(gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE))#substituting O8starts for O7 ends, etc
          O8starttable$code<-(gsub("END", "O8(end)", O8starttable$code, fixed = TRUE)) #substituting O8starts for O7 ends, etc
          
          O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
          
          O8ends<-O8starttable[1,]
          O8ends$time<-endtime
          O8ends$code<-"O8(end)"
          
          O8starttable<-rbind(O8starttable,O8ends)
          
        }
      } else { #IF O7 doesn't start the trial, O8 must
        if (isTRUE(all.equal(lastO7,endtime,tolerance = 0.0009))) {
          O8startsatzero<-O8starttable[1,]
          O8startsatzero$time<-TimeBirdFlewON
          O8startsatzero$code<-"O8(start)" 
          O8starttable<-rbind(O8startsatzero,O8starttable)#This (+3 lines above) adds an O8(start) at t=0
          
          O8starttable<-O8starttable[!(O8starttable$code == "O7(end)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]#Don't replace O7(end) if it's at the END
          O8starttable$code<-(gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE))
          O8starttable$code<-(gsub("O7(end)", "O8(start)", O8starttable$code, fixed = TRUE))
          
          O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
          
          
        } else { #When O7 neither starts the trial nor ends the trial, all subs can be made, END can be turned into O8(end), and O8(start) can be added at t=0
          O8startsatzero<-O8starttable[1,]
          O8startsatzero$time<-TimeBirdFlewON
          O8startsatzero$code<-"O8(start)" 
          O8starttable<-rbind(O8startsatzero,O8starttable)#This (+3 lines above) adds an O8(start) at t=0
          
          O8starttable$code<-(gsub("O7(end)", "O8(start)", gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE), fixed = TRUE)) #substituting O8starts for O7 ends, etc
          O8starttable$code<-(gsub("END", "O8(end)", O8starttable$code, fixed = TRUE)) #If 05 doesn't end the trial, O8 must
          
          O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
          
          O8ends<-O8starttable[1,]
          O8ends$time<-endtime
          O8ends$code<-"O8(end)"
          
          O8starttable<-rbind(O8starttable,O8ends)
          
          
        }
      }
    } else {
      if (isTRUE(all.equal(firstO7,0,tolerance = 0.0009))){
        if (isTRUE(all.equal(lastO7,endtime,tolerance = 0.0009))) { #If O7 begins and ends a trial, substiutions should not be made on those obs
          O8starttable<-O8starttable[!(O8starttable$code == "O7(start)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, 0, tolerance = 0.0009)))), ]
          O8starttable<-O8starttable[!(O8starttable$code == "O7(end)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]
          O8starttable$code<-(gsub("O7(end)", "O8(start)", O8starttable$code, fixed = TRUE)) 
          O8starttable$code<-(gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE))#substituting O8starts for O7 ends, etc
          
          O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
          
          
        } else { #If O7 begins a trial, substitutions should not be made on that first obs, and O8(end) should be at the end
          O8starttable<-O8starttable[!((O8starttable$code == "O7(end)" & abs(O8starttable$time - endtime) < TTolerance)), ]
          O8starttable$code<-(gsub("O7(end)", "O8(start)", O8starttable$code, fixed = TRUE)) 
          O8starttable$code<-(gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE))#substituting O8starts for O7 ends, etc
          O8starttable$code<-(gsub("END", "O8(end)", O8starttable$code, fixed = TRUE)) #substituting O8starts for O7 ends, etc
          
          O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
          
          O8ends<-O8starttable[1,]
          O8ends$time<-endtime
          O8ends$code<-"O8(end)"
          
          O8starttable<-rbind(O8starttable,O8ends)
          
        }
      } else { #IF O7 doesn't start the trial, O8 must
        if (isTRUE(all.equal(lastO7,endtime,tolerance = 0.0009))) {
          O8startsatzero<-O8starttable[1,]
          O8startsatzero$time<-0
          O8startsatzero$code<-"O8(start)" 
          O8starttable<-rbind(O8startsatzero,O8starttable)#This (+3 lines above) adds an O8(start) at t=0
          
          O8starttable<-O8starttable[!(O8starttable$code == "O7(end)" & sapply(O8starttable$time, function(x) isTRUE(all.equal(x, endtime, tolerance = 0.0009)))), ]#Don't replace O7(end) if it's at the END
          O8starttable$code<-(gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE))
          O8starttable$code<-(gsub("O7(end)", "O8(start)", O8starttable$code, fixed = TRUE))
          
          O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
          
          
        } else { #When O7 neither starts the trial nor ends the trial, all subs can be made, END can be turned into O8(end), and O8(start) can be added at t=0
          O8startsatzero<-O8starttable[1,]
          O8startsatzero$time<-0
          O8startsatzero$code<-"O8(start)" 
          O8starttable<-rbind(O8startsatzero,O8starttable)
          #This (+3 lines above) adds an O8(start) at t=0
          
          O8ends<-O8starttable[1,]
          O8ends$time<-endtime
          O8ends$code<-"O8(end)"
          
          O8starttable$code<-(gsub("O7(end)", "O8(start)", gsub("O7(start)", "O8(end)", O8starttable$code, fixed = TRUE), fixed = TRUE)) #substituting O8starts for O7 ends, etc
          O8starttable$code<-(gsub("END", "O8(end)", O8starttable$code, fixed = TRUE)) #If 05 doesn't end the trial, O8 must
          
          O8starttable$class<-(gsub("0", "5", O8starttable$class, fixed = TRUE))
          
          O8starttable<-rbind(O8starttable,O8ends)
          
          
        }
      }
    }
  } else {
    if(BirdWasOffatBeginning>0){
      O8starttable<-data.frame("time"=c(TimeBirdFlewON,TotalTime),"code"=c("O8(start)","O8(end)"),"class"=c(5,5))
    } else {
      O8starttable<-data.frame("time"=c(0,TotalTime),"code"=c("O8(start)","O8(end)"),"class"=c(5,5))
    }
  } 
  
  
  
  table5<-rbind(table,O8starttable)	
  
  table<-table5 #keeps table our working, modified dataframe
  table<-unique(table) #ELIMINATES DUPLICATE ROWS
  table<-table[order(table$time),]	
  
  
  #table<-fixOFFs(table)
  
  
  
  #############################################################	
  
  ##########REMOVES instances where there are multiple *starts* in a row
  # 	OrientationDurations<-c("O3","O4","O5","O6","O7","O8")	
  # 	
  # 	for (OD2 in OrientationDurations){  
  # 	  behav.observations <- table[grep(OD2, table$code), ]
  # 	  goo<-rle(behav.observations$code)
  # 	  
  # 	  firstone2<-behav.observations[(which(rep(goo$lengths>1,times=goo$lengths))),]
  # 	  firstone<-firstone2[1,]
  # 	  alternators<-behav.observations[-(which(rep(goo$lengths>1,times=goo$lengths))),]
  # 	  totalsequence<-rbind(firstone,alternators)
  # 	  totalsequence<-totalsequence[order(totalsequence$time),]    
  # 	  
  # 	
  # 	  kill<-firstone2[-1,]
  # 	  table<-table[!(table$time %in% kill$time & table$code %in% kill$code),]
  # 	  
  # 	}	
  ######################################################
  ####################################################
  #ELIMINATE FAULTILY CREATED DEPENDENT VARIABLES (e.g. O7(end) when 08(start) IF O8(start) = 0 )
  EndList<-c("BP1(end)","BP3(end)",
             "SS1(end)","SS2(end)",
             "O3(end)","O4(end)",
             "O5(end)","OFF(end)",
             "O6(end)","O7(end)","O8(end)",
             "OPMH.D(end)","OPMB.D(end)","OPMF.D(end)","OPMW.D(end)","OPMT.D(end)",
             "RM1(end)","RM2(end)","RM3(end)","RM4(end)",
             "MO2(end)","PU1(end)")
  
  for(xx in EndList){ #Nice little function that eliminates any "end" variables if they occur at t=0
    times2<-subset(table,code==xx)[,1]
    if(length(times2)>0){ 
      if(min(times2)==0){
        table<-table[-which(table$code == xx & table$time==0),]
      }
    }
  }
  
  for(xx in EndList){ #Nice little function that eliminates any "end" variables if they occur at the moment the bird flies on screen for the first time
    times<-subset(table,code==xx)[,1]
    if(length(times)>0){    
      if(min(times)==TimeBirdFlewON){
        table<-table[-which(table$code == xx & table$time==TimeBirdFlewON),]
      }
    }
  }	
  
  StartList2<-c("BP1(start)","BP3(start)",
                "SS1(start)","SS2(start)",
                "O3(start)","O4(start)",
                "O5(start)","OFF(end)",
                "O6(start)","O7(start)","O8(start)",
                "OPMH.D(start)","OPMB.D(start)","OPMF.D(start)","OPMW.D(start)","OPMT.D(start)",
                "RM1(start)","RM2(start)","RM3(start)","RM4(start)",
                "MO2(start)","PU1(start)")
  
  StartList<-c("BP1(start)","BP3(start)",
               "SS1(start)","SS2(start)",
               "O3(start)","O4(start)",
               "O5(start)",
               "O6(start)","O7(start)","O8(start)",
               "OPMH.D(start)","OPMB.D(start)","OPMF.D(start)","OPMW.D(start)","OPMT.D(start)",
               "RM1(start)","RM2(start)","RM3(start)","RM4(start)",
               "MO2(start)","PU1(start)")
  
  table<-fixOFFs(table)	#Fixes OFF issues
  
  
  for(xj in StartList){ #Nice little function that eliminates any "start" variables if they occur at the moment the bird flies off screen
    starttimes<-subset(table,code==xj)[,1]
    flyofftimes<-subset(table,code=="OFF(start)")[,1]
    if(length(starttimes)>0 & length(flyofftimes)){   
      table<-table[!(table$time %in% flyofftimes & table$code==xj),]        
    }
  }
  ############################################
  ###############################
  secondcheck<-table
  
  
  #############################################################  
  #OFF PERIODS---HOW MANY AND HOW MANY BEHAVIORS WERE MEASURED WHILE "OFF"
  
  happenedwhileoff<-data.frame()
  offstarts<-nrow(subset(table,code=="OFF(start)"))
  offstops<-nrow(subset(table,code=="OFF(end)"))
  off.ss.eq<-(offstarts==offstops)
  for (v in 1:nrow(subset(table,code=="OFF(start)"))){
    stopBegin<-subset(table,code=="OFF(start)")[v,1]
    stopEnd<-subset(table,code=="OFF(end)")[v,1]
    
    gap.wise.misses<-table[which(table$time < stopEnd & table$time > stopBegin),]
    happenedwhileoff<-rbind(happenedwhileoff,gap.wise.misses)
  }
  
  happenedwhileoff<-paste(happenedwhileoff$time,collapse="|")
  # 	
  #############################################################	
  table<-unique(table) #ELIMINATES DUPLICATE ROWS	
  
  

  #################
  #Additional summary measures
  #####################################
  if (BirdWasOffatBeginning>0){
       table$time<-table$time-TimeBirdFlewON+0.1
  }
  

  compositebehavioralmetrictable<-table[table$code!="END" & table$code!="Mpres" & table$code!="Fpres"&
                                          # table$code!="OFF(start)"& table$code!="OFF(end)" &
                                          table$code!="Fdisp",]
  
  compositebehavioralmetrictable$time<-round(compositebehavioralmetrictable$time,digit=1)
  
  #adds 0.1 sec to OFF(starts) to avoid overlap with behaviors that end as birds goes OFF
  compositebehavioralmetrictable[compositebehavioralmetrictable$code=="OFF(start)","time"]<-compositebehavioralmetrictable[compositebehavioralmetrictable$code=="OFF(start)","time"]+0.1
  
  #subtracts 0.1 sec to OFF(end)s to avoid overlap with behaviors that start as bird comes back on
  compositebehavioralmetrictable[compositebehavioralmetrictable$code=="OFF(end)","time"]<-compositebehavioralmetrictable[compositebehavioralmetrictable$code=="OFF(end)","time"]-0.1
  
 
  
########################### ###########################  ###########################  ###########################   
  DurationBehaviors.nosp.no078<-c("O5","O6","O3","O4",
                                  "BP1","BP3","OFF","SS1","SS2",
                                  "OPMH.D","OPMB.D","OPMF.D","OPMW.D","OPMT.D",
                                  #"RepWing","RepTail","RepHead","RepTorso",
                                  "MO2","PU1",
                                  "RM1","RM2","RM3","RM4")
  
  AllBehaviors<-c("O5","O6","O3","O4",
                  "BP1","BP2","BP3",
                  "SS1","SS2",
                  "OFF", "O2",
                  
                  #"O1","O2",
                  "OPMH","OPMB","OPMF","OPMW","OPMT",
                  "OPAH","OPAC1","OPAC2","OPAC3","OPAW","OPATW","OPAMW","OPAMTT",
                  "OPABH","OPABB","OPABF","OPABW","OPABT",
                  "PC1","PC2","PC3","PC4",
                  #"RepWing","RepTail","RepHead","RepTorso",
                  "RM1","RM2","RM3","RM4",
                  "MO1","MO2","PU1",
                  "SP1","SP2","SP3","SP4")
  
  
  
  allbehaviorstolatercombine<-length(AllBehaviors)

  EventBehaviors.Fill<-c("BP2","O1","O2",
                         "RM1","RM2","RM3","RM4",
                         "OPMH","OPMB","OPMF","OPMW","OPMT","OPAH",
                         "OPAC1","OPAC2","OPAC3","OPAW","OPATW","OPAMW","OPAMTT",
                         "OPABH","OPABB","OPABF","OPABW","OPABT",
                         "PC1","PC2","PC3","PC4",
                         "MO1","SP1","SP2","SP3","SP4")
  
  event.length<-length(EventBehaviors)
  
  ########
  compositebehaviorcategorizer<-function(dataframe,x){
    x2<-dataframe[,x,drop=FALSE]
    yyy<-as.data.frame(x2[complete.cases(x2),])
    behav<-paste(rownames(yyy),sep=".",collapse = ".")
    
  }  
  ########################################################  
  #Special class (0,1,2,3) determines which specials are treated as duration behaviors  
  if(specialclass==0){
    DurationBehaviors.toUse<-DurationBehaviors.nosp.no078
  }
  if(specialclass==1){
    DurationBehaviors.toUse<-c(DurationBehaviors.nosp.no078,"SP3")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP3"]<-"SP3(start)")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP4"]<-"SP3(end)")
  }
  if(specialclass==2){
    DurationBehaviors.toUse<-c(DurationBehaviors.nosp.no078,"SP1")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP1"]<-"SP1(start)")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP2"]<-"SP1(end)")
  }
  if(specialclass==3){
    DurationBehaviors.toUse<-c(DurationBehaviors.nosp.no078,c("SP1","SP3"))
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP1"]<-"SP1(start)")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP2"]<-"SP1(end)")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP3"]<-"SP3(start)")
    compositebehavioralmetrictable<-within(compositebehavioralmetrictable,code[code=="SP4"]<-"SP3(end)")
    
    
  }
  ###########################################################   
  ########################### ###########################  ###########################  ###########################   
  orderedtable<-compositebehavioralmetrictable[order(compositebehavioralmetrictable$time),]  
  
  tabletosplit<-orderedtable  
  row.n<-which(tabletosplit$code=="OFF(start)")  
  split.list<-list()
  
#IF THERE ARE ANY OFF PERIODS LEFT IN THE BEHAVIORAL TABLE (there should be no start or end OFFs)
  #THEN THIS LOOP SPLITS THE SECTIONS OUT INTO *split.list* as behavioral sequences
  if(sum(row.n)>0) {
    for(offsplit in 1:length(row.n)){
      
      thisrow<-which(tabletosplit$code=="OFF(start)")[1]  
      
      if(!is.na(thisrow)){
      splittable<-tabletosplit[c(1:thisrow),]#Splits out all behaviors before OFF(start)
      splittable<-splittable[c(!grepl("OFF",splittable[,2])),]#gets rid of OFF behaviors in split
      } else {
        splittable<-tabletosplit #Splits out last set of behaviors 
        splittable<-splittable[c(!grepl("OFF",splittable[,2])),]#gets rid of OFF behaviors in split
      }
      
      if(offsplit>1){
        splittable$time<-splittable$time-splittable$time[1]
      }
      
      TotalTime<-splittable[length(splittable[,1]),"time"]- splittable[1,"time"]#calculates total time of split

      split.info.start<-splittable[1,"time"]
      split.info.end<-splittable[length(splittable[,1]),"time"]
      split.info.period<-paste(split.info.start,"-",split.info.end,sep='')
      
      newlength<-round(TotalTime/0.1)+2 #Complete # of 0.1 sec intervals for the video (adds 1, b/c we add 0.1 to the first behavior so it didn't occur at 0 = 0th column)
      completebehaviormatrix<-matrix(data=NA,nrow=allbehaviorstolatercombine,ncol=newlength)
      row.names(completebehaviormatrix)<-AllBehaviors
      
      #Adds duration behaviors to completebehaviormatrix
      for (cb in 1:length(DurationBehaviors.toUse)){ #loop through all initially scored behaviors, but not Specials
        
        thisone<-DurationBehaviors.toUse[cb] #gives new name, "thisone" to the current behavior
        
        start.code<-paste(thisone,"(start)",sep="")
        end.code<-paste(thisone,"(end)",sep="")
        
        starts<-splittable[splittable$code==start.code, ] #df of starts
        ends<-splittable[splittable$code==end.code, ] #df of ends
        
        number.of.starts<-length(starts[,2])
        
        if(cb>16){ #chooses only RM1, RM2, RM3, RM4
          if(number.of.starts!=0){  #Only run the "fill-in loop" if the behavior was recorded at all 
            if((number.of.starts>0 & !length(ends[,2])==0)){  #oNLY IF THERE ARE START AND ENDS FOR RM1-4 DO WE FILL
              for (ooo in 1:number.of.starts)  { #This loop fills the matrix with 1s for every 1/10 sec where the behavior is occuring
                
                columnstofill<-10* (seq(starts[ooo,1],ends[ooo,1],by=0.1)   )
                
                if(grepl(".D",thisone)){
                  thisone<-substr(thisone,1,nchar(thisone)-2)
                } 
                completebehaviormatrix[rownames(completebehaviormatrix)==thisone,columnstofill]<-1 #fills columns while behavior is active
              }
            }
          }
        } else {
          if(number.of.starts!=0){  #Only run the "fill-in loop" if the behavior was recorded at all    
            for (ooo in 1:number.of.starts)  { #This loop fills the matrix with 1s for every 1/10 sec where the behavior is occuring
              
              if(!is.na(ends[ooo,1]) & ends[ooo,1]>starts[ooo,1]){ #Proceed if end time EXISTS AND is AFTER beginning time
                
                columnstofill<-10* (seq(starts[ooo,1],ends[ooo,1],by=0.1)   )
                
                if(grepl(".D",thisone)){
                  thisone<-substr(thisone,1,nchar(thisone)-2)
                } 
                
                completebehaviormatrix[rownames(completebehaviormatrix)==thisone,columnstofill]<-1 #fills columns while behavior is active
              } else {
                columnstofill<-10* (seq(starts[ooo,1],splittable[length(splittable[,1]),1],by=0.1)   )
                if(grepl(".D",thisone)){
                  thisone<-substr(thisone,1,nchar(thisone)-2)
                } 
                
                completebehaviormatrix[rownames(completebehaviormatrix)==thisone,columnstofill]<-1 #fills columns while behavior is active
              }
              
            }
          }
        }

      }  
      
      #Adds event behaviors to completebehaviormatrix
      for (ebs in 1:event.length){
        
        behaviortopull<-EventBehaviors.Fill[ebs]
        if(ebs>3 & ebs<7){ #for RM behaviors only
          
          start.code<-paste(behaviortopull,"(start)",sep="")
          end.code<-paste(behaviortopull,"(end)",sep="")
          
          starts<-splittable[splittable$code==start.code, ] #df of starts
          ends<-splittable[splittable$code==end.code, ] #df of ends
          
          number.of.starts<-length(starts[,2])
          if((number.of.starts>0 & length(ends[,2])==0)){ #If there are starts, and no ends, then each RM is an event
            behavs.to.fill<-splittable[splittable$code==start.code, ]
            column.to.add<-10*behavs.to.fill$time
            column.to.add<-ifelse(column.to.add==0,1,column.to.add)
            
            completebehaviormatrix[rownames(completebehaviormatrix)==behaviortopull,column.to.add]<-1 #fills columns while behavior is active
          }
        } else {
          
          
          
          behavs.to.fill<-splittable[splittable$code==behaviortopull, ]
          column.to.add<-10*behavs.to.fill$time
          column.to.add<-ifelse(column.to.add==0,1,column.to.add)
          
          completebehaviormatrix[rownames(completebehaviormatrix)==behaviortopull,column.to.add]<-1 #fills columns while behavior is active
        }
        
      }

      temporalbehaviorcategories<-matrix(data=NA,nrow=1,ncol=newlength)
      for(pdq in 1:newlength){
        temporalbehaviorcategories[1,pdq]<-compositebehaviorcategorizer(completebehaviormatrix,pdq)
      }
 
      #Creates *no.OFFs* matrix----and to be sure, gets rid of any "OFFs"
      no.OFFs<-as.matrix(temporalbehaviorcategories[,!grepl("OFF",temporalbehaviorcategories[1,])])#gets rid of behavior combos with OFF
      no.OFFs[no.OFFs == ""] <- NA
      no.OFFs<-na.omit(no.OFFs)
      ###########################################################################################################
      #Finds behaviors missing a body-axis (O5/O6) and adds axis from nearest (time-wise) axis
      deficient.behaviors<-no.OFFs[((!grepl("O5", no.OFFs[,1], fixed=TRUE)) & (!grepl("O6", no.OFFs[,1], fixed=TRUE)) ),]
      containing.behaviors<-no.OFFs[((grepl("O5", no.OFFs[,1], fixed=TRUE)) | (grepl("O6", no.OFFs[,1], fixed=TRUE)) ),]
      
      index.of.deficient<-which((!grepl("O5", no.OFFs[,1], fixed=TRUE)) & (!grepl("O6", no.OFFs[,1], fixed=TRUE)) )
      index.of.containing<-which((grepl("O5", no.OFFs[,1], fixed=TRUE)) | (grepl("O6", no.OFFs[,1], fixed=TRUE)) )
      
      ####
      #Little loop through all deficient (posture-absent) behaviors, adding posture from nearest time
      if(length(index.of.deficient)>0){
        for(defic in 1:length(index.of.deficient)){
          defbehav<-deficient.behaviors[defic]
          defind<-index.of.deficient[defic]
          
          replacementbehavior<-no.OFFs[which.min(abs(index.of.containing-defind)),]
          replacementposture<-ifelse(grepl("O5",replacementbehavior),"O5","O6")
          #replacementpostureDOT<-paste(replacementposture,".", sep ="")
          
          
          
          newbehaviorvalue<- paste(replacementposture,defbehav,sep=".")
          
          
          
          
          no.OFFs[defind,1]<-newbehaviorvalue
          
        }
      }  
      ####################################################
      #######################################################
      #no.OFFs<-as.matrix(gsub("O5.O6","O5",no.OFFs[,1]))#swaps out impossible O5/O6
      
      #need to deal with O5.O6s (maybe randomly choose one?)
      no.OFFs[((grepl("O5.O6", no.OFFs[,1], fixed=TRUE)) ),]
      simultaneousO5O6index<-grep("O5.O6", no.OFFs[,1], fixed=TRUE,value=FALSE)
      fixerO5O6er<-c("O5","O6")
      
      if(length(simultaneousO5O6index)>0){
        for(abc in 1:length(simultaneousO5O6index)){
          no.OFFs[simultaneousO5O6index[abc],1]<-as.matrix(gsub("O5.O6",fixerO5O6er[sample(1:2,1)],no.OFFs[simultaneousO5O6index[abc],1]))
        }
      } #randomly replaces any simultatneous O5.O6 with either O5 or O6
      
      FileisNamedSplit<-paste("split",offsplit,":",split.info.period,":",FileisNamed,sep="")
      
#       lab<-rbind(FileisNamedSplit,no.OFFs)
#       rownames(lab) <- c()
      
      split.list[[offsplit]]<-no.OFFs
      names(split.list)[offsplit]<-FileisNamedSplit
      
      tabletosplit<-tabletosplit[-c(1:thisrow),] #removes split that we just put into the new list
    }
  }  else {
    splittable<-orderedtable
    TotalTime<-splittable[length(splittable[,1]),"time"]- splittable[1,"time"]#calculates total time of split
    
    newlength<-round(TotalTime/0.1)+2 #Complete # of 0.1 sec intervals for the video (adds 1, b/c we add 0.1 to the first behavior so it didn't occur at 0 = 0th column)
    completebehaviormatrix<-matrix(data=NA,nrow=allbehaviorstolatercombine,ncol=newlength)
    row.names(completebehaviormatrix)<-AllBehaviors
    
    #Adds duration behaviors to completebehaviormatrix
    for (cb in 1:length(DurationBehaviors.toUse)){ #loop through all initially scored behaviors, but not Specials
      
      thisone<-DurationBehaviors.toUse[cb] #gives new name, "thisone" to the current behavior
      
      start.code<-paste(thisone,"(start)",sep="")
      end.code<-paste(thisone,"(end)",sep="")
      
      starts<-splittable[splittable$code==start.code, ] #df of starts
      ends<-splittable[splittable$code==end.code, ] #df of ends
      
      number.of.starts<-length(starts[,2])
      

      if(cb>16){ #chooses only RM1, RM2, RM3, RM4
        if(number.of.starts!=0){  #Only run the "fill-in loop" if the behavior was recorded at all 
          if((number.of.starts>0 & !length(ends[,2])==0)){  #oNLY IF THERE ARE START AND ENDS FOR RM1-4 DO WE FILL
            for (ooo in 1:number.of.starts)  { #This loop fills the matrix with 1s for every 1/10 sec where the behavior is occuring
              
              columnstofill<-10* (seq(starts[ooo,1],ends[ooo,1],by=0.1)   )
              
              if(grepl(".D",thisone)){
                thisone<-substr(thisone,1,nchar(thisone)-2)
              } 
              completebehaviormatrix[rownames(completebehaviormatrix)==thisone,columnstofill]<-1 #fills columns while behavior is active
            }
          }
        }
      } else {
              if(number.of.starts!=0){  #Only run the "fill-in loop" if the behavior was recorded at all    
                for (ooo in 1:number.of.starts)  { #This loop fills the matrix with 1s for every 1/10 sec where the behavior is occuring
                  
                  columnstofill<-10* (seq(starts[ooo,1],ends[ooo,1],by=0.1)   )
                  
                  if(grepl(".D",thisone)){
                    thisone<-substr(thisone,1,nchar(thisone)-2)
                  } 
                  completebehaviormatrix[rownames(completebehaviormatrix)==thisone,columnstofill]<-1 #fills columns while behavior is active
                }
              }
            }
    }  
    
    #Adds event behaviors to completebehaviormatrix
    for (ebs in 1:event.length){
      
      behaviortopull<-EventBehaviors.Fill[ebs]
      
      if(ebs>3 & ebs<7){ #for RM behaviors only
       
        start.code<-paste(behaviortopull,"(start)",sep="")
        end.code<-paste(behaviortopull,"(end)",sep="")
        
        starts<-splittable[splittable$code==start.code, ] #df of starts
        ends<-splittable[splittable$code==end.code, ] #df of ends
        
        number.of.starts<-length(starts[,2])
        if((number.of.starts>0 & length(ends[,2])==0)){ #If there are starts, and no ends, then each RM is an event
          behavs.to.fill<-splittable[splittable$code==start.code, ]
          column.to.add<-10*behavs.to.fill$time
          column.to.add<-ifelse(column.to.add==0,1,column.to.add)
          
          completebehaviormatrix[rownames(completebehaviormatrix)==behaviortopull,column.to.add]<-1 #fills columns while behavior is active
        }
      } else {
      
      
      
      behavs.to.fill<-splittable[splittable$code==behaviortopull, ]
      column.to.add<-10*behavs.to.fill$time
      column.to.add<-ifelse(column.to.add==0,1,column.to.add)
      
      completebehaviormatrix[rownames(completebehaviormatrix)==behaviortopull,column.to.add]<-1 #fills columns while behavior is active
      }
    }
    
    temporalbehaviorcategories<-matrix(data=NA,nrow=1,ncol=newlength)
    for(pdq in 1:newlength){
      temporalbehaviorcategories[1,pdq]<-compositebehaviorcategorizer(completebehaviormatrix,pdq)
    }
    
    #Creates *no.OFFs* matrix----and to be sure, gets rid of any "OFFs"
    no.OFFs<-as.matrix(temporalbehaviorcategories[,!grepl("OFF",temporalbehaviorcategories[1,])])#gets rid of behavior combos with OFF
    no.OFFs[no.OFFs == ""] <- NA
    no.OFFs<-na.omit(no.OFFs)
    ###########################################################################################################
    #Finds behaviors missing a body-axis (O5/O6) and adds axis from nearest (time-wise) axis
    deficient.behaviors<-no.OFFs[((!grepl("O5", no.OFFs[,1], fixed=TRUE)) & (!grepl("O6", no.OFFs[,1], fixed=TRUE)) ),]
    containing.behaviors<-no.OFFs[((grepl("O5", no.OFFs[,1], fixed=TRUE)) | (grepl("O6", no.OFFs[,1], fixed=TRUE)) ),]
    
    index.of.deficient<-which((!grepl("O5", no.OFFs[,1], fixed=TRUE)) & (!grepl("O6", no.OFFs[,1], fixed=TRUE)) )
    index.of.containing<-which((grepl("O5", no.OFFs[,1], fixed=TRUE)) | (grepl("O6", no.OFFs[,1], fixed=TRUE)) )
    
    ####
    #Little loop through all deficient (posture-absent) behaviors, adding posture from nearest time
    if(length(index.of.deficient)>0){
      for(defic in 1:length(index.of.deficient)){
        defbehav<-deficient.behaviors[defic]
        defind<-index.of.deficient[defic]
        
        replacementbehavior<-no.OFFs[which.min(abs(index.of.containing-defind)),]
        replacementposture<-ifelse(grepl("O5",replacementbehavior),"O5","O6")
        #replacementpostureDOT<-paste(replacementposture,".", sep ="")
        
        
        
        newbehaviorvalue<- paste(replacementposture,defbehav,sep=".")
        
        
        
        
        no.OFFs[defind,1]<-newbehaviorvalue
        
      }
    }  
    ####################################################
    #######################################################
    #no.OFFs<-as.matrix(gsub("O5.O6","O5",no.OFFs[,1]))#swaps out impossible O5/O6
    
    #need to deal with O5.O6s (maybe randomly choose one?)
    no.OFFs[((grepl("O5.O6", no.OFFs[,1], fixed=TRUE)) ),]
    simultaneousO5O6index<-grep("O5.O6", no.OFFs[,1], fixed=TRUE,value=FALSE)
    fixerO5O6er<-c("O5","O6")
    
    if(length(simultaneousO5O6index)>0){
      for(abc in 1:length(simultaneousO5O6index)){
        no.OFFs[simultaneousO5O6index[abc],1]<-as.matrix(gsub("O5.O6",fixerO5O6er[sample(1:2,1)],no.OFFs[simultaneousO5O6index[abc],1]))
      }
    } #randomly replaces any simultatneous O5.O6 with either O5 or O6
    
    split.info.period<-paste(splittable[1,"time"],splittable[length(splittable[,1]),"time"],sep='-')
    FileisNamedSplit<-paste("split",1,":",split.info.period,":",FileisNamed,sep="")
    #FileisNamedSplit<-paste(1,FileisNamed,sep="_")
#     lab<-rbind(FileisNamedSplit,no.OFFs)
#     rownames(lab) <- c()
    split.list[[1]]<-no.OFFs
    names(split.list)[1]<-FileisNamedSplit
  }
  
  

 
  return(split.list)
}


StreamLinedProcessBOP<-function(no.OFFs){


  #Counts occurences for each ID type
  summarytable<-table(no.OFFs)
  behaviorcounts<-table(no.OFFs)
  proportiontable<-prop.table(summarytable)
  uniquebehaviors<-length(summarytable)
  
  six<-no.OFFs
  #ALL BEHAVIORS and TRANSITIONS
  if(uniquebehaviors>1){
    
    ##This loop creates the denominator for the calculation of the inverse Simpson Index for behavioral diversity
    val2<-0
    for (u in 1:uniquebehaviors){
      proptosquare<-proportiontable[u]
      val <-proptosquare^2.0
      val2<-val+val2
    }			
    #inverse Simpson Index for behavioral diversity
    iSimpson<-1.0/as.numeric(val2)
    relativebehavioraldiversity<-iSimpson/(uniquebehaviors) #relative color diversity, accounting for # of color classes
    
    ##########################################################################################################
    #uses "createSequenceMatrix" function (from *markovchain*) to calculate transition matrix
    c.no.OFFs<-as.character(no.OFFs)
    
    
    TransitionMatrix.markov<-createSequenceMatrix(c.no.OFFs,toRowProbs = TRUE,sanitize=FALSE)
    TransitionMatrix.cum<-createSequenceMatrix(c.no.OFFs,sanitize=FALSE)
    
    #TransitionMatrix<-TransitionMatrix[-ncol(TransitionMatrix),-ncol(TransitionMatrix)]
    
    
    time.value<-sum(TransitionMatrix.cum)/10 #seconds watched
    
    
    noselfs<-TransitionMatrix.cum
    diag(noselfs)<-NA
    
    prop.table.excludeNAs<- function(x) {x/sum(x, na.rm=TRUE)}
    proportiontable.trans.Simpson<-prop.table.excludeNAs(noselfs) #get rid of diagonal/self-transitions for Simpson Analyses
    
    proportiontable.trans<-prop.table.excludeNAs(TransitionMatrix.cum) 
    
    
    if(uniquebehaviors>1){
      val3<-matrix(nrow=(((uniquebehaviors^2)-uniquebehaviors)/2),ncol=2)
      rownum<-1
      for (u in 1:uniquebehaviors){
        for(v in u:uniquebehaviors){#iteratively reduces columns analyzed in next loop to avoid double counting transitions 
          if (u!=v){
            val3[rownum,2]<-as.numeric(proportiontable.trans.Simpson[u,v]+proportiontable.trans.Simpson[v,u])
            val3[rownum,1]<-paste(row.names(proportiontable.trans.Simpson)[u],colnames(proportiontable.trans.Simpson)[v],sep="-")
            rownum<-rownum+1
          }
        }
      }
      val4<-as.data.frame(val3)
      val4[,2]<-as.numeric(as.character(val4[,2])) 
      
      simpsonT<-0
      
      for(z in 1:nrow(val4)){
        propsq<-val4[z,2]^2.0
        simpsonT<-simpsonT+propsq
      }
      
      
      iSimpsonT<-1.0/simpsonT #Calculates inverse Simpson index for horizontal transitions
      relativetransitiondiversity<-iSimpsonT/(uniquebehaviors*(uniquebehaviors-1)/2)  
    } else {
      
      iSimpsonT <- NA
      relativetransitiondiversity<- NA
    }
    
    
    PropSelfTransitioning<-sum(diag(TransitionMatrix.cum))/sum(TransitionMatrix.cum)
    
    network.TransitionMatrix.cum<-network(proportiontable.trans,directed=TRUE,loops=TRUE,hyper=FALSE,matrix.type="adjacency")
    
    
    
    showme<-graph_from_adjacency_matrix(proportiontable.trans,mode="directed",weighted=TRUE,diag=FALSE,add.colnames=TRUE)
    showme.undirected<-graph_from_adjacency_matrix(proportiontable.trans,mode="undirected",weighted=TRUE,diag=FALSE,add.colnames=TRUE)
    
    showme.looped<- graph_from_adjacency_matrix(proportiontable.trans,mode="directed",weighted=TRUE,diag=TRUE,add.colnames=TRUE)
    showme.looped.undirected<- graph_from_adjacency_matrix(proportiontable.trans,mode="undirected",weighted=TRUE,diag=TRUE,add.colnames=TRUE)
    
    
    
    E(showme)$width <- 1+E(showme)$weight/50
    V(showme)$nombres<-colnames(proportiontable.trans)
    
    E(showme.looped)$width <- 1+E(showme.looped)$weight/50
    V(showme.looped)$nombres<-colnames(proportiontable.trans)
    
    # Compute approx self-transition and use that to set node size:
    selfsizesscaled<-(diag(proportiontable.trans))*100+1
    V(showme)$selfsize<-sqrt(selfsizesscaled)
    V(showme)$size <- ((selfsizesscaled)^2)*9
    
    V(showme.looped)$selfsize<-log(selfsizesscaled)+1
    
    
    
    #####################################################
    # Network Clustering
    # l <- layout.fruchterman.reingold(showme.looped)
    # l <- layout_with_kk(showme.looped)
    # l <- layout_with_lgl(showme.looped)
    l <- layout_with_graphopt(showme.looped)
    l2 <-layout_with_graphopt(showme)
    #l <- layout_nicely(showme.looped)
    #behaviorplotlist
    
    
    #clp <- cluster_label_prop(showme.looped)
    clp <-cluster_label_prop(showme.looped.undirected) # input graph should be undirected to make sense.
    clp2<-cluster_label_prop(showme.undirected)
    # setwd("C:/Users/Rusty/Amazon Drive/BOP/Video_analysis")
    # source("PictoGrammerMega_Dec13ab.R", chdir = F)
    # #PictoGrammer<-function(proportiontable,hangingbird)
    #
    
    
    
    
    ################################################################################################################################

    
    
    
    
    
    
    
    
    # NETWORK SUMMARY VARIABLES
    
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    #6.1 Density
    
    #The proportion of present edges from all possible edges in the network.
    
    #densityscore<-ecount(showme.looped)/(vcount(showme.looped)*(vcount(showme.looped)-1)) #for a directed network
    densityscore<-ecount(showme.looped)/(vcount(showme.looped)^2)#for a directed network WITH self-loops
    
    #Mean degree
    deg <- degree(showme.looped,loops=FALSE)
    mean.degree<-mean(deg)
    
    #Mean path length
    mean_path_length<-mean_distance(showme.looped, directed=T)
    
    
    #network diameter is the longest geodesic distance (length of the shortest path between two nodes) in the network.
    network.diameter<-diameter(showme.looped, directed=T,weights=NA)
    
    
    #average clustering coefficient
    avg.cluster.coeff<-transitivity(showme.looped)
    
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
    
    printme<-data.frame(time.value,
                        uniquebehaviors,densityscore,mean.degree,mean_path_length,network.diameter,avg.cluster.coeff,
                        iSimpson,relativebehavioraldiversity,iSimpsonT,relativetransitiondiversity,
                        PropSelfTransitioning,Smallworldness)
    
    #CONNECTIVITY measures reflect degree to which network differs from complete network
    #Edge density: % of edges compared to maximum
    #Average degree: Avg. number of links
    #Average path length: avg of shortest pathes between reachable nodes
    #Network diameter: longest of shortest paths
    
    #CENTRALITY measures quantify heterogeneity in network structure
    #Average clustering coefficient:
    #Components
    
    
    
    
    
    
    
    
    
    
    
    #V(showme.looped)$size<-80
#     netname<-paste(name2,"Network.pdf",sep='_')
#     
#     setwd(location)
#     pdf(netname,width= 12, height= 12,family="NimbusRom")
#     op <- par(mfrow=c(1,1))# rows, columns
    # par(mar=c(0,0,0,0)) 
    # plot(clp, showme.looped,edge.curved=.1,edge.arrow.size = 0.45,vertex.size=18,layout=l,rescale=FALSE,
    #      xlim=range(l[,1]),ylim=range(l[,2]))
    
    
    # plot(clp, showme.looped,edge.curved=.1,edge.arrow.size = 0.45,vertex.size=18,layout=l,rescale=TRUE,
    #      xlim=c(-1,1),ylim=c(-1,1))
    
#     V(showme.looped)$raster <- behaviorplotlist
#     V(showme)$raster <- behaviorplotlist
    
    
    
    
    # print(  plot(clp, showme.looped,edge.curved=.1,layout=l,rescale=TRUE,vertex.shape="raster",
    #        vertex.label=NA, 
    #        #vertex.size=15*(V(showme.looped)$selfsize), # resizes vertices based on self-self transitions
    #        vertex.size=25,
    #        edge.color=col.tr,
    #        edge.arrow.size=0.35,
    #        
    #        edge.width=3,
    #        edge.arrow.width=1.6,
    #        edge.loop.angle=1.5,
    #        
    #        xlim=c(-1,1),ylim=c(-1,1)))
    
#     print(  plot(clp, showme,edge.curved=.1,layout=l,rescale=TRUE,vertex.shape="raster",
#                  vertex.label=NA, 
#                  #vertex.size=15*(V(showme.looped)$selfsize), # resizes vertices based on self-self transitions
#                  vertex.size=25,
#                  edge.color=col.tr,
#                  edge.arrow.size=0.35,
#                  
#                  edge.width=3,
#                  edge.arrow.width=1.6,
#                  edge.loop.angle=1.5,
#                  
#                  xlim=c(-1,1),ylim=c(-1,1)))
#     
#     print(text(0.9,-0.75,paste("Timescored(s) = ",time.value,sep=""),col="black",cex=0.75,adj=(0)))
#     print(text(0.9,-0.8,paste("behaviors = ",uniquebehaviors,sep=""),col="black",cex=0.75,adj=(0)))  
#     print(text(0.9,-0.85,paste("mean deg = ",round(mean.degree,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
#     print(text(0.9,-0.9,paste("density score = ",round(densityscore,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
#     print(text(0.9,-0.95,paste("mean path L = ",round(mean_path_length,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
#     print(text(0.9,-1,paste("network.diameter = ",network.diameter,sep=""),col="black",cex=0.75,adj=(0)))  
#     print(text(0.9,-1.05,paste("Smallworldness = ",Smallworldness,sep=""),col="black",cex=0.75,adj=(0)))
#     
#     
    #tkplot(clp, showme.looped)
    #####  par("usr")
    
    # library(plotrix)
    # 
    # scalelist<-list()
    # gettinsmaller<-seq(from=1,to=0.52, by=-0.02)
    # for (vc in 1:25){
    # 
    #   
    #   gettingsmaller<-gettinsmaller[vc]
    #   pointlocations<-cbind(rescale(l[,1],gettingsmaller*c(-1,1)),rescale(l[,2],gettingsmaller*c(-1,1)))
    #   par(new=TRUE)
    #   plot(pointlocations,xlab="", ylab="",xlim=c(-1,1),ylim=c(-1,1),pch=vc)
    # 
    # }
    # dev.off()
    
    
    
    # for(bbb in 1:uniquebehaviors){
    #   rasterImage(behaviorplotlist[[bbb]], xleft = pointlocations[bbb,1], xright=pointlocations[bbb,1]+0.15,
    #               ybottom=pointlocations[bbb,2], ytop=pointlocations[bbb,2]+0.15)  
    #   par(new=TRUE)
    #   plot(pointlocations,xlab="", ylab="",xlim=c(-1,1),ylim=c(-1,1))
    # }
    # 
    # pointlocations<-cbind(rescale(l[,1],0.92*c(-1,1)),rescale(l[,2],0.92*c(-1,1)))
    # pointlocations<-pointlocations
    # 
    # for(bbb in 1:uniquebehaviors){
    #   rasterImage(behaviorplotlist[[bbb]], xleft = pointlocations[bbb,1], xright=pointlocations[bbb,1]+0.15,
    #               ybottom=pointlocations[bbb,2], ytop=pointlocations[bbb,2]+0.15)  
    #   par(new=TRUE)
    #   
    # }
    
    # dev.off()
    
    #Another network program
    #geph
    
  }
  
  #FILTERED BEHAVIORS (minus behaviors that only happened once)
  no.OFFs.filtered1<-no.OFFs[which(duplicated(no.OFFs)),,drop=FALSE]
  
  summarytable.filtered1<-table(no.OFFs.filtered1)
  behaviorcounts.filtered1<-table(no.OFFs.filtered1)
  proportiontable.filtered1<-prop.table(summarytable.filtered1)
  uniquebehaviors.filtered1<-length(summarytable.filtered1)  
  
  seven<-no.OFFs.filtered1
  
  if(uniquebehaviors.filtered1>1){  
    ##This loop creates the denominator for the calculation of the inverse Simpson Index for behavioral diversity
    val2.filtered1<-0
    for (u in 1:uniquebehaviors.filtered1){
      proptosquare<-proportiontable[u]
      val <-proptosquare^2.0
      val2.filtered1<-val+val2.filtered1
    }			
    #inverse Simpson Index for behavioral diversity
    iSimpson.filtered1<-1.0/as.numeric(val2.filtered1)
    relativebehavioraldiversity.filtered1<-iSimpson.filtered1/(uniquebehaviors.filtered1) #relative color diversity, accounting for # of color classes
    
    ##########################################################################################################
    #uses "createSequenceMatrix" function (from *markovchain*) to calculate transition matrix
    c.no.OFFs<-as.character(no.OFFs.filtered1)
    
    
    TransitionMatrix.markov<-createSequenceMatrix(c.no.OFFs,toRowProbs = TRUE,sanitize=FALSE)
    TransitionMatrix.cum.filtered1<-createSequenceMatrix(c.no.OFFs,sanitize=FALSE)
    
    #TransitionMatrix<-TransitionMatrix[-ncol(TransitionMatrix),-ncol(TransitionMatrix)]
    
    time.value<-sum(TransitionMatrix.cum.filtered1)/10 #seconds watched
    
    
    
    noselfs<-TransitionMatrix.cum.filtered1
    diag(noselfs)<-NA
    
    prop.table.excludeNAs<- function(x) {x/sum(x, na.rm=TRUE)}
    proportiontable.trans.Simpson<-prop.table.excludeNAs(noselfs) #get rid of diagonal/self-transitions for Simpson Analyses
    
    proportiontable.trans<-prop.table.excludeNAs(TransitionMatrix.cum.filtered1) 
    
    
    if(uniquebehaviors.filtered1>1){
      val3<-matrix(nrow=(((uniquebehaviors.filtered1^2)-uniquebehaviors.filtered1)/2),ncol=2)
      rownum<-1
      for (u in 1:uniquebehaviors.filtered1){
        for(v in u:uniquebehaviors.filtered1){#iteratively reduces columns analyzed in next loop to avoid double counting transitions 
          if (u!=v){
            val3[rownum,2]<-as.numeric(proportiontable.trans.Simpson[u,v]+proportiontable.trans.Simpson[v,u])
            val3[rownum,1]<-paste(row.names(proportiontable.trans.Simpson)[u],colnames(proportiontable.trans.Simpson)[v],sep="-")
            rownum<-rownum+1
          }
        }
      }
      val4<-as.data.frame(val3)
      val4[,2]<-as.numeric(as.character(val4[,2])) 
      
      simpsonT<-0
      
      for(z in 1:nrow(val4)){
        propsq<-val4[z,2]^2.0
        simpsonT<-simpsonT+propsq
      }
      
      
      iSimpsonT.filtered1<-1.0/simpsonT #Calculates inverse Simpson index for horizontal transitions
      relativetransitiondiversity.filtered1<-iSimpsonT.filtered1/(uniquebehaviors.filtered1*(uniquebehaviors.filtered1-1)/2)  
    } else {
      
      iSimpsonT.filtered1 <- NA
      relativetransitiondiversity.filtered1<- NA
    }
    
    
    PropSelfTransitioning.filtered1<-sum(diag(TransitionMatrix.cum.filtered1))/sum(TransitionMatrix.cum.filtered1)
    
    network.TransitionMatrix.cum.filtered1<-network(proportiontable.trans,directed=TRUE,loops=TRUE,hyper=FALSE,matrix.type="adjacency")
    
    
    
    showme<-graph_from_adjacency_matrix(proportiontable.trans,mode="directed",weighted=TRUE,diag=FALSE,add.colnames=TRUE)
    showme.undirected<-graph_from_adjacency_matrix(proportiontable.trans,mode="undirected",weighted=TRUE,diag=FALSE,add.colnames=TRUE)
    
    showme.looped<- graph_from_adjacency_matrix(proportiontable.trans,mode="directed",weighted=TRUE,diag=TRUE,add.colnames=TRUE)
    showme.looped.undirected<- graph_from_adjacency_matrix(proportiontable.trans,mode="undirected",weighted=TRUE,diag=TRUE,add.colnames=TRUE)
    
    
    
    E(showme)$width <- 1+E(showme)$weight/50
    V(showme)$nombres<-colnames(proportiontable.trans)
    
    E(showme.looped)$width <- 1+E(showme.looped)$weight/50
    V(showme.looped)$nombres<-colnames(proportiontable.trans)
    
    # Compute approx self-transition and use that to set node size:
    selfsizesscaled<-(diag(proportiontable.trans))*100+1
    V(showme)$selfsize<-sqrt(selfsizesscaled)
    V(showme)$size <- ((selfsizesscaled)^2)*9
    
    V(showme.looped)$selfsize<-log(selfsizesscaled)+1
    
    
    
    #####################################################
    # Network Clustering
    # l <- layout.fruchterman.reingold(showme.looped)
    # l <- layout_with_kk(showme.looped)
    # l <- layout_with_lgl(showme.looped)
    l <- layout_with_graphopt(showme.looped)
    l2 <-layout_with_graphopt(showme)
    #l <- layout_nicely(showme.looped)
    #behaviorplotlist
    
    
    #clp <- cluster_label_prop(showme.looped)
    clp <-cluster_label_prop(showme.looped.undirected) # input graph should be undirected to make sense.
    clp2<-cluster_label_prop(showme.undirected)
    # setwd("C:/Users/Rusty/Amazon Drive/BOP/Video_analysis")
    # source("PictoGrammerMega_Dec13ab.R", chdir = F)
    # #PictoGrammer<-function(proportiontable,hangingbird)
    #
    
    
    
    

    # NETWORK SUMMARY VARIABLES
    
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    #6.1 Density
    
    #The proportion of present edges from all possible edges in the network.
    
    densityscore.filtered1<-ecount(showme.looped)/(vcount(showme.looped)*(vcount(showme.looped)-1)) #for a directed network
    
    #Mean degree
    deg <- degree(showme.looped,loops=FALSE)
    mean.degree.filtered1<-mean(deg)
    
    #Mean path length
    mean_path_length.filtered1<-mean_distance(showme.looped, directed=T)
    
    
    #network diameter is the longest geodesic distance (length of the shortest path between two nodes) in the network.
    network.diameter.filtered1<-diameter(showme.looped, directed=T,weights=NA)
    
    
    #average clustering coefficient
    avg.cluster.coeff.filtered1<-transitivity(showme.looped)
    
    ############################################
    #SMALL WORLDNESS --- video-wise, filtered
    #number of nodes/vertices in graph
    vertices<- uniquebehaviors.filtered1
    #number of edges in G(n,m) graph
    edges<- sum(TransitionMatrix.cum.filtered1!=0)
    
    rando.network.filtered1<-sample_gnm(n=vertices, m=edges, directed = TRUE, loops = TRUE)
    
    Trobserved<-avg.cluster.coeff.filtered1
    mean.Trrandom<-transitivity(rando.network.filtered1)
    
    SPobserved<-mean_path_length.filtered1
    mean.SPrandom<-mean_distance(rando.network.filtered1, directed=T)  
    
    Smallworldness.filtered1<- (Trobserved/mean.Trrandom)/(SPobserved/mean.SPrandom)
    ############################################
    
    
    printme<-data.frame(time.value,
                        uniquebehaviors,densityscore,mean.degree,mean_path_length,network.diameter,avg.cluster.coeff,
                        iSimpson,relativebehavioraldiversity,iSimpsonT,relativetransitiondiversity,
                        PropSelfTransitioning,Smallworldness,
                        
                        uniquebehaviors.filtered1,densityscore.filtered1,mean.degree.filtered1,mean_path_length.filtered1,network.diameter.filtered1,avg.cluster.coeff.filtered1,
                        iSimpson.filtered1,relativebehavioraldiversity.filtered1,iSimpsonT.filtered1,relativetransitiondiversity.filtered1,
                        PropSelfTransitioning.filtered1,Smallworldness.filtered1)
    
    
    
    
    
    
    #CONNECTIVITY measures reflect degree to which network differs from complete network
    #Edge density: % of edges compared to maximum
    #Average degree: Avg. number of links
    #Average path length: avg of shortest pathes between reachable nodes
    #Network diameter: longest of shortest paths
    
    #CENTRALITY measures quantify heterogeneity in network structure
    #Average clustering coefficient:
    #Components
    
    
    
    
    
    
    
    
    
    
    
    #V(showme.looped)$size<-80
#     netname<-paste(name2,"Network_filtered1.pdf",sep='_')
#     
#     setwd(location)
#     pdf(netname,width= 12, height= 12,family="NimbusRom")
#     op <- par(mfrow=c(1,1))# rows, columns
#     par(mar=c(0,0,0,0)) 
    # plot(clp, showme.looped,edge.curved=.1,edge.arrow.size = 0.45,vertex.size=18,layout=l,rescale=FALSE,
    #      xlim=range(l[,1]),ylim=range(l[,2]))
    
    
    # plot(clp, showme.looped,edge.curved=.1,edge.arrow.size = 0.45,vertex.size=18,layout=l,rescale=TRUE,
    #      xlim=c(-1,1),ylim=c(-1,1))
    
 #   V(showme.looped)$raster <- behaviorplotlist
  #  V(showme)$raster <- behaviorplotlist
    
    
    
    
    # print(  plot(clp, showme.looped,edge.curved=.1,layout=l,rescale=TRUE,vertex.shape="raster",
    #        vertex.label=NA, 
    #        #vertex.size=15*(V(showme.looped)$selfsize), # resizes vertices based on self-self transitions
    #        vertex.size=25,
    #        edge.color=col.tr,
    #        edge.arrow.size=0.35,
    #        
    #        edge.width=3,
    #        edge.arrow.width=1.6,
    #        edge.loop.angle=1.5,
    #        
    #        xlim=c(-1,1),ylim=c(-1,1)))
    
#     print(  plot(clp, showme,edge.curved=.1,layout=l,rescale=TRUE,vertex.shape="raster",
#                  vertex.label=NA, 
#                  #vertex.size=15*(V(showme.looped)$selfsize), # resizes vertices based on self-self transitions
#                  vertex.size=25,
#                  edge.color=col.tr,
#                  edge.arrow.size=0.35,
#                  
#                  edge.width=3,
#                  edge.arrow.width=1.6,
#                  edge.loop.angle=1.5,
#                  
#                  xlim=c(-1,1),ylim=c(-1,1)))
#     
#     print(text(0.9,-0.75,paste("Timescored(s) = ",time.value,sep=""),col="black",cex=0.75,adj=(0)))
#     print(text(0.9,-0.8,paste("behaviors = ",uniquebehaviors.filtered1,sep=""),col="black",cex=0.75,adj=(0)))  
#     print(text(0.9,-0.85,paste("mean deg = ",round(mean.degree.filtered1,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
#     print(text(0.9,-0.9,paste("density score = ",round(densityscore.filtered1,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
#     print(text(0.9,-0.95,paste("mean path L = ",round(mean_path_length.filtered1,digits=2),sep=""),col="black",cex=0.75,adj=(0)))  
#     print(text(0.9,-1,paste("network.diameter = ",network.diameter.filtered1,sep=""),col="black",cex=0.75,adj=(0)))  
#     print(text(0.9,-1.05,paste("Smallworldness = ",Smallworldness.filtered1,sep=""),col="black",cex=0.75,adj=(0)))
    
    #tkplot(clp, showme.looped)
    #####  par("usr")
    
    # library(plotrix)
    # 
    # scalelist<-list()
    # gettinsmaller<-seq(from=1,to=0.52, by=-0.02)
    # for (vc in 1:25){
    # 
    #   
    #   gettingsmaller<-gettinsmaller[vc]
    #   pointlocations<-cbind(rescale(l[,1],gettingsmaller*c(-1,1)),rescale(l[,2],gettingsmaller*c(-1,1)))
    #   par(new=TRUE)
    #   plot(pointlocations,xlab="", ylab="",xlim=c(-1,1),ylim=c(-1,1),pch=vc)
    # 
    # }
    # dev.off()
    
    
    
    # for(bbb in 1:uniquebehaviors){
    #   rasterImage(behaviorplotlist[[bbb]], xleft = pointlocations[bbb,1], xright=pointlocations[bbb,1]+0.15,
    #               ybottom=pointlocations[bbb,2], ytop=pointlocations[bbb,2]+0.15)  
    #   par(new=TRUE)
    #   plot(pointlocations,xlab="", ylab="",xlim=c(-1,1),ylim=c(-1,1))
    # }
    # 
    # pointlocations<-cbind(rescale(l[,1],0.92*c(-1,1)),rescale(l[,2],0.92*c(-1,1)))
    # pointlocations<-pointlocations
    # 
    # for(bbb in 1:uniquebehaviors){
    #   rasterImage(behaviorplotlist[[bbb]], xleft = pointlocations[bbb,1], xright=pointlocations[bbb,1]+0.15,
    #               ybottom=pointlocations[bbb,2], ytop=pointlocations[bbb,2]+0.15)  
    #   par(new=TRUE)
    #   
    # }
    
    # dev.off()
    
    #Another network program
    #geph
  } else {
    c.no.OFFs<-as.character(no.OFFs.filtered1)
    TransitionMatrix.cum.filtered1<-  createSequenceMatrix(c.no.OFFs,sanitize=FALSE)
    
    densityscore.filtered1<-mean.degree.filtered1<-mean_path_length.filtered1<-network.diameter.filtered1<-avg.cluster.coeff.filtered1<-iSimpson.filtered1<-relativebehavioraldiversity.filtered1<-iSimpsonT.filtered1<-relativetransitiondiversity.filtered1<-PropSelfTransitioning.filtered1<-Smallworldness.filtered1<-NA
    
    printme<-data.frame(time.value,
                        uniquebehaviors,densityscore,mean.degree,mean_path_length,network.diameter,avg.cluster.coeff,
                        iSimpson,relativebehavioraldiversity,iSimpsonT,relativetransitiondiversity,
                        PropSelfTransitioning,Smallworldness,
                        
                        uniquebehaviors.filtered1,densityscore.filtered1,mean.degree.filtered1,mean_path_length.filtered1,network.diameter.filtered1,avg.cluster.coeff.filtered1,
                        iSimpson.filtered1,relativebehavioraldiversity.filtered1,iSimpsonT.filtered1,relativetransitiondiversity.filtered1,
                        PropSelfTransitioning.filtered1,Smallworldness.filtered1)
  }
  
  

  
  
  ######################################################################################################################################  
  ######################################################################################################################################
  
  

  
  sevencomponents<-list(printme,TransitionMatrix.cum,TransitionMatrix.cum.filtered1,behaviorcounts,behaviorcounts.filtered1,six,seven)
  
  return(sevencomponents)
}


##############
## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#NGRAM AND ZIPF's LAW FNX
ZipNgram<-function(z,plotzips=FALSE){
  #N-gram part
  simpstring<-preprocess(concatenate(z))
    for(jk in 1:10){
      ng<-ngram(simpstring,n=jk)
      print(ng)
    }
  
  
  summarytable<-table(z)
  proportiontable<-prop.table(summarytable)
  uniquebehaviors<-length(summarytable)
  
  
  
  
  
  
  l.freqs<-log10(c(summarytable))
  l.ranks<-log10(c(rank(-summarytable,ties="random")))
  zipfreq<-lm(l.freqs~l.ranks)
  zipdata<-summary(zipfreq) 
  zipSlope<-zipdata$coefficients[2,1]
  
  if(plotzips==TRUE){ 
    plot(l.freqs~l.ranks)
    abline(zipfreq,lwd=2,col='blue')
  }
  
  return(zipSlope)
  
}

#Zipf's brevity
# aaa is an transitionmatrix, or composite character string
ZipBrevity<-function(aaa,spname,makeplots="no",behavior.or.words="words",
                     bootstrpped="no",dotcolor="gray",alphaval=0.1,dotoutlinealpha=0.2,
                     plotline="straight",
                     xxx=x.l,yyy=y.l,corrtype="pearson",plotiteration=1,plotbox="l"){
  
  if(class(aaa)=="matrix"){
    summ<-colSums(aaa)
    summ.df<-data.frame(summ)
  
  } 
  
  if(class(aaa)=="character"){
    summ<-table(aaa)
    summ.df<-data.frame(summ)
    row.names(summ.df)<-summ.df[,1]
    summ.df<-summ.df[,-1,drop=FALSE]
    colnames(summ.df)<-"summ"
  }
  
  
  compression<-row.names(summ.df)
  
  
  if(behavior.or.words=="behavior"){
  partslist<-paste(compression,".",sep='')
  units<-unlist(lapply(partslist,function(x) length(strsplit(x,".",fixed=TRUE)[[1]])))
  }
  
  if(behavior.or.words=="words"){
    units<-unlist(lapply(compression,function(x) nchar(x)))
  }
  
  summ.df$units<-units
  colnames(summ.df)[1]<-"Freq"
  
  if(makeplots=="yes"){
    if(bootstrpped=="yes"){
      
      
      #if plotiteration = 1, plot with axes, otherwise, do not replot axes
      if(plotiteration==1){
        plot(summ.df$Freq~jitter(summ.df$units,amount=.08),#main=spname,
             ylab="Frequency",xlab=paste("# of ", behavior.or.words," elements",sep=''),
             xlim=xxx,ylim=yyy,bty=plotbox,
             pch=21,bg=add.alpha(dotcolor,alphaval),col=add.alpha("blue",dotoutlinealpha),cex=1)
        title(spname,adj=0,col=add.alpha("darkgreen",.7))
      } else {
        plot(summ.df$Freq~jitter(summ.df$units,amount=.08),#main=spname,
             ylab='',xlab='',
             xlim=xxx,ylim=yyy,bty=plotbox,
             axes=FALSE,
             pch=21,bg=add.alpha(dotcolor,alphaval),col=add.alpha("blue",dotoutlinealpha),cex=1)
      }
      

    if(plotline=="loess"){
      loess.model <- loess(summ.df$units~summ.df$Freq,
                           span=.6,degree=0)
      L1<-predict(loess.model)
      
      ntem<-cbind(L1[order(summ.df$Freq)],summ.df$Freq[order(summ.df$Freq)])
      #ntem[is.na(ntem)]<-0
      lines(ntem[,1]~ntem[,2],#col="blue")
            col=add.alpha(dotcolor,.9))
      
    }
    if(plotline=="straight"){
      abline(lm(summ.df$Freq~summ.df$units),col=add.alpha(dotcolor,.9))
    }
    par(new=TRUE) 
      
    } else {
      plot(summ.df$Freq~jitter(summ.df$units,amount=.08),main=spname,
           ylab="Frequency",xlab=paste("# of ", behavior.or.words," elements",sep=''),
           bty=plotbox,
           pch=21,bg=add.alpha(dotcolor,alphaval),col=add.alpha("blue",dotoutlinealpha),cex=1)
      abline(lm(summ.df$Freq~summ.df$units),col=add.alpha(dotcolor,.9))
      
    }
    
    
  }
  
  library(Hmisc)
  xj<-rcorr(summ.df$units,summ.df$Freq, type=corrtype) # type can be pearson or spearman
  info<-t(data.frame(c(xj[[1]][1,2],xj[[3]][1,2])))
  colnames(info)<-c(corrtype,"pval")
  rownames(info)<-NULL
  info<-data.frame(info)
  info$sig<-ifelse(info$pval<0.05,"yes","no")
  info$totbeh<-sum(summ.df$Freq)
  info$uniqbev<-length(compression)
  ab<-summary(lm(summ.df$Freq~summ.df$units))
  info2<-cbind(info,t(data.frame(ab$coefficients[1,])),t(data.frame(ab$coefficients[2,])))
  colnames(info2)[c(6:13)]<-c("Int","Int.SE","Int.t","Int.p","Slope","Slope.SE","Slope.t","Slope.p")
  
  return(info2)
}

#Zipf's brevity, mean
#as above, but returns mean y-value (Frequency) for each x-value (# elements)
ZipBrevity.mean<-function(aaa,spname,makeplots="no",behavior.or.words="words",
                     bootstrpped="no",dotcolor="gray",alphaval=0.1,dotoutlinealpha=0.2,
                     plotline="straight",
                     xxx=x.l,yyy=y.l,corrtype="pearson",plotiteration=1,plotbox="l"){
  
  if(class(aaa)=="matrix"){
    summ<-colSums(aaa)
    summ.df<-data.frame(summ)
    
  } 
  
  if(class(aaa)=="character"){
    summ<-table(aaa)
    summ.df<-data.frame(summ)
    row.names(summ.df)<-summ.df[,1]
    summ.df<-summ.df[,-1,drop=FALSE]
    colnames(summ.df)<-"summ"
  }
  
  
  compression<-row.names(summ.df)
  
  
  if(behavior.or.words=="behavior"){
    partslist<-paste(compression,".",sep='')
    units<-unlist(lapply(partslist,function(x) length(strsplit(x,".",fixed=TRUE)[[1]])))
  }
  
  if(behavior.or.words=="words"){
    units<-unlist(lapply(compression,function(x) nchar(x)))
  }
  
  summ.df$units<-units
  colnames(summ.df)[1]<-"Freq"
  
  gruped.data<-as.data.frame(summ.df %>%
                  group_by(units) %>%
                    summarise(F.mean=mean(Freq)))
  
  
  if(makeplots=="yes"){
    if(bootstrpped=="yes"){
      
      
      #if plotiteration = 1, plot with axes, otherwise, do not replot axes
      if(plotiteration==1){
        plot(gruped.data$F.mean~jitter(gruped.data$units,amount=.00),#main=spname,
             ylab="Frequency",xlab=paste("# of ", behavior.or.words," elements",sep=''),
             xlim=xxx,ylim=yyy,bty=plotbox,
             pch=21,bg=add.alpha(dotcolor,alphaval),col=add.alpha("blue",dotoutlinealpha),cex=1)
        title(spname,adj=0,col=add.alpha("darkgreen",.7))
      } else {
        plot(gruped.data$F.mean~jitter(gruped.data$units,amount=.00),#main=spname,
             ylab='',xlab='',
             xlim=xxx,ylim=yyy,bty=plotbox,
             axes=FALSE,
             pch=21,bg=add.alpha(dotcolor,alphaval),col=add.alpha("blue",dotoutlinealpha),cex=1)
      }
      
      
      if(plotline=="loess"){
        loess.model <- loess(summ.df$units~summ.df$Freq,
                             span=.6,degree=0)
        L1<-predict(loess.model)
        
        ntem<-cbind(L1[order(summ.df$Freq)],summ.df$Freq[order(summ.df$Freq)])
        #ntem[is.na(ntem)]<-0
        lines(ntem[,1]~ntem[,2],#col="blue")
              col=add.alpha(dotcolor,.9))
        
      }
      if(plotline=="straight"){
        abline(lm(gruped.data$F.mean~gruped.data$units),col=add.alpha(dotcolor,.9))
      }
      par(new=TRUE) 
      
    } else {
      plot(gruped.data$F.mean~jitter(gruped.data$units,amount=.00),main=spname,
           ylab="Frequency",xlab=paste("# of ", behavior.or.words," elements",sep=''),
           bty=plotbox,
           pch=21,bg=add.alpha(dotcolor,alphaval),col=add.alpha("blue",dotoutlinealpha),cex=1)
      abline(lm(gruped.data$F.mean~gruped.data$units),col=add.alpha(dotcolor,.9))
      
    }
    
    
  }
  
  library(Hmisc)
  xj<-rcorr(summ.df$units,summ.df$Freq, type=corrtype) # type can be pearson or spearman
  info<-t(data.frame(c(xj[[1]][1,2],xj[[3]][1,2])))
  colnames(info)<-c(corrtype,"pval")
  rownames(info)<-NULL
  info<-data.frame(info)
  info$sig<-ifelse(info$pval<0.05,"yes","no")
  info$totbeh<-sum(summ.df$Freq)
  info$uniqbev<-length(compression)
  
  return(list(info,gruped.data))
}

#Zipf's brevity, as above, except with axes flipped (x-axis = frequency, y-axis= # elements)
# aaa is an transitionmatrix, or composite character string
ZipBrevity.dolphin<-function(aaa,spname,makeplots="no",behavior.or.words="words",
                     bootstrpped="no",dotcolor="gray",alphaval=0.1,dotoutlinealpha=0.2,
                     plotline="straight",
                     xxx=x.l,yyy=y.l,corrtype="pearson",plotiteration=1){
  
  if(class(aaa)=="matrix"){
    summ<-colSums(aaa)
    summ.df<-data.frame(summ)
    
  } 
  
  if(class(aaa)=="character"){
    summ<-table(aaa)
    summ.df<-data.frame(summ)
    row.names(summ.df)<-summ.df[,1]
    summ.df<-summ.df[,-1,drop=FALSE]
    colnames(summ.df)<-"summ"
  }
  
  
  compression<-row.names(summ.df)
  
  
  if(behavior.or.words=="behavior"){
    partslist<-paste(compression,".",sep='')
    units<-unlist(lapply(partslist,function(x) length(strsplit(x,".",fixed=TRUE)[[1]])))
  }
  
  if(behavior.or.words=="words"){
    units<-unlist(lapply(compression,function(x) nchar(x)))
  }
  
  summ.df$units<-units
  colnames(summ.df)[1]<-"Freq"
  
  if(makeplots=="yes"){
    if(bootstrpped=="yes"){
      
      #if plotiteration = 1, plot with axes, otherwise, do not replot axes
      if(plotiteration==1){
        plot(jitter(summ.df$units,amount=.08)~summ.df$Freq,main=spname,
           ylab=paste("# of ", behavior.or.words," elements",sep=''),xlab="Frequency",
           xlim=yyy,ylim=xxx,
           pch=21,bg=add.alpha(dotcolor,alphaval),col=add.alpha("blue",dotoutlinealpha),cex=2)
      } else {
        plot(jitter(summ.df$units,amount=.08)~summ.df$Freq,main=spname,
             ylab='',xlab='',
             xlim=yyy,ylim=xxx,
             axes=FALSE,
             pch=21,bg=add.alpha(dotcolor,alphaval),col=add.alpha("blue",dotoutlinealpha),cex=2)
      }
      
      
      
      
      if(plotline=="loess"){
        loess.model <- loess(summ.df$units~summ.df$Freq,
                             span=.6,degree=0)
        L1<-predict(loess.model)
        
        ntem<-cbind(L1[order(summ.df$Freq)],summ.df$Freq[order(summ.df$Freq)])
        #ntem[is.na(ntem)]<-0
        lines(ntem[,1]~ntem[,2],#col="blue")
              col=add.alpha(dotcolor,.9))
        
      }
      if(plotline=="straight"){
        abline(lm(summ.df$units~summ.df$Freq),col=add.alpha(dotcolor,.9))
      }
      par(new=TRUE) 
      
    } else {
      plot(jitter(summ.df$units,amount=.08)~summ.df$Freq,main=spname,
           ylab=paste("# of ", behavior.or.words," elements",sep=''),xlab="Frequency",
           pch=21,bg=add.alpha(dotcolor,alphaval),col=add.alpha("blue",dotoutlinealpha),cex=2)
      abline(lm(summ.df$units~summ.df$Freq),col=add.alpha(dotcolor,.9))
      
      
    }
    
    
  }
  
  library(Hmisc)
  xj<-rcorr(summ.df$Freq,summ.df$units, type=corrtype) # type can be pearson or spearman
  info<-t(data.frame(c(xj[[1]][1,2],xj[[3]][1,2])))
  colnames(info)<-c(corrtype,"pval")
  rownames(info)<-NULL
  info<-data.frame(info)
  info$sig<-ifelse(info$pval<0.05,"yes","no")
  info$totbeh<-sum(summ.df$Freq)
  info$uniqbev<-nrow(aaa)
  
  return(info)
}



#Zipf's brevity, as above, but with bootstrapping implemented in the function (rather than multiple calls to the ZipBrevity function)
ZipBrevity.Boot<-function(aaa=completestringforspecies,howmanyboots=nbooties,
                          sizeofeachsample=subsampleN,
                          spnameB="bigbird",
                          makeplots="yes",behavior.or.words="behavior",
                          bootstrpped="yes",dotcolor=use.color,alphaval=0.2,
                          dotoutlinealpha=0.05,
                          #plotline="loess",
                          plotline="straight",
                          xxx.B=x.yyy,yyy.B=y.xxx,corrtype=correlation.type,
                          plotbox="l"){
  
  for(u in 1:howmanyboots){
    oneuse<-sample(aaa,size=sizeofeachsample,replace = FALSE)
    
    plot.flag<-2 #plot.flag prevents subsequent plots from re-plotting axes (unnecessarily)
    if(u==1)
      plot.flag<-1
    

    thissample<-ZipBrevity(oneuse,spname=spnameB,makeplots=makeplots,behavior.or.words=behavior.or.words,
                           bootstrpped="yes",dotcolor=dotcolor,alphaval=alphaval,dotoutlinealpha=dotoutlinealpha,
                           plotline=plotline,
                           xxx=xxx.B,yyy=yyy.B,corrtype=corrtype,plotiteration=plot.flag,plotbox="l")
    
    if(u==1){
      bootdf<-thissample
    } else {
      bootdf<-rbind(bootdf,thissample)
    }
  }
  
  
  ok.plot<-function(){
 
    par(fig=c(0,1,0,1))
    par(mar=c(5,5,1,1)) #Margines of each plot (bottom, left, top, right)
    
    for(u in 1:howmanyboots){
      oneuse<-sample(aaa,size=sizeofeachsample,replace = FALSE)
      
      plot.flag<-2 #plot.flag prevents subsequent plots from re-plotting axes (unnecessarily)
      if(u==1)
        plot.flag<-1
      
      
      thissample<-ZipBrevity(oneuse,spname=spnameB,makeplots="yes",behavior.or.words=behavior.or.words,
                             bootstrpped="yes",dotcolor=dotcolor,alphaval=alphaval,dotoutlinealpha=dotoutlinealpha,
                             plotline=plotline,
                             xxx=xxx.B,yyy=yyy.B,corrtype=corrtype,plotiteration=plot.flag,plotbox="l")
      

    }
    
    par(fig=c(0.45,.95,0.7,.95), new=TRUE)
    par(mar=c(3,3,1,1)) #Margines of each plot (bottom, left, top, right)
    addthemhistograms(bootdf$spearman,1.75,1.75,xlabel="Correlation coefficient",
                      xlimit=c(-1,1),barcolors=rgb(0,1,1,0.4),vertbar=0,vertcol="darkgreen")
    
    par(fig=c(0.45,.95,0.45,.70), new=TRUE)
    par(mar=c(3,3,1,1)) #Margines of each plot (bottom, left, top, right)
    addthemhistograms(-log(bootdf$pval),1.75,1.75,xlabel="-log(p-value)",
                      xlimit=c(0,25),barcolors=rgb(1,1,0,0.4),vertbar=2.995732,vertcol="red")
  }
  
  #plotme<-ok.plot()
  
  plotvalues<-c(howmanyboots,sizeofeachsample,spnameB,dotcolor,xxx.B,yyy.B)
  
  
  return(list(bootdf,ok.plot))
  #resetPar() # reset to defaults (not necessarilly current values)
  #par(mfrow=c(8,4))# rows, columns
}

#Zipf's brevity, as above, but with bootstrapping implemented in the function 
#(rather than multiple calls to the ZipBrevity function)
# and it calls ZipBrevity.mean (not ZipBrevity)
ZipBrevity.Boot.mean<-function(aaa=completestringforspecies,howmanyboots=nbooties,
                          sizeofeachsample=subsampleN,
                          spnameB="bigbird",
                          makeplots="yes",behavior.or.words="behavior",
                          bootstrpped="yes",dotcolor=use.color,alphaval=0.2,
                          dotoutlinealpha=0.05,
                          #plotline="loess",
                          plotline="straight",
                          xxx.B=x.yyy,yyy.B=y.xxx,corrtype=correlation.type,
                          plotbox="l"){
  specific.values<-list()
  for(u in 1:howmanyboots){
    oneuse<-sample(aaa,size=sizeofeachsample,replace = FALSE)
    
    plot.flag<-2 #plot.flag prevents subsequent plots from re-plotting axes (unnecessarily)
    if(u==1)
      plot.flag<-1
    
    
    listsample<-ZipBrevity.mean(oneuse,spname=spnameB,makeplots=makeplots,behavior.or.words=behavior.or.words,
                           bootstrpped="yes",dotcolor=dotcolor,alphaval=alphaval,dotoutlinealpha=dotoutlinealpha,
                           plotline=plotline,
                           xxx=xxx.B,yyy=yyy.B,corrtype=corrtype,plotiteration=plot.flag,plotbox="l")
    thissample<-listsample[[1]]
    specific.values[[u]]<-listsample[[2]]
    
    if(u==1){
      bootdf<-thissample
    } else {
      bootdf<-rbind(bootdf,thissample)
    }
  }
  
  
  ok.plot<-function(){
    
    par(fig=c(0,1,0,1))
    par(mar=c(5,5,1,1)) #Margines of each plot (bottom, left, top, right)
    
    for(u in 1:howmanyboots){
      oneuse<-sample(aaa,size=sizeofeachsample,replace = FALSE)
      
      plot.flag<-2 #plot.flag prevents subsequent plots from re-plotting axes (unnecessarily)
      if(u==1)
        plot.flag<-1
      
      
      thissample<-ZipBrevity.mean(oneuse,spname=spnameB,makeplots="yes",behavior.or.words=behavior.or.words,
                             bootstrpped="yes",dotcolor=dotcolor,alphaval=alphaval,dotoutlinealpha=dotoutlinealpha,
                             plotline=plotline,
                             xxx=xxx.B,yyy=yyy.B,corrtype=corrtype,plotiteration=plot.flag,plotbox="l")
      
      
    }
    
    par(fig=c(0.45,.95,0.7,.95), new=TRUE)
    par(mar=c(3,3,1,1)) #Margines of each plot (bottom, left, top, right)
    addthemhistograms(bootdf$spearman,1.75,1.75,xlabel="Correlation coefficient",
                      xlimit=c(-1,1),barcolors=rgb(0,1,1,0.4),vertbar=0,vertcol="darkgreen")
    
    par(fig=c(0.45,.95,0.45,.70), new=TRUE)
    par(mar=c(3,3,1,1)) #Margines of each plot (bottom, left, top, right)
    addthemhistograms(-log(bootdf$pval),1.75,1.75,xlabel="-log(p-value)",
                      xlimit=c(0,25),barcolors=rgb(1,1,0,0.4),vertbar=2.995732,vertcol="red")
  }
  
  #plotme<-ok.plot()
  
  plotvalues<-c(howmanyboots,sizeofeachsample,spnameB,dotcolor,xxx.B,yyy.B)
  
  
  return(list(bootdf,ok.plot,specific.values))
  #resetPar() # reset to defaults (not necessarilly current values)
  #par(mfrow=c(8,4))# rows, columns
}

#Zipf's max info
ZipMaxInfo<-function(aaa,spname,behavior.or.words="words"){
  
  if(class(aaa)=="matrix"){
    summ<-colSums(aaa)
    summ.df<-data.frame(summ)
    
  } 
  
  if(class(aaa)=="character"){
    summ<-table(aaa)
    summ.df<-data.frame(summ)
    row.names(summ.df)<-summ.df[,1]
    summ.df<-summ.df[,-1,drop=FALSE]
    colnames(summ.df)<-"summ"
  }
  
  
  compression<-row.names(summ.df)
  
  
  if(behavior.or.words=="behavior"){
    partslist<-paste(compression,".",sep='')
    units<-unlist(lapply(partslist,function(x) length(strsplit(x,".",fixed=TRUE)[[1]])))
  }
  
  if(behavior.or.words=="words"){
    units<-unlist(lapply(compression,function(x) nchar(x)))
  }
  
  summ.df$units<-units
  colnames(summ.df)[1]<-"Freq"

  max.f.u<-c(max(summ.df$Freq),max(summ.df$units))
  return(max.f.u)
}

#############
#CODE TO ADD HISTOGRAMS TO EXISTING PLOTS
addthemhistograms<-function(x,xl,yl,xlabel,xlimit=c(-1,1),barcolors=rgb(0,1,1,0.4),vertbar=0,vertcol="darkgreen"){
  hist(x,main='',col=barcolors,xlab="",ylab="",tck=-.02,xlim=xlimit,axes=FALSE)
  axis(1,padj = -.65)
  axis(2,padj = .65)
  abline(v=vertbar,col=vertcol,lty=2,lwd=2)
  title(ylab="Freq",line=xl,cex=.8)
  title(xlab=xlabel,line=yl,cex=.8)
}







#Behavioral transition complexity measurer
# takes an input square matrix of behavioral transitions, where each row/column correspondings to a composite
# behavior (e.g. 05.BP1) and returns a matrix where every cell value now has an integer value for the
# 'complexity' of the transition between those two behavioral states, corresponding to how many behavioral states
# needed to change to facilitate the transition from compositebehavior1 to compositebehavior2
Alignment<-read.csv('C:/Users/Rusty/Amazon Drive/BOP/GeneSequences/BehavPosition.csv')
Codes<-c(rep(c(1,0),32))
BehavioralTransitionComplexity<-function(behavioraltransitionmatrix){
  filler<-matrix(nrow=nrow(behavioraltransitionmatrix),ncol=31,0)
  #need.strip<-cbind(row.names(behavioraltransitionmatrix),filler)
  need.strip<-filler
  row.names(need.strip)<-row.names(behavioraltransitionmatrix)
  #colnames(need.strip)<-c("orig",t(Alignment)[1,])
  colnames(need.strip)<-c(t(Alignment)[1,])
  
  for(m in 1:ncol(need.strip)){
    ppp<-m*2-1
    need.strip[,m]<-ifelse(grepl(colnames(need.strip)[m],row.names(need.strip)),Codes[ppp],Codes[ppp+1])
  }
  need.strip<-as.data.frame(need.strip)
 
  for(w in 1:nrow(behavioraltransitionmatrix)){
    ref.behav<-row.names(behavioraltransitionmatrix)[w]
    ref.string<-need.strip[which(row.names(need.strip)==ref.behav),]
    
    for(x in 1:nrow(behavioraltransitionmatrix)){
      comparison.behav<-row.names(behavioraltransitionmatrix)[x]
      comparison.string<-need.strip[which(row.names(need.strip)==comparison.behav),]
      behavioraldifferences<-sum(abs(ref.string-comparison.string))#sum of absolute differences between behaviors
      
      if(x==1){
        thiscolumn<-behavioraldifferences
      } else {
        thiscolumn<-rbind(thiscolumn,behavioraldifferences)
      }
      
    }
    
    if(w==1){
      completecomparisons<-thiscolumn
    } else {
      completecomparisons<-cbind(completecomparisons,thiscolumn)
    }
  }
  row.names(completecomparisons)<-colnames(completecomparisons)<-row.names(behavioraltransitionmatrix)
  
  return(completecomparisons)
}

#creates matrix with 'named' transitions, based on rownames and colnames of transition matrix
Transition.Namer<-function(thismatrix){
  for(nr in 1:nrow(thismatrix)){
    a.from<-rownames(thismatrix)[nr]
    for(nc in 1:ncol(thismatrix)){
      b.to<-colnames(thismatrix)[nc]
      transname<-paste(a.from,b.to,sep='-')
      
      if(nc==1){
        thisrowoftransitions<-transname
      } else {
        thisrowoftransitions<-c(thisrowoftransitions,transname)
      }
    }
    if(nr==1){
      namedts<-thisrowoftransitions
    } else {
      namedts<-rbind(namedts,thisrowoftransitions)
    }
  }
  return(namedts)
}


#BOOTSTRAPPER FOR COMPLEXITY ANALYSIS
transitional.Zips.Boots<-function(combovectorbig,actualprobmat,trans.to.use,transbooties,
                                  plot.trans.boots="yes",species,corrtype=correlation.type,
                                  dotcolor,alphaval,dotoutlinealpha,plotbox){
  
  if(trans.to.use>length(combovectorbig)){
    trans.to.use<-length(combovectorbig)
  }
  
  #quickly get a sense of reasonable x- and y-limits
  ###########################
  yy<-0
  xx<-0
  for(quicklim in 1:1000){
    rounder<-sample(combovectorbig,size=trans.to.use,replace=TRUE,prob=as.vector(actualprobmat))
    comparison.df<-data.frame(table(rounder))
    comparison.df$booted.trans.sample<-as.character(comparison.df$rounder)
    comparison.df$complx<-lapply(comparison.df$booted.trans.sample,function(x) {strsplit(x,"=")[[1]][2]})
    comparison.df$complx<-as.numeric(as.character(comparison.df$complx))
    
    max.info<-c(max(comparison.df$complx),max(comparison.df$Freq))
    
    if(max.info[2]>yy)
      yy<-max.info[2]
    
    if(max.info[1]>xx)
      xx<-max.info[1]
  }
  y.xxx<-c(0,xx)
  x.yyy<-c(0,yy)
  ########################
  allpulls<-list()
  for(u2 in 1:transbooties){
    #for each boot iteration, pull 
    booted.trans.sample<-sample(combovectorbig,size=trans.to.use,replace=TRUE,prob=as.vector(actualprobmat))
    comparison.df<-data.frame(table(booted.trans.sample))
    comparison.df$booted.trans.sample<-as.character(comparison.df$booted.trans.sample)
    comparison.df$complx<-lapply(comparison.df$booted.trans.sample,function(x) {strsplit(x,"=")[[1]][2]})
    comparison.df$complx<-as.numeric(as.character(comparison.df$complx))
    
    #correlation analysis requires 5 more values
    if(nrow(comparison.df)>4){
      library(Hmisc)
      xj<-rcorr(comparison.df$complx,comparison.df$Freq, type=corrtype) # type can be pearson or spearman
      info<-t(data.frame(c(xj[[1]][1,2],xj[[3]][1,2])))
      colnames(info)<-c(corrtype,"pval")
      rownames(info)<-NULL
      info<-data.frame(info)
      info$sig<-ifelse(info$pval<0.05,"yes","no")
      info$total.trans<-sum(comparison.df$Freq)
      info$uniquetrans<-nrow(comparison.df)
    } else {
      if(exists("info")){
       info2<-info  
      } else {
       info2<-t(as.data.frame(c(NA,NA,NA,NA,NA)))
       rownames(info2)<-NULL
       colnames(info2)<-c("spearman","pval","sig","total.trans","uniquetrans")
      }
       
      
      info2$spearman<-NA
      info2$pval<-NA
      info<-info2
    }
    
    allpulls[[u2]]<-comparison.df
    
    
    if(u2==1){
      
      correlation.summaries<-info
      
      if(plot.trans.boots=="yes"){
         plot(comparison.df$Freq~jitter(comparison.df$complx,amount=0.08),
             ylab="Frequency",xlab="Transitional complexity",
             xlim=y.xxx,ylim=x.yyy,bty=plotbox,
             pch=21,bg=add.alpha(dotcolor,alphaval),col=add.alpha("blue",dotoutlinealpha),cex=1)
        title(species,adj=0,col=add.alpha("darkgreen",.7))
        
        #if all zeros complexities, plot horizontal line
        if(!is.na(info$spearman)){
          abline(lm(comparison.df$Freq~comparison.df$complx),col=add.alpha(dotcolor,.9))
        } else {
          abline(h=0,col=add.alpha(dotcolor,.9))
        }
          

      }
      
    } else {
      
      correlation.summaries<-rbind(correlation.summaries,info)
      
      
      if(plot.trans.boots=="yes"){
        par(new=TRUE)
        
        plot(comparison.df$Fre~jitter(comparison.df$complx,amount=.08),#main=spname,
             ylab='',xlab='',
             xlim=y.xxx,ylim=x.yyy,bty=plotbox,
             axes=FALSE,
             pch=21,bg=add.alpha(dotcolor,alphaval),col=add.alpha("blue",dotoutlinealpha),cex=1)
        
        #if all zeros complexities, plot horizontal line
        if(!is.na(info$spearman)){
          abline(lm(comparison.df$Freq~comparison.df$complx),col=add.alpha(dotcolor,.9))
        } else {
          abline(h=0,col=add.alpha(dotcolor,.9))
        }
        
      }
    }
   
  }
  
  megainfo<-list(correlation.summaries,allpulls)
  
  return(megainfo)
  
}




#BOOTSTRAP.networksummarizer, which depends on 'networksummarizer' performs network summary, with bootstrap iterations
# to test robustness.
# Returns 'RemovalSummary' DF which gives a measure of 'robustness'
# (EffectiveGraphResistance)
#*******************************************
# Effective graph resistance measures and concept based on:
# Klein DJ, Randić M. Resistance distance. Journal of Mathematical Chemistry 1993; 12: 81–95.
# Ellens W, Kooij RE. 2013. Graph measures and network robustness. arXiv:1311.5064v1

# "Klein and Randić [30] found that the effective graph resistance of a connected network can be written as a function of all 
# non-zero Laplacian eigenvalues of the network." --- from Yang et al. 2016

# Yang et al. 2016. The Rationality of Four Metrics of Network Robustness: A Viewpoint of Robust Growth of 
# Generalized Meshes. PLOS ONE. https://doi.org/10.1371/journal.pone.0161077
#*************************************************************************************
# z = a nx1 matrix (1 column), containing behavioral observations
# n.bootstraps = n of repetitions per removal proportion
# pr.thresh = threshold for reduced data summary variable/unmanipulated summary variable 
# show. prog = logical (default TRUE) to print progress for each
BOOTSTRAP.networksummarizer<-function(z,n.bootstraps=100,
                                      #pr.thresh=0.25,which.var.ref="mean_path_length",
                                      show.prog=TRUE){
  
  
  
  OriginalFullNetwork<-networksummarizer(z)
  OriginalFullNetwork<-data.frame(1, OriginalFullNetwork);colnames(OriginalFullNetwork)[1]<-"HowMuch"
  
  remove.dis<-seq(0.01,0.99,.01)
  remove.dis<-c(remove.dis,0.995)
  
  
  #Bootstrap ---subset and recalculate
  md<-1
  
  pb<-txtProgressBar(min=0,max=n.bootstraps,char="/")
  
  while(md<length(remove.dis)){
    
    HowMuch<-1-remove.dis[md]
    
    
    
    
    for(xar in 1:n.bootstraps){
      #new.z<-sample(z,round(length(z)*HowMuch),replace=FALSE)
      
      samplestoreplace<-sample(length(z),length(z)-round(length(z)*HowMuch),replace=FALSE)
      
      tempuse<-z
      
      tempuse[samplestoreplace]<-NA
      new.z<-na.locf(tempuse,na.rm=FALSE)
      
      nnn<-as.matrix(new.z,ncol=1)
      
      SubNetwork<-networksummarizer(nnn)
      SubNetwork<-data.frame(HowMuch,SubNetwork)
      
      if(xar==1)
        SubFrame<-rbind(OriginalFullNetwork,SubNetwork)
      
      if(xar>1)
        SubFrame<-rbind(SubFrame,SubNetwork)
      
      
      if(xar==n.bootstraps){
        
        #relative<-t(matrix((unlist(apply(SubFrame, 1,function(x) x/OriginalFullNetwork))),ncol=n.bootstraps))
        #colnames(relative)<-colnames(SubFrame)
        
        Line.Item<-apply(SubFrame,2,function(x) median(x))
        SD.Lines<-((apply(SubFrame,2,sd)));names(SD.Lines)<-paste(names(SD.Lines),".SD",sep='')
        Line.Item2<-c(Line.Item,SD.Lines[c(2:29)])
      }
    }
    
    #Line.Item2
    
    #RatioForRemoved<-data.frame(Line.Item/OriginalFullNetwork)
    
    if(md==1)
      RemovalSummary<-Line.Item2
    
    if(md>1)
      RemovalSummary<-rbind(RemovalSummary,Line.Item2)
    
    md<-md+1
    #IF ratio for defined variable falls below threshold, this will stop the loop (if pr.thresh is active)
    
    
    # if(RatioForRemoved[,which.var.ref]<pr.thresh || is.na(RatioForRemoved[,which.var.ref])){
    #   md<-100
    #   print("run-done")
    # }
    
    
    if(show.prog==TRUE)
      setTxtProgressBar(pb, md)
    
  }
  
  return(RemovalSummary)
  
}
