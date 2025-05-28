library(scales)

pal37 = c('208 210 85',
         '100 100 80',
         '74 73 19',
         '52 18 8',
         '45 12 11',
         '147 34 30',
         '223 58 37',
         '242 179 85',
         '250 224 80',
         '196 41 44',
         '128 46 92',
         '129 68 122',
         '132 74 152',
         '96 44 134',
         '72 29 107',
         '60 11 58',
         '94 160 186',
         '32 33 69',
         '80 113 141',
         '45 22 42',
         '28 37 79',
         '31 64 114',
         '26 64 53',
         '21 51 18',
         '17 45 14',
         '103 149 81',
         '150 206 94',
         '89 79 56',
         '46 52 36',
         '179 216 99',
         '14 34 12',
         '120 149 45',
         '223 234 77',
         '157 74 91',
         '104 91 209',
         '74 152 200',
         '40 99 175')
length(pal37)
pal37H = sapply(strsplit(pal37, " "), function(x)
  rgb(x[1], x[2], x[3], maxColorValue=255))

#show_col(pal37H)



pal2 = c('186 137 65',
         '161 116 61',
         '156 125 84',
         '45 67 114',
         '60 81 121',
         '52 119 182',
         '17 48 96',
         '107 108 145',
         '82 102 145',
         '117 155 214',
         '209 48 99',
         '226 83 76',
         '239 155 80',
         '234 177 72',
         '214 49 106',
         '169 40 33',
         '191 90 39',
         '150 60 35',
         '165 51 77',
         '87 39 145',
         '80 41 133',
         '91 33 115',
         '140 61 153',
         '162 64 166',
         '138 203 237',
         '59 135 199',
         '34 77 137',
         '13 34 83',
         '148 33 20',
         '83 50 25',
         '227 110 46',
         '234 152 57',
         '241 194 75',
         '122 75 130')
length(pal2)

pal34H = sapply(strsplit(pal2, " "), function(x)
  rgb(x[1], x[2], x[3], maxColorValue=255))

#show_col(pal34H)

#' Adds transparency to colour
#'
#' @param cols Vector of colours.
#' @param alphas Single value or vector of alphas
#' @param ... Passed to rgb
#' @return rgb colours with transparency set.
colAlpha = function(cols,alphas,...) {
  if(length(alphas)==1)
    alphas = rep(alphas,length(cols))
  tmp = col2rgb(cols)
  sapply(seq_len(ncol(tmp)),function(e) rgb(tmp[1,e],tmp[2,e],tmp[3,e],alphas[e]*255,maxColorValue=255,...))
}


col25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  #"black", 
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)



col22=c("#9cb169","#a75acb","#59b648","#636edd","#99b534","#d5439a","#51bf7f","#d73f52","#4cc8c6","#cd5d2a","#609dd8","#d0a63d","#5e63a9","#687327","#bd8dd9","#377b40","#9c4c89","#479e80","#df80ae","#a8793e","#a74a57","#e28674")


theme_classic_2 = theme_classic(base_size = 15) + theme(panel.border = element_rect(fill=F),axis.line = element_blank())


annotateBTB = function(btb.fp,PDID,minSegLen=1e6,subCl.minSegLen=1e7,tgtChrs=c(1:23,'X'),longFormat=T,removeBalancedSegs = FALSE,method = c('totalCN','allelicRatio')){
  chromInfo = read.delim('/lustre/scratch117/casm/team274/mt22/chrom_abspos_kb.txt',sep = '\t')
  #if(is.null(data) & !is.null(btb)){
  #  data = read.csv(btb,header=F,
  #                  col.names = c('Idx','Chr','Start','Stop','normTot','normMin','tumTot','tumMin'))  
  #}else if(is.null(data) & is.null(btb)){
  #  stop('Please specify the correct btb file path')
  #}else if (!is.null(data)){
  #  data=data
  #}
  
  #### Read in CN profile for both major and minor clone from Battenberg ####
  #PDID = unique(projMani$PDID)[1]
  #donorMani = projMani[projMani$PDID == PDID,]
  #btb.fp = unique(donorMani$battenbergFp)[!is.na(unique(donorMani$battenbergFp))]
  
  btb = gsub('summary.csv','subclones.txt.gz',btb.fp)
  btb = read.delim(btb,sep = '\t')
  btb$posID = paste0(btb$chr,'_',btb$startpos,'_',btb$endpos)
  btb$Idx = c(1:nrow(btb))
  
  # Keep CN segments on chr of interest only
  btb = btb[btb$chr %in% tgtChrs,]
  btb$chr = ifelse(btb$chr == 'X',23,as.numeric(btb$chr))
  
  # remove segments <= 1Mb
  btb$segLen = btb$endpos - btb$startpos + 1
  if(sum(btb$segLen <= minSegLen) >0 ){
    message(sprintf('Removing %d short segment for sample %s',sum(btb$segLen <= minSegLen),PDID))  
  }
  btb = btb[btb$segLen > minSegLen,]
  
  # Get major clone vs minor clone CN profile
  btb.majCl = tibble()
  btb.subCl = tibble()
  
  for(i in 1:nrow(btb)){
    if(is.na(btb$frac2_A[i])){
      tmp = btb[i,c('Idx','chr','startpos','endpos','nMaj1_A','nMin1_A','frac1_A','posID','segLen')]
      colnames(tmp) = c('Idx','Chr','Start','Stop','matNum','patNum','frac','posID','segLen')
      btb.majCl = rbind(btb.majCl,tmp)
    }else{
      if(btb$frac1_A[i] > btb$frac2_A[i]){
        tmp.maj = btb[i,c('Idx','chr','startpos','endpos','nMaj1_A','nMin1_A','frac1_A','posID','segLen')]
        tmp.min = btb[i,c('Idx','chr','startpos','endpos','nMaj2_A','nMin2_A','frac2_A','posID','segLen')]
      }else{
        tmp.maj = btb[i,c('Idx','chr','startpos','endpos','nMaj2_A','nMin2_A','frac2_A','posID','segLen')]
        tmp.min = btb[i,c('Idx','chr','startpos','endpos','nMaj1_A','nMin1_A','frac1_A','posID','segLen')]
      }
      colnames(tmp.maj) = c('Idx','Chr','Start','Stop','matNum','patNum','frac','posID','segLen')
      colnames(tmp.min) = c('Idx','Chr','Start','Stop','matNum','patNum','frac','posID','segLen')
      btb.majCl = rbind(btb.majCl,tmp.maj)
      btb.subCl = rbind(btb.subCl,tmp.min)
    }
  }
  
  
  
  btb.majCl$type = 'maj'
  btb.subCl$type = 'sub'
  message(sprintf('%s - Lowest frac for Major clone is %f',PDID,min(btb.majCl$frac)))
  message(sprintf('%s - Lowest frac for Sub clone is %f',PDID,min(btb.subCl$frac)))
  
  # Remove subclone CNA with < 0.1% cells and shorter than 50Mb
  message(sprintf('%s - Removing %d subclone CNA fragment due to low frac',PDID,sum(btb.subCl$frac<0.2)))
  message(sprintf('%s - Removing %d subclone CNA fragment due to short length < 20Mb',PDID,sum(btb.subCl$segLen<2e7)))
  message(sprintf('%s - Removing %d subclone CNA fragment due to short length < 50Mb',PDID,sum(btb.subCl$segLen<5e7)))
  btb.subCl = btb.subCl[btb.subCl$frac >=0.1 & btb.subCl$segLen >= subCl.minSegLen,]
  
  
  # Merge clonal CNprofile (if any) together
  all = rbind(btb.subCl,btb.majCl)
  all$tumTot = all$matNum+all$patNum
  all$tumFrac = all$matNum/(all$tumTot)
  all$tot2min = paste0(all$tumTot,':',all$patNum)
  all$newIdx = all$Idx
  
  # Process primary CNprofile (major clone)
  data = all[all$type == 'maj',]
  
  # Smoothing step:
  modSegs = data.frame()
  toRemove = c()
  
  for(chr in unique(data$Chr)){
    chromLen = chromInfo[chromInfo$chrom == chr & chromInfo$arm == 'q',]$end*1000
    chr.data = data[data$Chr == chr,]
    
    
    if(method=='totalCN'){ # to plot total copy number changes (eg. copyKat)
      for(config in unique(chr.data$tot2min)){
        cf.len.perc = sum(chr.data[chr.data$tot2min == config,]$segLen)*100/chromLen 
        #config = unique(chr.data[chr.data$tumFrac == tumFrac,]$config)
        message(sprintf('\nConfig %s for chr %s occupies %f perc',config,chr,cf.len.perc))
        
        if(cf.len.perc >= 90){
          message(sprintf('Merging config %s for chr %s',config,chr))
          breakFlag=1
          removeFlag=0
          toRemove = c(toRemove,unique(chr.data$Idx))
          
          tmp = chr.data[chr.data$tot2min == config,][1,]
          tmp$Idx = paste(chr.data[chr.data$tot2min == config,]$Idx,collapse = '_')
          tmp$newIdx = paste0(tmp$newIdx,'.', 1)
          tmp$Stop = chromLen
          tmp$Start = 1
          tmp$tumFrac = paste(unique(chr.data[chr.data$tot2min == config,]$tumFrac),collapse = '_') 
          
        }else if(cf.len.perc <= 10){ # remove if it's <10 % of Chr
          removeFlag=1
          message(sprintf('Removing config %s for chr %s',config, chr))
          tmp=NULL
          #tmp = data.frame(Idx=0,Chr=chr,Start=1,Stop=chromLen,normTot=2,normMin=1,tumTot = 2,tumMin = 1,segLen=chromLen, config='2:1',idx=0)
          breakFlag=0
        }else{
          # Merge fragments
          a = chr.data[chr.data$tot2min == config,]
          rmIdx.a = c()
          rmIdx.a.all = c()
          newSeg.a = tibble()
          
          if(nrow(a) > 1){
            for(f in 2:nrow(a)){
              if((a$Start[f] - a$Stop[f-1]) <= minSegLen){
                message('Merging fragments')
                rmIdx.a = c(rmIdx.a,a$Idx[c(f,(f-1))])
                rmIdx.a.all = c(rmIdx.a.all,a$Idx[c(f,(f-1))])
              }
              
              # Do the removal
              
              # If 2 fragments are far apart, or if this is the last segment --> stop and perform merger
              if(((a$Start[f] - a$Stop[f-1]) > minSegLen) || f == nrow(a)){
                if(length(rmIdx.a) > 1){
                  a.toRemove = a[a$Idx %in% unique(rmIdx.a),]
                  # Write new merged segment
                  segLen = (max(a.toRemove$Stop) - min(a.toRemove$Start)) + 1
                  patNum = as.numeric(sapply(strsplit(config,split = ':'),'[',2))
                  tumTot = as.numeric(sapply(strsplit(config,split = ':'),'[',1))
                  matNum = tumTot-patNum
                  tmp.a = data.frame(Idx=paste(unique(rmIdx.a),collapse = '_'),
                                     Chr=chr,Start=min(a.toRemove$Start),Stop=max(a.toRemove$Stop),
                                     matNum = matNum,tumTot = tumTot,patNum = patNum,segLen=segLen,type='maj',
                                     posID = '0',frac = paste(a.toRemove$frac,collapse = '_'),
                                     tumFrac = matNum/tumTot,tot2min = paste0(tumTot,':',patNum),newIdx = paste0(unique(rmIdx.a)[1],'.1'))
                  newSeg.a=rbind(newSeg.a,tmp.a)
                  rmIdx.a = c()
                }else if(length(rmIdx.a) == 1){
                  message('OOPSS')
                  print(rmIdx.a)
                  stop()
                }
              }
            } 
            if(length(rmIdx.a.all) > 1){
              toRemove = c(toRemove,unique(rmIdx.a.all))
              tmp = newSeg.a    
            }else{
              tmp=NULL
            }
          }else{
            tmp = NULL
          }
          
          
          breakFlag=0
          removeFlag=0
        }
        if(!is.null(tmp)){
          #tmp$tumFrac = tmp$matNum/(tmp$tumTot)  
          #if(sum(is.na(tmp$matNum))>0 & unique(tmp$patNum)==0){
          #  tmp$tumFrac = 1
          #}else{
          #  
          #}
          
          modSegs = rbind(modSegs,tmp)
        }
        if(removeFlag == 1){
          toRemove = c(toRemove,c(chr.data[chr.data$tot2min == config,]$Idx))  
        }
        if(breakFlag == 1){
          break
        }
        
      }  
      
      
      
      
      
    }else if (method == 'allelicRatio'){ 
      for(tumFrac in unique(chr.data$tumFrac)){
        tumFrac.len.perc = sum(chr.data[chr.data$tumFrac == tumFrac,]$segLen)*100/chromLen 
        config = unique(chr.data[chr.data$tumFrac == tumFrac,]$tot2min)
        message(sprintf('\ntumFrac %f for chr %s occupies %f perc, with config %s',tumFrac,chr,tumFrac.len.perc,config))
        
        if(tumFrac.len.perc >= 90){
          message(sprintf('Merging tumFrac %f for chr %s',tumFrac,chr))
          breakFlag=1
          removeFlag=0
          toRemove = c(toRemove,unique(chr.data$Idx))
          if(length(config) > 1){
            message(sprintf('%s - >1 config detected for tumFrac %f',PDID,tumFrac))
            #tumMin = sapply(strsplit(config,split = ':'),'[',2)
            #tumTot = sapply(strsplit(config,split = ':'),'[',1)
            #if(unique(tumMin) == 0){
            #  tmp = data.frame(Idx=0,Chr=chr,Start=1,Stop=chromLen,normTot=2,normMin=1,tumTot = paste(tumTot,collapse=','),tumMin = 0,segLen=chromLen, config=config,idx=0)
            #}else{
            #  stop('WHOOPS')
            #}
          }
          #}else if(length(config) == 1){
          tmp = chr.data[chr.data$tumFrac == tumFrac,][1,]
          tmp$Idx = paste(chr.data[chr.data$tumFrac == tumFrac,]$Idx,collapse = '_')
          tmp$newIdx = paste0(tmp$newIdx,'.', 1)
          tmp$Stop = chromLen
          tmp$Start = 1
          #}
          
        }else if(tumFrac.len.perc <= 10){ # remove if it's <10 % of Chr
          removeFlag=1
          message(sprintf('Removing config %s for chr %s',tumFrac, chr))
          tmp=NULL
          #tmp = data.frame(Idx=0,Chr=chr,Start=1,Stop=chromLen,normTot=2,normMin=1,tumTot = 2,tumMin = 1,segLen=chromLen, config='2:1',idx=0)
          breakFlag=0
        }else{
          # Merge fragments
          a = chr.data[chr.data$tumFrac == tumFrac,]
          rmIdx.a = c()
          rmIdx.a.all = c()
          newSeg.a = tibble()
          
          if(nrow(a) > 1){
            for(f in 2:nrow(a)){
              if((a$Start[f] - a$Stop[f-1]) <= minSegLen){
                message('Merging fragments')
                rmIdx.a = c(rmIdx.a,a$Idx[c(f,(f-1))])
                rmIdx.a.all = c(rmIdx.a.all,a$Idx[c(f,(f-1))])
              }
              if(((a$Start[f] - a$Stop[f-1]) > minSegLen) || f == nrow(a)){
                if(length(rmIdx.a) > 1){
                  a.toRemove = a[a$Idx %in% unique(rmIdx.a),]
                  # Write new merged segment
                  segLen = (max(a.toRemove$Stop) - min(a.toRemove$Start)) + 1
                  patNum = paste(unique(a.toRemove$patNum),collapse = '_')
                  tumTot = paste(unique(a.toRemove$tumTot),collapse = '_')
                  matNum = paste(unique(a.toRemove$matNum),collapse = '_')
                  #patNum = as.numeric(sapply(strsplit(unique(a.toRemove$),split = ':'),'[',2))
                  #tumTot = as.numeric(sapply(strsplit(config,split = ':'),'[',1))
                  #matNum = tumTot-patNum
                  tmp.a = data.frame(Idx=paste(unique(rmIdx.a),collapse = '_'),
                                     Chr=chr,Start=min(a.toRemove$Start),Stop=max(a.toRemove$Stop),
                                     matNum = matNum,tumTot = tumTot,patNum = patNum,segLen=segLen,type='maj',
                                     posID = '0',frac = paste(a.toRemove$frac,collapse = '_'),
                                     tumFrac = tumFrac,tot2min = paste0(tumTot,':',patNum),newIdx = paste0(unique(rmIdx.a)[1],'.1'))
                  newSeg.a=rbind(newSeg.a,tmp.a)
                  rmIdx.a = c()
                  
                }else if(length(rmIdx.a) == 1){
                  message('OOPSS')
                  print(rmIdx.a)
                  stop()
                }
              }
            } 
            if(length(rmIdx.a.all) > 1){
              toRemove = c(toRemove,unique(rmIdx.a.all))
              tmp = newSeg.a    
            }else{
              tmp=NULL
            }
          }else{
            tmp = NULL
          }
          
          
          breakFlag=0
          removeFlag=0
        }
        if(!is.null(tmp)){
          #tmp$patNum = as.numeric(tmp$tumMin)
          #tmp$matNum = as.numeric(tmp$tumTot) - as.numeric(tmp$tumMin)
          #if(sum(is.na(tmp$matNum))>0 & unique(tmp$patNum)==0){
          #  tmp$tumFrac = 1
          #}else{
          #  tmp$tumFrac = tmp$matNum/(tmp$patNum+tmp$matNum)  
          #}
          
          modSegs = rbind(modSegs,tmp)
        }
        if(removeFlag == 1){
          toRemove = c(toRemove,c(chr.data[chr.data$tumFrac == tumFrac,]$Idx))  
        }
        if(breakFlag == 1){
          break
        }
        
      }  
    }
    
  }
  
  
  if((!is.null(toRemove)) & (nrow(modSegs)>0)){
    data = data[!data$Idx %in% toRemove,]
    data = rbind(data, modSegs)
  }
  
  
  data = data[order(data$Chr),]
  
  
  #### Generate final output ####
  new_data = data.frame()
  # Add missing segments
  for(chr in unique(chromInfo$chrom)){
    if(!chr %in% tgtChrs){
      next
    }
    
    chromLen = max(chromInfo[chromInfo$chrom == chr,]$end)*1000
    
    
    if(!chr %in% unique(data$Chr)){
      tmp = data.frame(Idx=0,Chr=chr,Start=1,Stop=chromLen,matNum=1,patNum=1,frac=1,segLen=chromLen,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=paste0('0.',chr))
      tmp$posID = paste0(tmp$Chr,'_',tmp$Start,'_',tmp$Stop)
    }else{
      tmp = data[data$Chr==chr,]
    }
    
    # Add missing segments at the beginning of each Chr
    if(min(tmp$Start) > 1){
      new.tmp = data.frame(Idx=0,Chr=chr,Start=1,Stop=min(tmp$Start)-1,matNum=1,patNum=1,frac=1,segLen=min(tmp$Start)-1,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=paste0('0.',chr))
      new.tmp$posID = paste0(new.tmp$Chr,'_',new.tmp$Start,'_',new.tmp$Stop)
      tmp = rbind(new.tmp,tmp)
    }
    
    # Add missing segments at the end of each Chr
    if(max(tmp$Stop) < chromLen){
      new.tmp = data.frame(Idx=0,Chr=chr,Start=max(tmp$Stop)+1,Stop=chromLen,matNum=1,patNum=1,frac=1,segLen=chromLen,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=paste0('0.1.',chr))
      new.tmp$posID = paste0(new.tmp$Chr,'_',new.tmp$Start,'_',new.tmp$Stop)
      tmp = rbind(new.tmp,tmp)
    }else if(max(tmp$Stop) > chromLen){
      message('Weird STOP...')
      print(chromLen)
      print(tmp[tmp$Stop == max(tmp$Stop),])
      tmp[tmp$Stop == max(tmp$Stop),]$Stop = chromLen
    }
    tmp = tmp[order(tmp$Start),]
    final_tmp=tmp
    for(i in 1:nrow(tmp)){
      if(tmp$Stop[i] == max(tmp$Stop)){
        message(sprintf('%s - Chr %d : DONE!',PDID,chr))
      }else if(tmp$Stop[i] > tmp$Start[i+1]){
        message(sprintf('%s - Chr %d : Overlapping segments: %s',PDID,chr,tmp[i:i+1,]))
      }else if(tmp$Stop[i] < tmp$Start[i+1]-1){
        message(sprintf('%s - Chr %d : Adding segment...',PDID,chr))
        new.tmp = data.frame(0,chr,tmp$Stop[i]+1,tmp$Start[i+1]-1,2,1,2,1,1,1,1,0.5,'2:1',0)
        colnames(new.tmp) = colnames(tmp)
        final_tmp=rbind(new.tmp,final_tmp)
      }else if(tmp$Stop[i] == tmp$Start[i+1]-1){
        message(sprintf('%s - Chr %d : GREAT!',PDID,chr))
      }else{
        message(sprintf('%s: What is happening? chr %s, i = %d',PDID,chr,i))
        print(tmp[i:(i+1),])
        stop()
      }
    }
    
    new_data = rbind(new_data,final_tmp)
  }
  
  
  
  
  ### Add subclone CNAs ###
  subCl.CNA = all[all$type == 'sub',]
  for(chr in unique(subCl.CNA$Chr)){
    maxChromLen = as.vector((chromInfo[(chromInfo$chrom == chr) & (chromInfo$arm == 'q'),]$end)*1000)
    idxs = subCl.CNA[(subCl.CNA$Chr == chr) & (subCl.CNA$Stop > maxChromLen),]$Idx
    if(max(subCl.CNA[subCl.CNA$Chr == chr,]$Stop) > maxChromLen){
      subCl.CNA$Stop[subCl.CNA$Idx %in% idxs] <- rep(maxChromLen,length(idxs))
    }
  }
  new_data = rbind(new_data,subCl.CNA)
  
  ### Process output further, transforming to longer format, add abspos
  if(!'X' %in% tgtChrs){
    new_data$Chr = as.numeric(new_data$Chr)  
  }
  
  #if(PDID == 'PD46693'){
  #  segs[segs$Chr ==4,]$tumFrac = 2/3
  #}else if(PDID == 'PD43255'){
  #  segs[segs$Chr ==2,]$tumFrac = 2/3
  #}else if(PDID == 'PD37104'){
  #  segs[segs$Chr == 9,]$tumFrac = 0.5
  #}
  if(longFormat == T){
    out = pivot_longer(new_data,cols = c('Start','Stop'),names_to = 'posType',values_to = 'pos')
    
    out$idx = c(1:nrow(out))
    # Match max length for each chromosome
    for(chrom in unique(out$Chr)){
      #current_maxChromLen = max(data[data$Chr == chrom,]$pos)
      maxChromLen = as.vector((chromInfo[(chromInfo$chrom == chrom) & (chromInfo$arm == 'q'),]$end)*1000)
      row = out[(out$Chr == chrom) & (out$pos > maxChromLen),]$idx
      if(length(row) > 0){
        print(sprintf('Heey! %d',chrom))
        #out$pos[row] <- rep(maxChromLen,length(row))
      }
    }
    
    # Get absolute genomic position
    out$abspos_kb = out$pos/1000 # if chromosome 1, abspos = pos
    for(r in 1:nrow(out)){
      chrom = out$Chr[r]
      if (chrom > 1){
        out$abspos_kb[r] = out$abspos_kb[r] + (chromInfo[(chromInfo$chrom == (chrom-1)) & (chromInfo$arm == 'q'),]$abspos_kb)
      }
    }
    
  }else{
    out = new_data
  }
  
  # Remove all lines with segLen < 1e6
  out = out[out$segLen >= 1e6,]
  
  # Remove balanced segment
  if(removeBalancedSegs==T){
    out = out[out$tumFrac != 0.5,]
  }
  
  return(out)
}








  
  
  
  
  
  
  
  
  
  
  

  
