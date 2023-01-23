library(stringr)
library(readr)

rd = function(v,o) { 100*(v-o)/o }

## normalize names
normalize_names = function(names) {
    names %>% str_replace("Battarra,etal/","") %>% str_replace("setII/","") %>% str_replace("setIII/","")
}

## get test result file
getFile = function(file) {
  if (is.numeric(file)) {
      test=paste(base,"/test/test",sprintf("%03d",file),"/model.log",sep="")
      if (!file.exists(test)) {
          test=paste(base,"/test/test",sprintf("%03d",file),"/solver.log",sep="")
          if (!file.exists(test)) {
              stop(test," does not exist.\n")
          }
      }
      file=test
  } else if (!file.exists(file)) {
      test=paste(base,"/test/test",sprintf("%03d",as.numeric(file)),"/solver.log",sep="")
      if (!file.exists(test)) {
          stop("Neither ",file," nor ",test," exists.\n")
      }
      file=test
  }
  file
}

## get lbnames for table `t`
getLBnames = function(t) {
    f=file(t)
    tok=str_split(readLines(f,n=1),"\\s+")[[1]]
    close(f)
    p=which(tok=="SUMM")
    if (length(p)==0) {
        p=which(tok=="ILS")
        if (length(p)==0) {
            p = length(tok)+1
        }
    }
    nlb=p-which(tok=="LB")
    if ((nlb-5)==1) {
        "lb3"
    } else {
        k=(nlb-5)/2
        paste(rep(c("lb3","lb3t"),k),as.numeric(gl(k,2))-1,sep=".")
    }
}

## get number of STAT fields
stat_names=c("tag1","inst1","seed1","mem","expansions")

getStat = function(t) {
    f=file(t)
    tok=str_split(readLines(f,n=1)," +")[[1]]
    close(f)
    n=length(tok)
    p=which(tok=="STAT")
    if (length(p)==0) {
        0
    } else {
        n-p+1
    }
}

## read lower bound results (with a variable number of lb3)
readLB = function(name) {
    file=getFile(name)
    lbnames=getLBnames(file)
    cnames=c("tag0","inst0","seed0","lb1","lb2",lbnames)
    d=read.table(file,col.names=cnames) %>% select(-tag0) %>% rename("inst"=inst0)
    d$inst=d$inst %>% normalize_names()
    d %>% as_tibble()
}

## read full results, including lower bounds (with a variable number of lb3)
readResults = function(name,dropLB=F) {
    file=getFile(name)
    ns=getStat(file)
    lbnames=getLBnames(file)
    cnames=c("tag0","inst0","seed0","lb1","lb2",lbnames,"tag","inst","seed","time","ini","heur","htime","ip","ipstat","iptime")
    dropnames=c("tag0","inst0","seed0","tag")
    if (ns>0) {
        cnames=c(cnames,stat_names[1:ns])
        dropnames=c(dropnames,stat_names[1:3])
    }
    d=read.table(file,col.names=cnames) %>% select(-all_of(dropnames))
    d$inst=d$inst %>% normalize_names()
    if (dropLB) {
        d=d%>% select(-all_of(lbnames),-lb1,-lb2)
    }
    d %>% as_tibble()
}

## simple function for test105, models m2i' and m3i
readModels = function(name) {
    file=getFile(name)
    cnames=c("tag0","model0","inst0","lb4","value0","status0","tlb4","tag1","model1","inst","lb5","value1","status1","tlb5")
    d=read.table(file,col.names=cnames) %>% select(-all_of(matches("[01]$"))) %>% relocate(inst,.before=lb4)
    d$inst=d$inst %>% normalize_names()
    d %>% as_tibble()
}
           
## read full results, including lower bounds (with a variable number of lb3); ILS version
readILS = function(name,dropLB=F) {
    file=getFile(name)
    ns=getStat(file)
    lbnames=getLBnames(file)
    cnames=c("tag0","inst0","seed0","lb1","lb2",lbnames,"tag","inst","seed","time","iter","htime","heur","stat")
    dropnames=c("tag0","inst0","seed0","tag")
    if (ns>0) {
        cnames=c(cnames,stat_names[1:ns])
        dropnames=c(dropnames,stat_names[1:3])
    }
    d=read.table(file,col.names=cnames) %>% select(-all_of(dropnames))
    d$inst=d$inst %>% normalize_names()
    if (dropLB) {
        d=d%>% select(-lbnames,-lb1,-lb2)
    }
    d %>% as_tibble()
}

## read full results, including lower bounds (with a variable number of lb3); ILS version
readHeur = function(name,dropLB=F) {
    file=getFile(name)
    cnames=c("tag0","inst0","seed0","lb1","lb2","lb3","tag","inst","seed","sample","multi","uni","tmean","tstd","tmin")
    dropnames=c("tag0","inst0","seed0","tag")
    d=read.table(file,col.names=cnames) %>% select(-all_of(dropnames))
    d$inst=d$inst %>% normalize_names()
    if (dropLB) {
        d=d%>% select(-lb1,-lb2,-lb3)
    }
    d %>% as_tibble()
}

readBounds = function() {
    read.table(paste(base,"dat/Battarra,etal.dat",sep="/"),h=T) %>% as_tibble()
}

readInfo = function(set) {
    df=read_table(paste(base,"/dat/",set,".dat",sep=""),comment="#",col_types=cols())
    if (set=="setI") {
        df = df %>% separate(inst,into=c("type","in"),"(?<=[a-zA-Z])(?=[0-9])",remove=F)
    } else {
        df = df %>% separate(inst,into=c("type","in"),"(?<=p)-",remove=F)  %>% mutate(type=str_to_upper(type))
    }
    df
}
