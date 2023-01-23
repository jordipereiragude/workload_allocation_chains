library(dplyr)
library(ggplot2)
library(tidyr)
library(here)
library(GGally)
library(broom)
base=here()

source(paste(base,"src/scripts/common.R",sep="/"))
source(paste(base,"src/scripts/functions.R",sep="/"))

options(xtable.math.style.negative=F) ## turn off, since NA.string "---" gets set in math mode

setI.opt = readResults(30) %>% select(inst,ip) %>% rename("opt"=ip)
setII.opt = readResults(34) %>% select(inst,ip) %>% rename("opt"=ip)
setIII.opt = readResults(104) %>% select(inst,ip) %>% rename("opt"=ip)

epsrd = function(v,o) {
    reldev = rd(v,o)
    ifelse(abs(reldev)<1e-4,0,reldev)
}

######################################################################
## Table: Overview of instances
######################################################################
df=readInfo("setI") %>% rbind(readInfo("setII")) %>% rbind(readInfo("setIII"))
ef=df %>% group_by(m) %>% summarize(N=length(m),n=mean(n),B=mean(B),s=mean(sumStd))
ot=cbind(ef[1:7,],ef[8:14,],ef[15:21,])

ot.head=c(
    "\\toprule",
    "\\multicolumn{5}{c}{Set I} & \\multicolumn{5}{c}{Set II} & \\multicolumn{5}{c}{Set III}\\\\",
    "\\cmidrule(lr){1-5}\\cmidrule(lr){6-10}\\cmidrule(lr){11-15}"
)
ot.caption="Overview of instances. For each group with the same number of workers $m$, we report number of instances $N$, average number of tasks $n$, average batch size $B$, and average sum of operation times $\\overline T$."
ot.file=paste(base,"doc/tables/instances.tex",sep="/")
ot.label="tab:instances"
colnames(ot)=rep(c("$m$","$N$","$n$","$B$","$\\overline T$"),3)
print(
    xtable(ot,digits=c(0,rep(c(0,0,1,1,1),3)),caption=ot.caption,label=ot.label),
    file=ot.file,
    hline.after=c(0,nrow(ot)),
    size="\\small\n\\medskip\\setlength{\\tabcolsep}{0.8ex}",
    sanitize.colnames=function(x) { return(x) },
    add.to.row=list(pos=list(-1,-1,-1),command=ot.head),
    )

######################################################################
## Tables: comparison of lower bounds
######################################################################

lb3_value="lb3\\.[1-9][0-9]*"
lb3_time="lb3t\\.[1-9][0-9]*"
lbb0=readLB(60) %>% left_join(readInfo("setI")) %>% left_join(readResults(30) %>% select(-starts_with("lb")),by="inst") %>% select(-starts_with("seed"),-sumStd,-B) %>% mutate(across(lb1|lb2|matches(lb3_value),~-epsrd(.,ip))) 
lbb = lbb0 %>% group_by(m) %>% summarize(N=length(m),across(lb1|lb2|matches(lb3_value),list(mv=~mean(.),ct=~100*sum(abs(.)<1e-4)/N)),across(matches(lb3_time),~mean(.,na.rm=T))) %>% relocate(all_of(ends_with("ct")),.after=last_col()) %>% select(-N)
   
### Part 1: bounds
lbb1 = lbb %>% select(m:lb3.12_mv) %>% bind_rows(lbb %>% summarize(across(lb1_mv:lb3.12_mv,mean)))
lbb1.head=c(
    "\\toprule",
    " & & & \\multicolumn{12}{c}{lb3}\\\\",
    "\\cmidrule(lr){4-15}"
)
lbb1.caption="Relative deviations (in percent) from the optimal values of lower bounds lb1, lb2, and lb3 on instance set I."
lbb1.file=paste(base,"doc/tables/lb-deviations.tex",sep="/")
lbb1.label="tab:lb:deviations"
colnames(lbb1)=c("$m$","lb1","lb2",1:12)
print(
    xtable(lbb1,digits=c(0,0,rep(2,14)),caption=lbb1.caption,label=lbb1.label),
    file=lbb1.file,
    hline.after=c(0,nrow(lbb1)-1,nrow(lbb1)),
    size="\\small\n\\medskip\\setlength{\\tabcolsep}{0.75ex}",
    sanitize.colnames=function(x) { return(x) },
    add.to.row=list(pos=list(-1,-1,-1),command=lbb1.head),
    NA.string="Avg."
    )

### Part 2: times
lbb2 = lbb %>% select(m,lb3t.3:lb3t.12,lb3.3_ct:lb3.12_ct)
lbb2 = lbb2 %>% bind_rows(lbb2 %>% summarize(across(lb3t.3:lb3.12_ct,function(x) { mean(x,na.rm=T) })))
lbb2.head=c(
    "\\toprule",
    " & \\multicolumn{10}{c}{Time (s)} & \\multicolumn{10}{c}{Percentage optimal}\\\\",
    "\\cmidrule(lr){2-11}\\cmidrule(lr){12-21}"
)
lbb2.caption="Percentage of optimally solved instances and computation for lb3 controlling $3$ to $12$ workers on instance set I."
lbb2.file=paste(base,"doc/tables/lb-performance.tex",sep="/")
lbb2.label="tab:lb:performance"
colnames(lbb2)=c("$m$",3:12,3:12)
print(
    xtable(lbb2,digits=c(0,0,rep(1,10),rep(0,10)),caption=lbb2.caption,label=lbb2.label),
    file=lbb2.file,
    hline.after=c(0,nrow(lbb2)-1,nrow(lbb2)),
    size="\\small\n\\medskip\\setlength{\\tabcolsep}{0.35ex}",
    sanitize.colnames=function(x) { return(x) },
    add.to.row=list(pos=list(-1),command=str_c(lbb2.head,collapse="")),
    NA.string="-"
    )

### Part 3: comparison to Battarra et al. (2020)
corffile= paste(base,"dat/Battarra,etal.dat",sep="/")
lb3v= lbb0 %>% select(c(inst,starts_with("lb3.")))  %>% pivot_longer(starts_with("lb3."), names_prefix="lb3.", names_to="nw",values_to="lb3")  %>% mutate(nw=as.integer(nw))
corf = read_table(corffile,comment="#") %>% left_join(readResults(30) %>% select(inst,ip)) %>% mutate(LB=-rd(LB,ip)) %>% left_join(lb3v)

corf.file=paste(base,"doc/figures/lb-scatterplot.pdf",sep="/")
pdf(corf.file,9,4)
ggplot(data=corf %>% filter(nw==4|nw==8|nw==12),aes(x=LB,y=lb3))+geom_abline(color="red")+geom_point()+facet_wrap(vars(nw),ncol=1)+xlim(0,60)+ylim(0,20)+labs(x="LB (rel. dev. in %)",y="lb3 (rel. dev. in %)")
dev.off()

### Part 4: lower bounds from IP models
lbb4=readModels(105) %>% left_join(readInfo("setI")) %>% left_join(readResults(30) %>% select(-starts_with("lb")),by="inst") %>% mutate(across(lb4|lb5,~-epsrd(.,ip))) 
lbb4s = lbb4 %>% group_by(m) %>% summarize(N=length(m),across(lb4|lb5,list(mv=~mean(.))),across(tlb4|tlb5,~mean(.,na.rm=T))) %>% relocate(lb5_mv,.after=tlb4) %>% select(-N)

lbb4.head=c(
    "\\toprule",
    " & \\multicolumn{2}{c}{M2I'} & \\multicolumn{2}{c}{M3I}\\\\",
    "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}"
)
lbb4.caption="Relative deviations (in percent) of lower bounds from the linear relaxations of models M2I' and M3I from the optimum and computation times on instance set I."
lbb4.file=paste(base,"doc/tables/lb-ipmodels.tex",sep="/")
lbb4.label="tab:lb:ipmodels"

colnames(lbb4s)=c("$m$","r.d.","t(s)","r.d.","t(s)")

print(
    xtable(lbb4s,digits=c(0,0,rep(1,4)),caption=lbb4.caption,label=lbb4.label),
    file=lbb4.file,
    hline.after=c(0,nrow(lbb4s)),
    size="\\small\n\\medskip\\setlength{\\tabcolsep}{6pt}",
    sanitize.colnames=function(x) { return(x) },
    add.to.row=list(pos=list(-1),command=str_c(lbb4.head,collapse="")),
    NA.string="-"
    )

######################################################################
## Tables: comparison of heuristic methods
######################################################################

bh=readHeur(61) %>% left_join(readInfo("setI")) %>% left_join(setI.opt) %>% mutate(across(lb1:lb3,~-epsrd(.,opt))) %>% mutate(gmin=pmin(sample,multi,uni,tmean,tstd,tmin)) %>% relocate(gmin,.after=tmin) %>% mutate(across(sample:gmin,~epsrd(.,opt))) 

bh.m=bh %>% group_by(m) %>% summarize(across(sample:gmin,mean)) %>% bind_rows(bh %>% summarize(across(sample:gmin,mean)))

ils=readILS(63) %>% left_join(readInfo("setI")) %>% left_join(setI.opt) %>% mutate(heur=epsrd(heur,opt))

ils.m=ils %>% group_by(m) %>% summarize(heur=mean(heur),htime=mean(htime)) %>% bind_rows(ils %>% summarize(heur=mean(heur),htime=mean(htime)))

bils=readBounds() %>% left_join(readInfo("setI")) %>% left_join(setI.opt) %>% mutate(heur=epsrd(bils.avg,opt))
bils.m=bils %>% group_by(m) %>% summarize(bils=mean(heur)) %>% bind_rows(bils %>% summarize(bils=mean(heur)))

hc = bind_cols(bh.m,ils.m %>% rename("m0"=m)) %>% bind_cols(bils.m %>% rename("m1"=m)) %>% select(-m0,-m1) %>% relocate(htime,.after=last_col())
colnames(hc)=c("$m$","Smp.","MS","Uni.","Mean","Std.","Min","Basic","r.d.","BILS","t(s)")

hc=hc %>% select(`$m$`,colnames(hc)[sort(as.numeric(hc[8,2:10]),index.return=T)$ix+1],`t(s)`) %>% relocate(`t(s)`,.after=r.d.)

hc.head=c(
    "\\toprule",
    " & \\multicolumn{2}{c}{ILS} & & & \\multicolumn{6}{c}{Basic heuristics}\\\\",
    "\\cmidrule(lr){2-3}\\cmidrule(lr){6-11}"
)
hc.caption="Relative deviations from the optimum for several heuristics on instance set I."
hc.file=paste(base,"doc/tables/heuristics.tex",sep="/")
hc.label="tab:heuristics"
print(
    xtable(hc,digits=c(0,0,rep(1,10)),caption=hc.caption,label=hc.label),
    file=hc.file,
    hline.after=c(0,nrow(hc)-1,nrow(hc)),
    size="\\small\n\\medskip\\setlength{\\tabcolsep}{0.75ex}",
    sanitize.colnames=function(x) { return(x) },
    add.to.row=list(pos=list(-1),command=str_c(hc.head,collapse="")),
    NA.string="Avg."
)

summarize_results = function(t) {
    t %>% group_by(type,m) %>% summarize(value=mean(heur),time=mean(time)) %>% bind_rows(t %>% summarize(value=mean(heur),time=mean(time)))
}

beam1=readResults(67)  %>% left_join(readInfo("setI"))   %>% left_join(setI.opt)   %>% mutate(heur=epsrd(heur,opt))
beam2=readResults(71)  %>% left_join(readInfo("setII"))  %>% left_join(setII.opt)  %>% mutate(heur=epsrd(heur,opt))
beam3=readResults(106) %>% left_join(readInfo("setIII")) %>% left_join(setIII.opt) %>% mutate(heur=epsrd(heur,opt))
beam1.m=beam1 %>% summarize_results()
beam2.m=beam2 %>% summarize_results()
beam3.m=beam3 %>% summarize_results()
cbfs1=readResults(83)  %>% left_join(readInfo("setI"))   %>% left_join(setI.opt)   %>% mutate(heur=epsrd(heur,opt))
cbfs2=readResults(87)  %>% left_join(readInfo("setII"))  %>% left_join(setII.opt)  %>% mutate(heur=epsrd(heur,opt))
cbfs3=readResults(107) %>% left_join(readInfo("setIII")) %>% left_join(setIII.opt) %>% mutate(heur=epsrd(heur,opt))
cbfs1.m=cbfs1 %>% summarize_results()
cbfs2.m=cbfs2 %>% summarize_results()
cbfs3.m=cbfs3 %>% summarize_results()

ah=bind_cols(beam1.m,cbfs1.m,beam2.m,cbfs2.m,beam3.m,cbfs3.m) %>% select(-c(type...5,m...6,type...9,type...13,m...14,type...17,type...21,m...22))

colnames(ah)=c("Type",rep(c("$m$",rep(c("r.d.","t(s)"),2)),3))
ah[nrow(ah),1]="Avg."

ah.head=c(
    "\\toprule",
    " & \\multicolumn{5}{c}{Set I} & \\multicolumn{5}{c}{Set II} & \\multicolumn{5}{c}{Set III}\\\\",
    " \\cmidrule(lr){3-6}\\cmidrule(lr){8-11}\\cmidrule(lr){13-16}",
    " & & \\multicolumn{2}{c}{Beam} & \\multicolumn{2}{c}{tBBR} & & \\multicolumn{2}{c}{Beam} & \\multicolumn{2}{c}{tBBR} & & \\multicolumn{2}{c}{Beam} & \\multicolumn{2}{c}{tBBR}\\\\",
    "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){8-9}\\cmidrule(lr){10-11}\\cmidrule(lr){13-14}\\cmidrule(lr){15-16}"
)
ah.caption="Relative deviations from the optimum for beam search and truncated BBR on instance sets I, II, and III."
ah.file=paste(base,"doc/tables/adv-heuristics.tex",sep="/")
ah.label="tab:adv:heuristics"
print(
    xtable(ah,digits=c(0,0,rep(c(0,rep(c(2,1),2)),3)),caption=ah.caption,label=ah.label),
    file=ah.file,
    hline.after=c(0,nrow(ah)-1,nrow(ah)),
    size="\\small\n\\medskip\\setlength{\\tabcolsep}{0.8ex}",
    sanitize.colnames=function(x) { return(x) },
    add.to.row=list(pos=list(-1),command=str_c(ah.head,collapse="")),
    NA.string="-"
)

######################################################################
## Tables: comparison of exact methods
######################################################################
summarize_exact_tim = function (t) {
    cat((t %>% summarize(opt=100*sum(ipstat=="opt")/length(ipstat)))$opt,"% optimal\n")
    t %>% group_by(type,m) %>% summarize(time=mean(time)) %>% bind_rows(t %>% summarize(time=mean(time)))
}
summarize_exact_tim_mem = function (t) {
    cat((t %>% summarize(opt=100*sum(ipstat=="opt")/length(ipstat)))$opt,"% optimal\n")
    t %>% group_by(type,m) %>% summarize(time=mean(time),mem=mean(mem)/1024/1024) %>% bind_rows(t %>% summarize(time=mean(time),mem=mean(mem)/1024/1024))
}
summarize_exact_opt_tim_mem = function (t) {
    cat((t %>% summarize(opt=100*sum(ipstat=="opt")/length(ipstat)))$opt,"% optimal\n")
    t %>% group_by(type,m) %>% summarize(opt=100*sum(ipstat=="opt")/length(ipstat),time=mean(time),mem=mean(mem)/1024/1024) %>% bind_rows(t %>% summarize(opt=100*sum(ipstat=="opt")/length(ipstat),time=mean(time),mem=mean(mem)/1024/1024))
}

## beam search heuristic, sadpF exact; omitted
beam1=readResults(90) %>% left_join(readInfo("setI"))  %>% left_join(setI.opt)  %>% mutate(ip=epsrd(ip,opt))
beam2=readResults(91) %>% left_join(readInfo("setII")) %>% left_join(setII.opt) %>% mutate(ip=epsrd(ip,opt))
beam3=readResults(95) %>% left_join(readInfo("setIII"))
beam1.m=beam1 %>% summarize_exact_tim()
beam2.m=beam2 %>% summarize_exact_tim_mem()
beam3.m=beam3 %>% summarize_exact_opt_tim_mem()
## cbfs heuristic, sadpF exact
cbfsF1=readResults(88) %>% left_join(readInfo("setI"))  %>% left_join(setI.opt)  %>% mutate(ip=epsrd(ip,opt))
cbfsF2=readResults(89) %>% left_join(readInfo("setII")) %>% left_join(setII.opt) %>% mutate(ip=epsrd(ip,opt))
cbfsF3=readResults(94) %>% left_join(readInfo("setIII"))
cbfsF1.m=cbfsF1 %>% summarize_exact_tim()
cbfsF2.m=cbfsF2 %>% summarize_exact_tim_mem()
cbfsF3.m=cbfsF3 %>% summarize_exact_opt_tim_mem()
## cbfs heuristic, cbfsF exact, 5M sadpF limit
cbfsC1=readResults(99)  %>% left_join(readInfo("setI"))  %>% left_join(setI.opt)  %>% mutate(ip=epsrd(ip,opt))
cbfsC2=readResults(100) %>% left_join(readInfo("setII")) %>% left_join(setII.opt) %>% mutate(ip=epsrd(ip,opt))
cbfsC3=readResults(101) %>% left_join(readInfo("setIII"))
cbfsC1.m=cbfsC1 %>% summarize_exact_tim()
cbfsC2.m=cbfsC2 %>% summarize_exact_tim_mem()
cbfsC3.m=cbfsC3 %>% summarize_exact_opt_tim_mem()
## cbfs heuristic, cbfsF exact, 50M sadpF limit
cbfsD3=readResults(104) %>% left_join(readInfo("setIII"))
cbfsD3.m=cbfsD3 %>% summarize_exact_opt_tim_mem()

ae=bind_cols(cbfsF2.m,cbfsC2.m,cbfsF3.m,cbfsC3.m,cbfsD3.m) %>% rename("type"=type...1) %>% select(-all_of(starts_with("type..."))) %>% select(-c(m...6,m...15,m...20))
ae[nrow(ae),1]="Avg."

colnames(ae)=c("Type",rep(c("$m$",rep(c("t","Mm."),2)),1),"$m$",rep(c("opt","t","Mm."),3))
ae.head=c(
    "\\toprule",
    " & \\multicolumn{5}{c}{Set II} & \\multicolumn{10}{c}{Set III}\\\\",
    " \\cmidrule(lr){3-6}\\cmidrule(lr){7-16} ",
    " & & \\multicolumn{2}{c}{SADP} & \\multicolumn{2}{c}{BBR} & & \\multicolumn{3}{c}{SADP} & \\multicolumn{3}{c}{BBR} & \\multicolumn{3}{c}{BBR (high mem.)}\\\\",
    "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){8-10}\\cmidrule(lr){11-13}\\cmidrule(lr){14-16}"
)
ae.caption="Performance of SADP and BBR on instance II and III."
ae.file=paste(base,"doc/tables/exact.tex",sep="/")
ae.label="tab:exact"
print(
    xtable(ae,digits=c(0,0,rep(c(0,rep(c(1,1),2)),1),0,rep(c(1,1,1),3)),caption=ae.caption,label=ae.label),
    file=ae.file,
    hline.after=c(0,nrow(ae)-1,nrow(ae)),
    size="\\small\n\\medskip\\setlength{\\tabcolsep}{0.5ex}",
    sanitize.colnames=function(x) { return(x) },
    add.to.row=list(pos=list(-1),command=str_c(ae.head,collapse="")),
    NA.string="-"
)
