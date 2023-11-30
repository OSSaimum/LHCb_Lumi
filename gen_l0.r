# the code below is designed for the analysis of mu scan of LHCb.
# however, in principle, the codes can be used for other types of runs as well, e.g., VdM scan.
# the key functions for calculating the mu are:
#   a. the cut() function to determine the cut that characterizes the ee signal.
#   b. the lz() function that calculates the mu using the log zero method.
#   c. the sp_lz() function that calculates the mu using the log zero method for counters that require 'sp'ecial consideration(s).
#   d. the mu_calc() function that calculates the mu for the given data located in a local or virtual directory.


# call the required libraries
library(ggplot2)
library(fst)
library(data.table)
library(bit64)
library(dplyr)

txt.8 <- element_text(size = 8)
txt.16 <- element_text(size = 16)

# define function q() for {easier} plotting
q <- function (d, ...) qplot(data = d, ...)
# define the function to get the time 
T = function(time) as.POSIXct(time, origin = "1970-1-1", tz = "CET")


dir.to.hist.files <- 'C:\\Users\\saimum\\R\\plume\\mu_scan_june\\re'


### characterizing the bunch crossings, i.e., bb, be, eb, ee

# set the directory to the location where the hist files are
setwd("C:/Users/saimum/R/plume/mu_scan_june/re")

jdt <- data.table(file = list.files(dir.to.hist.files, pattern='^run.*gz')
                  )[, {
                        . <- fread(file, col.names =  c('minute','name','bx','x','y'))
                    }, by = .(file)]

# set the directory to the location where the ee files are
# only execute if specific ee bunch crossing info is needed
dir.to.ee.files <- 'C:\\Users\\saimum\\R\\plume\\june_ee'

jeesdt <- data.table(file = list.files(dir.to.ee.files, pattern='^run.*gz')
                  )[, {
                        . <- fread(file, col.names =  c('minute','name','bx','x','y'))
                    }, by = .(file)]


# set the directory to the location where the bb and/or be files are

setwd("C:/Users/saimum/R/plume/june_bb")

bb_dt <- fread('run_268205_bb.gz', col.names = c('minute','name','bx','x','y'))

# obtain the bb bunch crossings
bbs <- unique(bb_dt$bx)

be_dt <- fread('run_268205_be.gz', col.names = c('minute','name','bx','x','y'))

# obtain the be bunch crossings
bes <- unique(be_dt$bx)

# obtain the eb bunch crossings; note that usually all ees are stored in one bxID labeled as 10000
ebs <- list()
ebs <- setdiff(unique(jdt$bx), c(bbs,bes,10000))

# obtain the ee bunch crossings (if needed)
ees <- list()
ees <- setdiff(unique(jeesdt$bx), c(bbs, bes, ebs))



### cut function

setwd("C:/Users/saimum/R/plume/mu_scan_june/re")

# set of counters to be used (add counters as may require)
counters = list("ECalET", "ECalETInnerBottom",    "ECalETInnerTop",      
  "ECalETMiddleBottom",   "ECalETMiddleTop",      "ECalETOuterBottom",
  "ECalETOuterTop",       "ECalEtot",             "RawECalEInnerBottom",       
  "RawECalEInnerTop",     "RawECalEMiddleBottom", "RawECalEMiddleTop",
  "RawECalEOuterBottom",  "RawECalEOuterTop",     "SciFiT1M123",
  "SciFiT1M4",            "SciFiT2M123",          "SciFiT2M4",
  "SciFiT3M123",          "SciFiT3M45",           "FiducialVeloVertices",      
  "MuonHitsM2R1",         "MuonHitsM2R2",         "MuonHitsM2R3",
  "MuonHitsM2R4",         "MuonHitsM3R1",         "MuonHitsM3R2",
  "MuonHitsM3R3",         "MuonHitsM3R4",         "MuonHitsM4R1",        
  "MuonHitsM4R2",         "MuonHitsM4R3",         "MuonHitsM4R4",
  "MuonTracks",           "VeloFiducialTracks",   "VeloTracks",
  "VeloTracksEtaBin0",    "VeloTracksEtaBin1",    "VeloTracksEtaBin2",
  "VeloTracksEtaBin3",    "VeloTracksEtaBin4",    "VeloTracksEtaBin5",
  "VeloTracksEtaBin6",    "VeloTracksEtaBin7",    "VeloVertices")

# list of different types of counters that may be required
ECals = list('ECalETOuterTop', 'ECalETOuterBottom', 'ECalETMiddleTop', 'ECalETMiddleBottom', 'ECalETInnerTop', 'ECalETInnerBottom')

SciFis = list('SciFiT1M123','SciFiT1M4','SciFiT2M123','SciFiT2M4','SciFiT3M123','SciFiT3M45')

Raws = list('RawECalEOuterTop', 'RawECalEOuterBottom', 'RawECalEMiddleTop', 'RawECalEMiddleBottom', 'RawECalEInnerTop', 'RawECalEInnerBottom')

ETots = list('ECalET','ECalEtot')

# copy paste the counters to the names variable to indicate for which counters the cut is desired to be defined.
names = list('ECalET','ECalEtot','ECalETOuterTop', 'ECalETOuterBottom', 'ECalETMiddleTop', 'ECalETMiddleBottom', 'ECalETInnerTop', 'ECalETInnerBottom','RawECalEOuterTop', 'RawECalEOuterBottom', 'RawECalEMiddleTop', 'RawECalEMiddleBottom', 'RawECalEInnerTop', 'RawECalEInnerBottom')

# copy paste the counters to the names variable to indicate for which counters 'sp'ecial attention is needed.
sp_counts = list('SciFiT1M123','SciFiT1M4','SciFiT2M123','SciFiT2M4','SciFiT3M123','SciFiT3M45')


# create an initial data table w/ the desired counters for the cut function
dt1 <- data.table(file = list.files('.',pattern='^run.*gz$')
                 )[, mu :=  as.numeric(sub('run_[0-9]+_mu_', '', sub('.gz$', '', file)))
                   ][, {
                       . <- fread(file, col.names =  c('minute','name','bx','x','y'))[name %in% names]
                       .[,bx.type := ifelse(bx %in% bbs, 'bb', ifelse(bx %in% bes, 'be', ifelse(bx %in% ebs, 'eb', 'ee')))]
                   }, by = .(file, mu)]

# keep only ees (can be omitted; see lines 130-133)
in_dt1 <- subset(in_dt1, bx.type=='ee')

# create an initial data table w/o the desired counters for the cut function
dt2 <- data.table(file = list.files('.',pattern='^run.*gz$')
                 )[, mu :=  as.numeric(sub('run_[0-9]+_mu_', '', sub('.gz$', '', file)))
                   ][, {
                       . <- fread(file, col.names =  c('minute','name','bx','x','y'))[!name %in% names&!name %in% sp_counts]
                       .[,bx.type := ifelse(bx %in% bbs, 'bb', ifelse(bx %in% bes, 'be', ifelse(bx %in% ebs, 'eb', 'ee')))]
                   }, by = .(file, mu)]

# keep only ees
in_dt2 <- subset(in_dt2, bx.type=='ee')

# function to get the cut for the desired counters
# input: a. bin0: initial binning set for ECals in the root file
#        b. bin0: initial binning set for Raws in the root file
#        c. per: set how much of the cumulative sum of the Gaussian to be considered
#        d. jmp: set how far from the distance between the mean and the 'per'th percentile to be added
#        e. names: counters to be considered for cut
cut <- function(bin0,rawbin0,per,jmp,names) {
    # select only ees; note: may seem redundant since it is performed in line 110.
    # however, after running the first time, the in_dt1 data table will be modified.
    # as such, any subsequent execution of the cut function is likely to cause an error.
    # note: line 110 may be removed but kept since it may be needed for analysis.
    in_dt1 <- subset(dt1, bx.type=='ee')
    n1 <- in_dt1[, {
                   # for each counter, obtain the 'combined' histogram
                   . <- in_dt1[,.(y=sum(y)),by=.(name,x)]
                   # for each counter, obtain the cumulative sum
                   .[,ycs:=cumsum(y),by=.(name)]
                   # obtain the median of the histogram
                   .[,peak:=max(y),by=.(name)]
                   # get the x value of the median
                   .[,ipeak:=ifelse(y == peak, x, 0),by=.(name)]
                   .[,peakx:=ifelse(ipeak == 0, max(ipeak), ipeak), by=.(name)]
                   # get the y value up to which threshold the cut is desired
                   .[,ycut:=ceiling(ycs[length(ycs)]*per),by=.(name)]
                   # obtain the first value above that cut
                   .[,cscut:=ycs[ycs>ycut][1],by=.(name)]
                   # get the x value for that cut; note: y value is taken and then the x which may be redundant;
                   # but it's done for easier debugging
                   .[,icut:=ifelse(ycs == cscut, x, 0),by=.(name)]
                   .[,icut:=ifelse(icut == max(icut), icut, -1),by=.(name)]
                   # to make sure the ee signal doesn't affect, we move the 'jmp'*difference between the cut and the peak distance away
                   .[icut!=-1,.(xcut=icut+jmp*max(abs(icut-peakx),1)),by=.(name)]
               },]
    # multiply the initial binning (done when converting from the root files)
    n1 <- n1[,.(cut=ifelse(name %in% ETots, xcut*bin0, ifelse(name %in% ECals, xcut*bin0, ifelse(name %in% Raws, xcut*rawbin0, xcut)))),by=name]
    # for the non-desired counters, set the cut to simply zero
    n2 <- in_dt2[,.(cut=0),by=name]
    # add the two data table objects and return
    all <- dplyr::bind_rows(n1, n2)
    data.table(all)
}

# get the cut specified by the parameters
get_cut <- cut(20,100,.95,5,names)

# from my analysis so far, 95%, and the 5 'jmp' looks ok.
# it is advised to change 'jmp' instead of 'per'.
# harsher cut may be useful depending on the data but please do check if the cut is 'too' much.
# the following function may be executed to check the ratio of the ee signal that falls ABOVE the cut:
ratio <- function(bin0, rawbin0) {
    # create an empty data table with two rows w/ column names 'name' and 'ratio'
    rt <- data.table()
    rt[,name := character(0)]
    rt[,ratio:= numeric(0)]
    # loop over all the counters counters
    for (i in 1:length(counters)){
        # obtain the number of events ABOVE the cut
        right <- ifelse(grepl('^ECal',counters[i]),sum(in_dt1[name==counters[i]&bx.type=='ee'&x*bin0>get_cut[name==counters[i]]$cut]$y),
            ifelse(counters[i]%in%Raws,sum(in_dt1[name==names[i]&bx.type=='ee'&x*rawbin0>get_cut[name==counters[i]]$cut]$y),
            sum(in_dt2[name==counters[i]&bx.type=='ee'&x>get_cut[name==counters[i]]$cut]$y)))
        # obtain the total number of events
        all <- ifelse(counters[i]%in%names,sum(in_dt1[name==counters[i]&bx.type=='ee']$y),sum(in_dt2[name==counters[i]&bx.type=='ee']$y))
        # obtain the ratio
        rto <- right/all
        # add a row to the data table rt with the obtained ratio corresponding to the ith counter
        rt <- rbind(rt,c(counters[i],rto))
    }
    data.table(rt)
}
rto <- ratio(20,100)


# note: one can repeat the ratio check for 'sp'ecial cuts specified below as well.
# since the 'special' counters were not considered at the end of my analysis, i haven't repeated it.


# obtain the 'sp'ecial cuts for the 'sp'ecified counters:
# the same algorithm as before except the cut is obtained per file/mu for each counter instead of just per counter
sp_cut <- function(bin0,rawbin0,per,jmp,sp_count) {
    n1 <- data.table(file = list.files('.',pattern='^run.*gz$')
                 )[, mu :=  as.numeric(sub('run_[0-9]+_mu_', '', sub('.gz$', '', file)))
                   ][, {
                       . <- fread(file, col.names =  c('minute','name','bx','x','y'))[name %in% sp_count]
                       .[,bx.type := ifelse(bx %in% bbs, 'bb', ifelse(bx %in% bes, 'be', ifelse(bx %in% ebs, 'eb', 'ee')))]
                       . <- subset(., bx.type=='ee')
                       # next line may be omitted w/o any issues since this code is run per file
                       .[,.(y=sum(y)),by=.(name,x)]
                       .[,ycs:=cumsum(y),by=.(name,bx.type)]
                       .[,peak:=max(y),by=.(name,bx.type)]
                       .[,ipeak:=ifelse(y == peak, x, 0),by=.(name,bx.type)]
                       .[,peakx:=ifelse(ipeak == 0, max(ipeak), ipeak), by=.(name,bx.type)]
                       .[,ycut:=round(ycs[length(ycs)]*per),by=.(name,bx.type)]
                       .[,cscut:=ycs[ycs>ycut][1],by=.(name)]
                       .[,icut:=ifelse(ycs == cscut, x, 0),by=.(name,bx.type)]
                       .[,icut:=ifelse(icut == max(icut), icut, -1),by=.(name,bx.type)]
                       .[icut!=-1,.(xcut=icut+jmp*max(abs(icut-peakx),1)),by=.(name,bx.type)]
                   }, by =.(file,mu)]
}

# get the 'sp'ecial cut specified by the parameters
get_sp_cut <- sp_cut(1,1,.95,5,sp_counts)

## note: ideally, 'per' and 'jmp' should be the same for cut and sp_cut



### Define a function to calculate the mu using the log zero method

lz <- function(x, y, nm) {
    # obtain the 'sp'ecified cut for the 'sp'ecific counter
    sp_cut <- get_cut[name==nm]$cut
    # obtain the number of empty events
    emp_ev <- sum(y[x<=sp_cut])
    # obtain the number of all events
    all_ev <- sum(y)
    # calculate the mu using the log zero method
    lz_mu <- if(emp_ev==0){log(sum(y))
        } else (-log(emp_ev/all_ev)-0.5*(1/emp_ev-1/all_ev))
    # lz_err <- if(emp_ev==0){1
    #     } else if (emp_ev==all_ev) {
    #         1/all_ev
    #     } else sqrt(1/emp_ev-1/all_ev)
    list(lz_mu)
}


sp_lz <- function(x, y, nm, sp_mu) {
    # obtain the 'sp'ecified cut for the 'sp'ecific counter AND mu
    sp_cut <- get_sp_cut[name==nm&sp_mu==mu]$cut
    # obtain the number of empty events
    emp_ev <- sum(y[x<=sp_cut])
    # obtain the number of all events
    all_ev <- sum(y)
    # calculate the mu using the log zero method
    lz_mu <- if(emp_ev==0){log(sum(y))
        } else (-log(emp_ev/all_ev)-0.5*(1/emp_ev-1/all_ev))
    # lz_err <- if(emp_ev==0){1
    #     } else if (emp_ev==all_ev) {
    #         1/all_ev
    #     } else sqrt(1/emp_ev-1/all_ev)
    list(lz_mu)
}



# function to calculate the mu for the desired data set
# input: a. dir: the directory where the hist files are
#        b. bin0: initial binning set for ECals in the root file
#        c. rawbin0: initial binning set for Raws in the root file

mu_calc <- function(dir,bin0,rawbin0){
    setwd(dir)
    dt <- data.table(file = list.files('.',pattern='^run.*gz$')
                 )[, mu :=  as.numeric(sub('run_[0-9]+_mu_', '', sub('.gz$', '', file)))
                    ][, {
                    print(file)
                       . <- fread(file, col.names =  c('minute','name','bx','x','y'),colClasses=c('integer','character','integer','integer','integer'))
                       # for the concerned counters, multiply x with appropriate initial binnings
                       . <- .[grepl('^ECal',name), x:=x*bin0]
                       . <- .[grepl('^Raw', name), x:=x*rawbin0]
                       # calculate the mu per minute, counter, and bunch crossings based on whether the counter needs 'sp'ecial attention
                       . <- .[,.(my_mu = ifelse(!name %in% sp_counts,as.numeric(lz(x,y,name)),as.numeric(sp_lz(x,y,name,mu))),xa = sum(x*y)/sum(y),n=sum(y)), by=.(minute,name,bx)]
                       # characterize the bunch crossing types
                       .[,bx.type := ifelse(bx %in% bbs, 'bb', ifelse(bx %in% bes, 'be', ifelse(bx %in% ebs, 'eb', 'ee')))]
                       # take the average of the obtained mu values per counter, per minute, AND per bunch crossing type
                       .[, .(lz_mu=sum(my_mu*n)/sum(n),avg_mu = sum(xa*n)/sum(n), n = sum(n)), by = .(minute,name,bx.type)]
                   }, by = .(file,mu)]
}

# calculate the mu values for the specified data set
jmu <- mu_calc("C:/Users/saimum/R/plume/mu_scan_june/re",20L,100L)

# store the data locally
write.fst(jmu, "C:/Users/saimum/R/plume/mu_scan_june/re/jmu_113023.fst")

# read the locally stored data
jmu <- read.fst('C:/Users/saimum/R/plume/mu_scan_june/re/jmu_113023.fst',as.data.table=TRUE)

# read the TIMBER data supplied centrally by LHC
l <- fread('C:\\Users\\saimum\\R\\plume\\mu_scan_june\\TIMBER_data.csv',col.names=c('t','lumi'))

# add the time
l[, t := t/1e6][,time := T(t)][, t.l := t]

# concatenate with the given data based on the time information; note that it may be needed to modify depending on the data and the run.
# one is advised to check to what extent the time indicated in the TIMBER data from LHC aligns with the data to be analyzed.
# for example, if the TIMBER data has a 10 seconds delay compared to the data to be analyzed, one should run the following line of code instead:
# cmp <- l[copy(jmu)[, t.dt :=  minute+10], on=.(t=t.dt),roll='nearest'][abs(t-minute-10)<2]
# however, it should not cause any major issues in the analysis.
cmp <- l[copy(jmu)[, t.dt :=  minute], on=.(t=t.dt),roll='nearest'][abs(t-minute)<2]


# arrange the dataset based on the bunch crossing type for each counter column-wise for the subsequent operations
jmu1 <- dcast(cmp, file+mu+lumi+minute+name~bx.type, value.var=c('lz_mu','avg_mu','n'),sep='.')

# fix the bias from be, eb, and ee to get the proper mu value using log zero method
jmu1[, lz := lz_mu.bb - lz_mu.be - lz_mu.eb + lz_mu.ee]

# calclatte the total number of events (optional, not likely to be needed)
# jmu1[, nt:=n.bb+n.be+n.eb+n.ee]

# fix the bias from be, eb, and ee to get the proper mu value using the average method (optional)
jmu1[, avg := avg_mu.bb - avg_mu.be - avg_mu.eb + avg_mu.ee]

# arrange the dataset based on the counters column-wise
jmu2 <- dcast(jmu1, file+mu+minute+lumi+n.bb~name, value.var=c('lz'))

# add the plume mu values to the dataset (optional)
jmu2[,mu.Plume:=lumi/length(bbs)*63.4/11.245]
setnames(jmu2, 'mu', 'mu.runDB')

# disregard the minutes when the number of bb events is less than 30% of the median
jmu3 <- jmu2[n.bb > 0.3 * median(n.bb)]

## optional
# if the mu values from the average method were also obtained, for easier calculation, take the avg ones out (or log zero ones if the avg method is desired)
jmu3 <- subset(jmu1, select=-avg)
jmu4 <- dcast(jmu3, file+mu+minute+lumi~name, value.var='lz')
jmu4[,mu.Plume:=lumi/length(bbs)*63.4/11.245]
setnames(jmu4, 'mu', 'mu.runDB')



### plotting

# not really needed; just that i prefer to save files by date and the version on each date for ease.
date = 113023
ver = 1

# rearrange the data table row-wise for plotting (the last argument of id.vars is the counter w.r.t which we desire to do the plots, in this case, 'VeloVertices')
pt4 <- melt(jmu3, id.vars = c('file','mu','minute','lumi','n.bb','VeloVertices'),variable.name='name',value.name='x')
print(q(pt4, VeloVertices, x, size = I(0.4)) + facet_wrap(~name, scale='free')+ geom_smooth(method=lm, formula=y~x+0, linewidth = 0.5) +
    labs(x = 'VeloVertices', y = 'Average per minute'))
ggsave(paste0(date,'_',ver,'_mu_VV.pdf'),width = 297, height = 210, units = "mm")

# obtain a data table for plotting the relative residuals for each counter
rel_res <- pt4[name !=  'mu', {
        # get the linear model  w.r.t. to the chosen counter (in this case, VeloVertices)
        # set the intersection to zero
        fit <- lm(x ~ VeloVertices+0)
        # get the values obtained from the linear model
        vl <- fit$fitted.values
        # calculate the relative residual
        .(rel_res = (x-vl)/vl, x = x, VeloVertices)
    }, by = name]

# plot the relative residuals
print(q(jmu_res, VeloVertices, rel_res, size = I(0.6)) + facet_wrap(~name, scale='free')+ geom_smooth(method=lm, formula=y~x+0) +
    labs(y = paste0('Residual of linear fit (Counter = a*VeloVertices + 0) / fit value')))
ggsave(paste0(date,'_',ver,'_rel_res_VV.pdf'),width = 297, height = 210, units = "mm")


## re-do the process for other desired counters
# example done for another counter - 'VeloTracks'

pt5 <- melt(jmu3, id.vars = c('file','mu','minute','lumi','n.bb','VeloTracks'),variable.name='name',value.name='x')
print(q(pt5, VeloTracks, x, size = I(0.4)) + facet_wrap(~name, scale='free')+ geom_smooth(method=lm, formula=y~x+0, linewidth = 0.5) +
    labs(x = 'VeloTracks', y = 'Average per minute'))
ggsave(paste0(date,'_',ver,'_mu_VT.pdf'),width = 297, height = 210, units = "mm")

rel_res_vt <- pt5[name !=  'mu', {
        fit <- lm(x ~ VeloTracks+0)
        vl <- fit$fitted.values
        .(rel_res = (x-vl)/vl, x = x, VeloTracks)
    }, by = name]

print(q(rel_res_vt, VeloTracks, rel_res, size = I(0.6)) + facet_wrap(~name, scale='free')+ geom_smooth(method=lm, formula=y~x+0) +
    labs(y = paste0('Residual of linear fit (Counter = a*VeloTracks + 0) / fit value')))
ggsave(paste0(date,'_',ver,'_rel_res_VT.pdf'),width = 297, height = 210, units = "mm")



