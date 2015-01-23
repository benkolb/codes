
## Code by B. Kolb for Banking Union project, September 2014

rm(list=ls()) # clear workspace
ILoveMATLABS <- TRUE # switch generation of .mat file on
ILoveHTML <- TRUE    # switch generation of .html output on
DisturbMe <- TRUE    # switch some graphs and checks on
Detrend2 <- FALSE    # use polynomial of order 2 instead of demeaning for series

library(rdatamarket)  # required for download from datamarket.com
library(xts)          # extended time series
library(pdfetch)      # for downloading ECB series
library(pracma)       # for detrending via polynomial
if(ILoveMATLABS) library(R.matlab)     # (optional) generate .mat file
if(ILoveHTML)    library(R2HTML)       # (optional) generate HTML table
if(DisturbMe)    library(tikzDevice)   # (optional) save plot of final series as tikz

options(RCurlOptions = list(cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")))


### Housekeeping
if (Sys.info()[["nodename"]] == "SHAKAKHAN") setwd("C:/Users/Benedikt/Desktop/bu/data")

### 1) Download raw series (which contain series information)

## population
popES_raw <- dmlist("https://datamarket.com/data/set/1a6v/population-and-employment-quarterly-data#!ds=1a6v!qx2=1:qx4=4:qx5=n:7yh1=1")
popGY_raw <- dmlist("https://datamarket.com/data/set/1a6v/population-and-employment-quarterly-data#!ds=1a6v!qx2=1:qx4=4:qx5=d:7yh1=1")

## check for revisions
popG <- as.xts(popGY_raw$Value, as.yearqtr(popGY_raw$Quarter))
if (DisturbMe) plot(diff(popG)/popG)
popS <- as.xts(popES_raw$Value, as.yearqtr(popES_raw$Quarter))
if (DisturbMe) plot(diff(popS)/popS)


## output
# realgdpES_raw <- dmlist("https://datamarket.com/data/set/191o/quarterly-national-accounts-by-6-branches-volumes#!ds=191o!ngj=1:ngl=3:ngn=l:7z0n=9")
# realgdpGY_raw <- dmlist("https://datamarket.com/data/set/191o/quarterly-national-accounts-by-6-branches-volumes#!ds=191o!ngj=1:ngl=3:ngn=c:7z0n=9")
gdpES_raw <- dmlist("https://datamarket.com/data/set/18qt/gdp-and-main-components-current-prices#!ds=18qt!mh2=3.1:mh3=8:mh4=b:mh5=o")
gdpGY_raw <- dmlist("https://datamarket.com/data/set/18qt/gdp-and-main-components-current-prices#!ds=18qt!mh2=1:mh3=8:mh4=b:mh5=f")

## consumption
# realconES_raw <- dmlist("https://datamarket.com/data/set/1ad9/final-consumption-aggregates-volumes#!ds=1ad9!ybd=4:ybe=o:2e9i=2:7z42=7")
# realconGY_raw <- dmlist("https://datamarket.com/data/set/1ad9/final-consumption-aggregates-volumes#!ds=1ad9!ybd=4:ybe=e:2e9i=2:7z42=7")
conES_raw <- dmlist("https://datamarket.com/data/set/18qt/gdp-and-main-components-current-prices#!ds=18qt!mh2=3.1:mh3=8:mh4=a:mh5=o")
conGY_raw <- dmlist("https://datamarket.com/data/set/18qt/gdp-and-main-components-current-prices#!ds=18qt!mh2=1:mh3=8:mh4=a:mh5=f")

## investment
# realinvES_raw <- dmlist("https://datamarket.com/data/set/1aj3/gross-fixed-capital-formation-by-6-asset-types-volumes#!ds=1aj3!z9g=o:1zn6=1:2eoz=2:7z5t=7")
# realinvGY_raw <- dmlist("https://datamarket.com/data/set/1aj3/gross-fixed-capital-formation-by-6-asset-types-volumes#!ds=1aj3!z9g=e:1zn6=1:2eoz=2:7z5t=7")
invES_raw <- dmlist("https://datamarket.com/data/set/18qt/gdp-and-main-components-current-prices#!ds=18qt!mh2=3.1:mh3=8:mh4=g:mh5=o")
invGY_raw <- dmlist("https://datamarket.com/data/set/18qt/gdp-and-main-components-current-prices#!ds=18qt!mh2=1:mh3=8:mh4=g:mh5=f")

## GDP deflator
defES_raw <- dmlist("https://datamarket.com/data/set/18qy/gdp-and-main-components-price-indices#!ds=18qy!mhk=1:mhm=a:mhn=o:7z4p=3")
defGY_raw <- dmlist("https://datamarket.com/data/set/18qy/gdp-and-main-components-price-indices#!ds=18qy!mhk=1:mhm=a:mhn=f:7z4p=3")


## hours worked
hrsES_raw <- dmlist("https://datamarket.com/data/set/191i/quarterly-national-accounts-by-6-branches-employment-data#!ds=191i!nfx=1:nfy=2:nfz=3:ng0=3:ng1=m")
hrsGY_raw <- dmlist("https://datamarket.com/data/set/191i/quarterly-national-accounts-by-6-branches-employment-data#!ds=191i!nfx=1:nfy=2:nfz=3:ng0=3:ng1=d")

## inflation
infES_raw <- dmlist("https://datamarket.com/data/set/1a6e/hicp-2005-100-monthly-data-annual-rate-of-change#!ds=1a6e!qvc=1m:qvd=o")
infGY_raw <- dmlist("https://datamarket.com/data/set/1a6e/hicp-2005-100-monthly-data-annual-rate-of-change#!ds=1a6e!qvc=1m:qvd=e")

infG <- as.xts(infGY_raw$Value, as.yearmon(infGY_raw$Month))
if (DisturbMe) plot(infG)

## leverage of MFIs
# total assets/liabilities of MFIs
assES_ecb <- pdfetch_ECB("BSI.M.ES.N.A.T00.A.1.Z5.0000.Z01.E")
assGY_ecb <- pdfetch_ECB("BSI.M.DE.N.A.T00.A.1.Z5.0000.Z01.E")
assES_raw <- setNames(xts(assES_ecb, as.yearmon(seq(from=as.Date("1997-09-01"),length.out=length(assES_ecb),by="month", format="%Ym%m"))), "assES")
assGY_raw <- setNames(xts(assGY_ecb, as.yearmon(seq(from=as.Date("1997-09-01"),length.out=length(assGY_ecb),by="month", format="%Ym%m"))), "assGY")
# capital and reserves of MFIs
capES_ecb <- pdfetch_ECB("BSI.M.ES.N.A.L60.X.1.Z5.0000.Z01.E")
capGY_ecb <- pdfetch_ECB("BSI.M.DE.N.A.L60.X.1.Z5.0000.Z01.E")
capES_raw <- setNames(xts(capES_ecb, as.yearmon(seq(from=as.Date("1997-09-01"),length.out=length(capES_ecb),by="month", format="%Ym%m"))), "capES")
capGY_raw <- setNames(xts(capGY_ecb, as.yearmon(seq(from=as.Date("1997-09-01"),length.out=length(capGY_ecb),by="month", format="%Ym%m"))), "capGY")
lefES_raw <- assES_raw/capES_raw
lefGY_raw <- assGY_raw/capGY_raw
names(lefES_raw) <- "lefES"
names(lefGY_raw) <- "lefGY"


## loans by MFIs
loaES_ecb <- pdfetch_ECB("BSI.M.ES.N.A.A20.A.1.U2.1000.Z01.E")
loaGY_ecb <- pdfetch_ECB("BSI.M.DE.N.A.A20.A.1.U2.1000.Z01.E")
loaES_raw <- setNames(xts(loaES_ecb, as.yearmon(seq(from=as.Date("1997-09-01"),length.out=length(loaES_ecb),by="month", format="%Ym%m"))), "loaES")
loaGY_raw <- setNames(xts(loaGY_ecb, as.yearmon(seq(from=as.Date("1997-09-01"),length.out=length(loaGY_ecb),by="month", format="%Ym%m"))), "loaGY")


## deposits of MFIs
depES_ecb <- pdfetch_ECB("BSI.M.ES.N.A.L20.A.1.U2.1000.Z01.E")
depGY_ecb <- pdfetch_ECB("BSI.M.DE.N.A.L20.A.1.U2.1000.Z01.E")
depES_raw <- xts(depES_ecb, as.yearmon(seq(from=as.Date("1997-09-01"),length.out=length(depES_ecb),by="month", format="%Ym%m")))
depGY_raw <- xts(depGY_ecb, as.yearmon(seq(from=as.Date("1997-09-01"),length.out=length(depGY_ecb),by="month", format="%Ym%m")))
names(depES_raw) <- "depES"
names(depGY_raw) <- "depGY"

## leverage firms
levES_raw <- pdfetch_ECB("IEAQ.Q.ES.N.V.LE.RFPG.S11.A1.S.2.N.F.Z")
levGY_raw <- pdfetch_ECB("IEAQ.Q.DE.N.V.LE.RFPG.S11.A1.S.2.N.F.Z")
names(levES_raw) <- "levES"
names(levGY_raw) <- "levGY"
# alternatives:
# 1) loans-to-GDP ratios
# levES_raw <- pdfetch_ECB("IEAQ.Q.ES.N.V.LE.RFPG.S11.A1.S.2.N.F.Z")
# levGY_raw <- pdfetch_ECB("IEAQ.Q.DE.N.V.LE.RFPG.S11.A1.S.2.N.F.Z")
# 2) debt-to-GDP ratios
# levES_raw <- pdfetch_ECB("IEAQ.Q.ES.N.V.LE.RF4A.S1M.A1.S.2.N.F.Z")
# levGY_raw <- pdfetch_ECB("IEAQ.Q.DE.N.V.LE.RF4A.S1M.A1.S.2.N.F.Z")


## create ECB policy rate by hand
fixtend <- pdfetch_ECB("FM.B.U2.EUR.4F.KR.MRR_FR.LEV")  # fixed tender, used 99-00 and 08-14
vartend <- pdfetch_ECB("FM.B.U2.EUR.4F.KR.MRR_MBR.LEV") # vari. tender, used 00-08
intcd <- append(fixtend,vartend) # changes in policy rate
changedate <- as.Date(strptime(index(intcd), '%Y-%m-%d')) # extract dates of pol. changes
intm <- na.locf(merge(intcd, zoo(, seq(min(changedate), as.Date(Sys.Date()), "days")))) # fill dates of pol. inaction with the previous policy rate value

spreads <- read.csv2("http://www.banque-france.fr/fileadmin/user_upload/banque_de_france/Economie_et_Statistiques/GM_monthly_EA_Credit_risk_indicators.csv")
sprES_raw <- setNames(xts(spreads$spr_bk_bund_sp,  as.yearmon(spreads$date, format="%Ym%m")),  "sprES")
sprGY_raw <- setNames(xts(spreads$spr_bk_bund_de,  as.yearmon(spreads$date, format="%Ym%m")),  "sprGY")

## 2) Convert raw series to xts objects

countries <- c('ES','GY')
vars <- c('pop','gdp','con','inv','hrs', 'def')
series <- setNames(paste0(vars, rep(countries, each=length(vars)), "_raw"), 
                   paste0(vars, rep(countries, each=length(vars))))
L <- lapply(series, get)
L <- lapply(L, function(ser) xts(ser$Value, as.yearqtr(ser$Quarter))) 
rdat <- do.call(merge, L)


# add averages of monthly inflation and spread series
infES <- apply.quarterly(as.zoo(xts(infES_raw$Value, as.yearmon(infES_raw$Month))), mean)
index(infES) <- as.yearqtr(index(infES))
infGY <- apply.quarterly(as.zoo(xts(infGY_raw$Value, as.yearmon(infGY_raw$Month))), mean)
index(infGY) <- as.yearqtr(index(infGY))
sprES <- apply.quarterly(sprES_raw, mean)
index(sprES) <- as.yearqtr(index(sprES))
sprGY <- apply.quarterly(sprGY_raw, mean)
index(sprGY) <- as.yearqtr(index(sprGY))
loaES <- apply.quarterly(loaES_raw, mean)
index(loaES) <- as.yearqtr(index(loaES))
loaGY <- apply.quarterly(loaGY_raw, mean)
index(loaGY) <- as.yearqtr(index(loaGY))
depES <- apply.quarterly(depES_raw, mean)
index(depES) <- as.yearqtr(index(depES))
depGY <- apply.quarterly(depGY_raw, mean)
index(depGY) <- as.yearqtr(index(depGY))
lefES <- apply.quarterly(lefES_raw, mean)
index(lefES) <- as.yearqtr(index(lefES))
lefGY <- apply.quarterly(lefGY_raw, mean)
index(lefGY) <- as.yearqtr(index(lefGY))
levES <- levES_raw
levGY <- levGY_raw
index(levES) <- as.yearqtr(index(levES))
index(levGY) <- as.yearqtr(index(levGY))

rdat <- merge(rdat,infES,infGY,sprES,sprGY,loaES,loaGY,depES,depGY,levES,levGY,lefES,lefGY)

# add averages of monthly policy rate
int <- apply.quarterly(intm,'mean') # fill with previous value and average for quarters
index(int) <- as.yearqtr(index(int))
names(int) <- "int"
rdat <- merge(rdat,int)


## 3) Data treatment (per-capitalisation, detrending) 
dat <- rdat # keep rdat containing raw data

# Per-capita values for gdp, cons, inv, hrs
dat[,c("gdpES","conES","invES","hrsES","loaES","depES")] <- dat[,c("gdpES","conES","invES","hrsES","loaES","depES")]/rep(dat[, "popES"],6)
dat[,c("gdpGY","conGY","invGY","hrsGY","loaGY","depGY")] <- dat[,c("gdpGY","conGY","invGY","hrsGY","loaGY","depGY")]/rep(dat[, "popGY"],6)

# Deflate nominal series for GDP, cons. and inv. by GDP deflator (see Pfeifer 2014)
dat[,c("gdpES","conES","invES")] <- dat[,c("gdpES","conES","invES")]/rep(dat[, "defES"],3)
dat[,c("gdpGY","conGY","invGY")] <- dat[,c("gdpGY","conGY","invGY")]/rep(dat[, "defGY"],3)

# Divide inflation by 4 for quarterised values
dat[,c("infES","infGY")] <- 1+dat[,c("infES","infGY")]/400

# Create appropriate gross quarterly interest rates and corresponding spread series
dat[,c("int")] <- 1+dat[,c("int")]/400 # approximate geometric mean
dat[,c("sprES","sprGY")] <- dat[,c("sprES","sprGY")]/400

dat <- na.omit(dat[,!(names(dat) %in% c("popES","popGY","defES","defGY"))]) # carry on with dataset w/o NaNs and w/o population & deflator series


## take logs BEFORE detrending (see Pfeifer, 2014) - except for inflation series
dat.detr <- log(dat)

if (Detrend2) {
## detrend via polynomial of order two
polyorder <- 2
trendregs <- seq(1,length(dat[,1]))

for (kol in colnames(dat.detr)) {
  pY <- polyfit(trendregs, dat.detr[,kol],polyorder)
  if (nrow(pY)<ncol(pY)) {pY <- t(pY)}
  Y.trend <- polyval(as.vector(pY), trendregs)
  dat.detr[, kol] <- dat.detr[, kol]-Y.trend
}
} else {
  for (kol in colnames(dat.detr))
  # Calculate the column means
    N_mean<-apply(dat.detr,2,mean)
  # (2) Subtract the columen
    dat.detr<-sweep(dat.detr,2,N_mean)
}

# check results:
if (DisturbMe){
for (kol in 1:length(dat.detr[1,])) {print(mean(dat.detr[,kol]))} # mean of columns (ie series)
plot.zoo(dat.detr,type='l') 
}


## 4) save as .csv, .mat and .html file
# csv
write.zoo(dat.detr, "bu_dat.csv", quote = FALSE, sep = ",")

if (ILoveMATLABS){
writeMat("GK_dat.mat", 
         YFobs=as.numeric(dat.detr$gdpES),    CFobs=as.numeric(dat.detr$conES),   
         IFobs=as.numeric(dat.detr$invES),    LpFobs=as.numeric(dat.detr$hrsES), 
         inflFobs=as.numeric(dat.detr$infES), sprFobs=as.numeric(dat.detr$sprES),
         BFobs=as.numeric(dat.detr$loaES),    DFobs=as.numeric(dat.detr$depES), 
         phiFobs=as.numeric(dat.detr$levES), 
         Yobs=as.numeric(dat.detr$gdpGY),     Cobs=as.numeric(dat.detr$conGY),    
         Iobs=as.numeric(dat.detr$invGY),     Lpobs=as.numeric(dat.detr$hrsGY),  
         inflobs=as.numeric(dat.detr$infGY),  sprobs=as.numeric(dat.detr$sprGY), 
         Bobs=as.numeric(dat.detr$loaGY),     Dobs=as.numeric(dat.detr$depGY),  
         phiobs=as.numeric(dat.detr$levGY),
         Robs=as.numeric(dat.detr$int))

writeMat("GNSS_dat.mat", 
         YFobs=as.numeric(dat.detr$gdpES),    CFobs=as.numeric(dat.detr$conES),   
         IFobs=as.numeric(dat.detr$invES),    LFobs=as.numeric(dat.detr$hrsES), 
         PiFobs=as.numeric(dat.detr$infES), sprFobs=as.numeric(dat.detr$sprES),
         BFobs=as.numeric(dat.detr$loaES),    DFobs=as.numeric(dat.detr$depES), 
         phiFobs=as.numeric(dat.detr$lefES), 
         Yobs=as.numeric(dat.detr$gdpGY),     Cobs=as.numeric(dat.detr$conGY),    
         Iobs=as.numeric(dat.detr$invGY),     Lobs=as.numeric(dat.detr$hrsGY),  
         Piobs=as.numeric(dat.detr$infGY),    sprobs=as.numeric(dat.detr$sprGY), 
         Bobs=as.numeric(dat.detr$loaGY),     Dobs=as.numeric(dat.detr$depGY),  
         phiobs=as.numeric(dat.detr$lefGY),
         Robs=as.numeric(dat.detr$int))
}


# html
if (ILoveHTML) {
htmltable <- HTMLInitFile(getwd(), "htmltable", CSSFile=paste0(getwd(),"/table.css"), useGrid=F)
HTML(as.title("Some Spanish & German macro variables"), file=htmltable)
HTML(as.data.frame(dat.detr), Border=NULL, file=htmltable, digits = 2, nsmall=2)
}

### Plot the series
if (DisturbMe){
  tikz('bu_obs_series.tikz',
       bareBone=T,
       height=3, 
       width=5,
       standAlone=F # I do not want to compile the tikz file
  )
  
nam <- c('Y','C','I','L','infl','spr','B','D','lev')
dat.plot <- dat.detr[, c("gdpES","conES","invES","hrsES","infES","sprES","loaES","depES","levES",
                         "gdpGY","conGY","invGY","hrsGY","infGY","sprGY","loaGY","depGY","levGY")]
plotd <- plot.zoo(dat.plot[, !names(dat.plot) %in% c("int")],
         col = 'blue',
         ylab = nam,
         xlab = NULL,
         nc = 2,
         main = "Transformed data set (left: Spain, right: Germany)", # figure name
)
  
  dev.off()
  
}

if (ILoveHTML){
htmlplot <- HTMLInitFile(getwd(), "plotd", CSSFile=paste0(getwd(),"/table.css"), useGrid=F)
HTML(as.title("Transformed data set"), file=htmlplot)
HTML(plotd, Border=NULL, file=htmlplot)
}
