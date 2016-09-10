
library(ape)
#install.packages("rcolgem", repos="http://R-Forge.R-project.org")
library(rcolgem)

treetxt <- "((((((TX_22-Jul-2013_1_KJ645639.1_2013.56:0.02294579,KS_15-Jul-2013_1_KJ645637.1_2013.54:0.00294579):0.41363576,IA_01-May-2013_86_KF468753.1_2013.34:0.21658155):0.03618120,(KS_01-Sep-2013_1_KJ645643.1_2013.67:0.54658155,((OK_26-Aug-2013_1_KJ645642.1_2013.65:0.00000002,TX_12-Nov-2013_1_KJ645697.1_2013.87:0.22000002):0.48638763,OK_07-Sep-2013_1_KJ645644.1_2013.69:0.52638766):0.04019389):0.03618120):0.00000000,((CO_22-Jul-2013_1_KJ645638.1_2013.56:0.47273389,(((OK_24-Jul-2013_1_KJ645640.1_2013.57:0.34563379,KS_27-Jan-2014_1_KJ645701.1_2014.07:0.84563379):0.00001443,(KS_23-Oct-2013_1_KJ645650.1_2013.81:0.17914293,CO_23-Oct-2013_1_KJ645651.1_2013.81:0.17914293):0.40650528):0.13707126,(NC_05-Oct-2013_1_KJ645646.1_2013.76:0.37603835,((MN_14-Oct-2013_1_KJ645647.1_2013.79:0.07622385,(MN_10-Dec-2013_1_KJ645687.1_2013.94:0.20422629,(((MN_22-Jan-2014_1_KM077139.1_2014.06:0.15001443,MN_28-Nov-2013_1_KJ645673.1_2013.91:0.00001443):0.17419744,(MN_04-Dec-2013_1_KJ645678.1_2013.93:0.07669231,(MN_20-Nov-2013_1_KJ645661.1_2013.89:0.00000000,NC_23-Nov-2013_1_KJ645662.1_2013.90:0.01000000):0.03669231):0.11751955):0.00001443,(((MN_02-Dec-2013_1_KJ645677.1_2013.92:0.07698571,IA_24-Nov-2013_1_KJ645666.1_2013.90:0.05698571):0.10268761,MN_30-Nov-2013_1_KJ645674.1_2013.92:0.17967331):0.00455297,(MN_07-Dec-2013_1_KJ645679.1_2013.94:0.05126423,MN_10-Dec-2013_1_KJ645686.1_2013.94:0.05126423):0.15296206):0.00000000):0.00000000):0.02199756):0.03822558,IA_11-Dec-2013_1_KJ645688.1_2013.95:0.27444943):0.29158892):0.29668111):0.00001443):0.00001443,CO_10-May-2013_1_KF272920.1_2013.36:0.27274832):0.00001443):0.09769420,Mexico_03-Nov-2013_1_KJ645708.1_2013.84:0.85045694):0.79111263,((IN_21-Aug-2013_1_KJ645641.1_2013.64:0.99120571,((((MO_11-Dec-2013_1_KJ645685.1_2013.95:0.03559679,MO_11-Dec-2013_1_KJ645684.1_2013.95:0.03559679):0.00000000,MO_15-Dec-2013_1_KJ645693.1_2013.96:0.04559679):0.11494596,MN_02-Feb-2014_1_KJ645703.1_2014.09:0.29054275):0.36200548,(((MN_25-Nov-2013_1_KJ645668.1_2013.90:0.05721679,((((MN_02-Dec-2013_1_KJ645676.1_2013.92:0.01001387,MN_28-Nov-2013_1_KJ645672.1_2013.91:0.00001387):0.00998614,MN_23-Nov-2013_1_KJ645663.1_2013.90:0.00000000):0.02158519,MN_04-Dec-2013_1_KJ645682.1_2013.93:0.05158520):0.03561716,MN_04-Dec-2013_1_KJ645681.1_2013.93:0.08720236):0.00001443):0.02731198,MN_18-Nov-2013_1_KJ645658.1_2013.88:0.06452876):0.26187848,MN_18-Nov-2013_1_KJ645705.1_2013.88:0.32640724):0.11614099):0.78865748):0.00001443,((TN_04-Nov-2013_1_KJ645654.1_2013.84:1.08388321,(((OH_14-Jun-2013_1_KJ584361.1_2013.46:0.50983262,(((OH_27-Nov-2013_1_KJ645670.1_2013.91:0.03998875,OH_24-Nov-2013_1_KJ645665.1_2013.90:0.02998875):0.00001443,OH_12-Nov-2013_1_KJ645657.1_2013.87:0.00000318):0.30606791,(WI_03-Nov-2013_1_KJ645653.1_2013.84:0.20202941,IA_15-Dec-2013_1_KJ645694.1_2013.96:0.32202941):0.07404167):0.61376154):0.13658648,(WI_27-Nov-2013_1_KJ645669.1_2013.91:0.02903195,(((MO_15-Dec-2013_1_KJ645692.1_2013.96:0.03437055,IL_11-Dec-2013_1_KJ645690.1_2013.95:0.02437055):0.02577397,MN_30-Nov-2013_1_KJ645671.1_2013.92:0.02014452):0.00000000,IL_03-Dec-2013_1_KJ645675.1_2013.92:0.02014452):0.01888743):1.06738716):0.02491402,(IA_29-Apr-2013_1_KF452322.1_2013.33:0.54131871,((IN_07-May-2013_1_KF452323.1_2013.35:0.00000000,IN_16-May-2013_1_KF650370.1_2013.38:0.03000000):0.39469507,(((((((IL_04-Dec-2013_1_KJ645680.1_2013.93:0.00193567,IL_11-Dec-2013_1_KJ645689.1_2013.95:0.02193567):0.07634045,OH_24-Nov-2013_1_KJ645664.1_2013.90:0.04827613):0.44175250,(((MN_02-Nov-2013_1_KJ645652.1_2013.84:0.08159194,OH_01-Jan-2014_1_KJ645698.1_2014.00:0.24159194):0.00001442,MN_19-Nov-2013_1_KJ645660.1_2013.89:0.13160637):0.03119502,MN_20-Oct-2013_1_KJ645648.1_2013.80:0.07280139):0.31722724):0.20066886,IA_01-May-2013_86_KF468754.1_2013.34:0.13069749):0.03763314,((((OH_23-Jan-2014_1_KJ408801.1_2014.06:0.51813510,((MN_28-Nov-2013_1_KJ645667.1_2013.91:0.12175345,OH_15-Dec-2013_1_KJ645699.1_2013.96:0.17175345):0.11535848,Mexico_22-Jan-2014_1_KJ645700.1_2014.06:0.38711193):0.13102317):0.00001443,MN_11-Nov-2013_1_KJ645656.1_2013.86:0.31814952):0.35045905,IA_14-Jul-2013_1_KJ645636.1_2013.54:0.34860858):0.00001448,TX_28-Sep-2013_1_KJ645645.1_2013.74:0.54862306):0.01970756):0.00001441,(IA_06-Jun-2013_1_KF650373.1_2013.43:0.25833062,(((MN_25-Nov-2013_1_KJ645706.1_2013.90:0.14620915,MN_14-Dec-2013_1_KJ645691.1_2013.96:0.20620915):0.24359635,MN_03-Dec-2013_1_KJ645707.1_2013.92:0.40980550):0.33851069,MN_19-Jun-2013_1_KF468752.1_2013.47:0.29831619):0.00001443):0.00001442):0.21633561,IL_19-Nov-2013_1_KJ645659.1_2013.89:0.93468065):0.00001443):0.16662363):0.00001443):0.03255007):0.08068519,NC_08-Dec-2013_1_KJ645683.1_2013.94:1.26456840):0.02665174):0.45034944);"

tree <- read.tree(file="tdtree.txt")
rem <- regexpr(".......$", tree$tip.label)
tmpf <- function(x, y) substring(x, first=y)
sampleTimes <- mapply(tmpf, tree$tip.label, rem)
sampleTimes <- as.numeric(sampleTimes)

names(sampleTimes) <- tree$tip.label
bdt <- binaryDatedTree( tree, sampleTimes=sampleTimes)

births <- c( I = 'parms$beta * S * I' )
deaths <- c( I = 'parms$gamma * I' )
nonDemeDynamics <- c(S = '-parms$beta * S * I')
x0 <- c(I=1, S= unname(parms_truth$S0) )
t0 <- bdt$maxSampleTime - max(bdt$heights) -1

print(
    system.time(
        print(
            coalescent.log.likelihood(
                bdt
                , births, deaths, nonDemeDynamics
                , t0, x0
                , parms=parms_truth
                , fgyResolution = 1000
                , integrationMethod = 'rk4')
            )))

obj_fun <- function(lnbeta, lnI0){
    beta <- exp(lnbeta)
    I0 <- exp(lnI0)
    parms <- parms_truth
    parms$beta <- beta
    x0 <- c(I=unname(I0), S=unname(parms$S0))
    nll <- -coalescent.log.likelihood(
        bdt, births, deaths, nonDemeDynamics, t0, x0, parms=parms,
        fgyResolution =1000, integrationMethod='rk4')
    print(paste(nll, beta, I0))
    nll
}
library(bbmle)
fit <- mle2(obj_fun, start=list(lnbeta=log(parms_truth$beta*.75), lnI0=log(1)),
            method='Nelder-Mead', optimizer='optim', control=list(trace=6, reltol=1e-8))


## Now with two demes

beta <- 2e-4
parms_truth <- c(beta11=beta, beta12=beta, beta21=beta, beta22=beta,
                 gamma=1, S1_0=9999, S2_0=9999, t0=2012, I1_0=1, I2_0=1)
INFECTEDNAMES <- c('I1', 'I2')
births <- rbind(c('parms$beta11 * I1 * S1', 'parms$beta12 * I1 * S2'),
                c('parms$beta21 * I2 * S1', 'parms$beta22 * I2 * S2'))
rownames(births) <- colnames(births) <- INFECTEDNAMES

deaths <- c(I1='parms$gamma * I1', I2='parms$gamma * I2')

migrations <- matrix('0', 2, 2, dimnames=list(INFECTEDNAMES, INFECTEDNAMES))

nonDemeDynamics <- c(S1=paste0('-', births[1, 1], '-', births[2, 1]),
                     S2=paste0('-', births[1, 2], '-', births[2, 2]))

sampleStates <- matrix(0, nrow=length(sampleTimes), ncol=length(INFECTEDNAMES))
colnames(sampleStates) <- INFECTEDNAMES
rownames(sampleStates) <- names(sampleTimes)

isIA <- grepl('IA', names(sampleTimes))
sampleStates[!isIA, 1] <- 1
sampleStates[isIA, 2] <- 1

bdt <- binaryDatedTree(phylo=tree, sampleTimes=sampleTimes, sampleStates=sampleStates)

coalescent.log.likelihood(bdt, births, deaths, nonDemeDynamics, t0=2011,
                          x0=c(I1=1, I2=1, S1=9999, S2=9999), migrations=migrations,
                          parms=as.list(parms_truth), fgyResolution=1000, integrationMethod='rk4')

obj_fun <- function(lnbw, lnprop){
    parms <- as.list(parms_truth)
    parms$beta11 <- parms$beta22 <- exp(lnbw)
    parms$beta12 <- parms$beta21 <- exp(lnbw) * exp(lnprop)
    t0max <- bdt$maxSampleTime - bdt$maxHeight
    t0 <- t0max - 0.1
    nll <- -coalescent.log.likelihood(bdt, births, deaths, nonDemeDynamics, t0=t0,
                          x0=c(I1=1, I2=1, S1=9999, S2=9999), migrations=migrations,
                                      parms=parms, fgyResolution=1000, integrationMethod='rk4')
    #print(c(nll, exp(x[1]), exp(x[2])))
    nll
}

pars <- c(log(2e-4), log(1))
foo <- optim(pars, fn=obj, method="Nelder-Mead")

fit <- mle2(obj_fun, start=list(lnbw=log(3e-4), lnprop=log(0.2)),
            method='Nelder-Mead', optimizer='optim', control=list(trace=6, reltol=1e-8))

## simulations

birthsSim <- gsub("parms\\$", "", births)
deathsSim <- gsub("parms\\$", "", deaths)
nddSim <- gsub("parms\\$", "", nonDemeDynamics)
demo.model <- build.demographic.process(births=birthsSim, nonDemeDynamics=nddSim, migrations=migrations, deaths=deathsSim,
                                        parameterNames=c("beta11", "beta12", "beta21", "beta22", "gamma"), rcpp=TRUE, sde=TRUE)

beta <- 2e-4
theta <- c(beta11=beta, beta12=beta, beta21=beta, beta22=beta, gamma=1)
t0 <- 2012
t1 <- 2014
x0 <- c(S1=9999, S2=9999, I1=1, I2=1)
show.demographic.process(demo.model, theta=theta, x0=x0, t0=t0, t1=t1)
