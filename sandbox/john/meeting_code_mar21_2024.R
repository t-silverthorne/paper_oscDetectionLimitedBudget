# Parameters
Outloop = 100
Amp     = 1
phase   = 3
period  = 24
Nmeas   = 24
mt      = c(1:Nmeas)/Nmeas -1/Nmeas

# Option 1: Simphony
# .... generates simulated Simphony data

# Option 2: Manual white noise
Xdat = Amp*cos(2*pi*mt/period - phase) + matrix(rnorm(Outloop*Nmeas),ncol=Nmeas)

# useful to use this for fixed-period harmonic regression 
# Xdat should be matrix with same number of columns as there are measurements
# same number of rows as the size of the outer loop
# period is specified by user
DiscoRhythm::lmCSmat(Xdat,mt, period)

# for comparing to method from paper
require(devtools)
load_all()
# since you are using units of hours, you need 24/period instead of 1/period
costfun_svdpower(mt,freqs=c(24/period),Amp=Amp,alpha=.05,cfuntype='power')