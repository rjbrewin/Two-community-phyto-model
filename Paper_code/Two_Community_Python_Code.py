################# START CODE ##########################

### STEP 1: DEFINE ALL FUNCTIONS USED IN THE SCRIPT ###
# Define function for diffuse attenuation coefficient
def fcn2min_Kds(params1, X1, Y1):
    G1 = params1['G1']
    KD = params1['KD']
    model = G1+(-KD*X1)
    return(model-Y1)
# Define function for linear relationship
def fcn2min_LINREGRESS(params1, X, Y):
    P1 = params1['P1']
    P2 = params1['P2']
    model = P1+P2*X
    return(model-Y)
# Define function for Community 1 normalised Chl-a
def fcn2min_1pop(params1, X, Y):
    P1 = params1['P1']
    P2 = params1['P2']
    model = 1 - 1./(1+np.exp(-(P1/P2)*(Y-P2)))
    return(model-X)
# Define function for Community 1 and 2 normalised Chl-a
def fcn2min_2pop(params2, X, Y):
    P1 = params2['P1']
    P2 = params2['P2']
    P3 = params2['P3']
    P4 = params2['P4']
    P5 = params2['P5']
    MLD_pop = 1 - 1./(1+np.exp(-(P1/P2)*(Y-P2)))
    DCM_pop = P3*np.exp(-((Y - ((P4+(P5*3.000))))/(P5))**2.)
    model = MLD_pop + DCM_pop
    return(model-X)
# Define bbp function for Community 1
def fcn3min_bbp_norm_1pop(params2, X1, A1):
    P2 = params2['P2']
    model   = (1.-P2)*A1 + P2
    return(model-X1)
# Define bbp function for Community 1 and 2
def fcn3min_bbp_norm(params1, X1, A1, A2):
    P1 = params1['P1']
    P2 = params1['P2']
    model   = (1.-P2)*A1 + P1*A2 + P2
    return(model-X1)

### STEP 2: IMPORT PACKAGES REQUIRED FOR SCRIPT ###
import sys
import numpy as np
import numpy.ma as ma
import julian
import scipy
import scipy.integrate as spi
import matplotlib as mpl
import pandas as pd
import seawater as sw # Import CSIRO seawater package
import scipy.stats as stats
import matplotlib.ticker as tick
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from math import nan
from netCDF4 import Dataset
from matplotlib import pyplot as plt
from scipy import interpolate
from scipy import signal
from scipy.stats import anderson
from lmfit import Minimizer, Parameters
from scipy.signal import savgol_filter
from holteandtalley import HolteAndTalley
from scipy import stats
from suntime import Sun
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
# Supresses outputs when trying to "divide by zero" or "divide by NaN" or "treatment for floating-point overflow"
np.seterr(all='ignore');

### STEP 3: READ IN DATA FROM FLOAT ###
# Define float number
FLOAT_NUMBER = '6901573'
# Define path to float data (***NB: Needs to be edited depending on where data is****)
PATH_TO_FLOAT = '/Users/bobbrewin/Jupyter/Ideas/Two_comp_vertical_Chl/'
# Read the BGC-Argo netCDF file into Python
BGC_DATA = PATH_TO_FLOAT+FLOAT_NUMBER+'_Sprof.nc'
fh = Dataset(BGC_DATA, mode='r')
# Read key information on the float
PROJECT_NAME  = fh.variables['PROJECT_NAME'][:]
PI_NAME       = fh.variables['PI_NAME'][:]
PLATFORM_TYPE = fh.variables['PLATFORM_TYPE'][:]
CYCLE_NUMBER  = fh.variables['CYCLE_NUMBER'][:]
# Read time and locations of profiles
JULD      = fh.variables['JULD'][:] #Julian day of the profile.
LONGITUDE = fh.variables['LONGITUDE'][:] #Longitude degrees East
LATITUDE  = fh.variables['LATITUDE'][:] #Latitude degress North
# Read netCDF CORE ARGO  (Adjusted no data)
PRESSURE  = fh.variables['PRES'][:] #Pressure (Dbar) we assume equivalent to depth (m)
TEMP      = fh.variables['TEMP'][:] #Temperature (degrees C)
PSAL      = fh.variables['PSAL'][:] #Salinity (PSU)
# Read netCDF BGC ARGO ADJUSTED (s-profile)
CHLA      = fh.variables['CHLA_ADJUSTED'][:]   #Chlorophyll-a
BBP       = fh.variables['BBP700_ADJUSTED'][:] #Backscattering at 700nm (m^-1)
DISS_OXY  = fh.variables['DOXY_ADJUSTED'][:]   #Dissolved Oxygen (micro mole kg^-3)
PAR       = fh.variables['DOWNWELLING_PAR'][:] #PAR
# Close netCDF
fh.close()

### STEP 4: CONVERT TIME TO DATETIME ###
# Convert time from Julian day to Year-Month-Day T HOUR:MINUTE:SECOND
REFERENCE_DATE = 2433282.50000  # JULIAN DATE for reference 1950-01-01 00:00:00
TIME_BGC = np.empty(len(JULD), dtype='datetime64[s]')
for i in range(len(JULD)):
    TIME_BGC[i] = julian.from_jd(JULD[i] + REFERENCE_DATE, fmt='jd')
# Create a time matrix needed later for the contour plot (providing time for each depth level)
TIME_BGC_MATRIX = np.empty([len(TEMP[:, 0]), len(TEMP[0, :])], dtype='datetime64[s]')
for i in range(len(TEMP[0, :])):
    TIME_BGC_MATRIX[:, i] = TIME_BGC

### STEP 5: COMPUTATION OF DAILY INTEGRATED PAR ###
# Compute sunrise, sunset and daylength
Sun_rise = np.empty(len(JULD), dtype='datetime64[s]')
Sun_set = np.empty(len(JULD), dtype='datetime64[s]')
for i in range(len(JULD)):
    sun = Sun(LATITUDE[i], LONGITUDE[i])
    Sun_rise[i] = sun.get_sunrise_time(pd.to_datetime(TIME_BGC[i]))
    Sun_set[i] = sun.get_sunset_time(pd.to_datetime(TIME_BGC[i]))
DH_TIME_BGC = pd.DatetimeIndex(TIME_BGC, dtype='datetime64[ns]', name='date_time', freq=None)
DH_TIME_BGC = DH_TIME_BGC.hour + DH_TIME_BGC.minute / 60
DH_SUNRISE  = pd.DatetimeIndex(Sun_rise, dtype='datetime64[ns]', name='date_time', freq=None)
DH_SUNRISE  = DH_SUNRISE.hour + DH_SUNRISE.minute / 60
DH_SUNSET   = pd.DatetimeIndex(Sun_set, dtype='datetime64[ns]', name='date_time', freq=None)
DH_SUNSET   = DH_SUNSET.hour + DH_SUNSET.minute / 60
DH_DAYLEN   = DH_SUNSET - DH_SUNRISE
# NaN PAR data collected before sunrise and after sunset
INDEX          = np.empty([len(JULD)]) + nan
SCALE_TIME_BGC = DH_TIME_BGC - DH_SUNRISE
AFTER_TIME_BGC = DH_SUNSET - DH_TIME_BGC
ads = np.where(SCALE_TIME_BGC > 0.)
INDEX[ads] = 1.0
ads = np.where(AFTER_TIME_BGC < 0.)
INDEX[ads] = 1.0
# Convert to seconds
DH_DAYLEN = DH_DAYLEN * (60 * 60.)  # convert hours to seconds
SCALE_TIME_BGC = SCALE_TIME_BGC * (60 * 60.)  # convert hours to seconds
# Loop nan PAR data collected before sunrise and after sunset convert PAR to daily integral
NOON_PAR = PRESSURE * 0
DAY_PAR = PRESSURE * 0
for i in range(len(JULD)):
    if np.isnan(INDEX[i]):
        PAR[i, :] = PAR[i, :] * nan
        NOON_PAR[i, :] = PAR[i, :] * nan
        DAY_PAR[i, :] = PAR[i, :] * nan
    else:
        NOON_PAR[i, :] = PAR[i, :] / np.sin((np.pi * SCALE_TIME_BGC[i]) / DH_DAYLEN[i])
        DAY_PAR[i, :] = ((NOON_PAR[i, :] * 2. * DH_DAYLEN[i]) / np.pi) / 1000000.  # 1000000 converts to Einstein
# Rewrite PAR as daily integrated PAR
PAR = DAY_PAR

### STEP 6: COMPUTATION OF DEPTH, DENSITY AND BRUNT–VäISäLä USING SEAWATER PACKAGE ###
# Compute depth and density
DEPTH = PRESSURE * 0
DENSITY = PRESSURE * 0
for i in range(len(JULD)):
    DEPTH[i, :] = sw.dpth(PRESSURE[i, :], LATITUDE[i])  # use "sw.dpth" function to calculate depth from pressure
    DENSITY[i, :] = sw.dens0(PSAL[i, :], TEMP[i, :])
# Compute Brunt–Väisälä (buoyancy) frequency
BRUNT = PRESSURE * nan
for i in range(len(JULD)):
    BRUNT_T = sw.bfrq(PSAL[i, :], TEMP[i, :], PRESSURE[i, :], LATITUDE[i])
    BRUNT_T1 = BRUNT_T[0]
    BRUNT_T1 = np.resize(BRUNT_T1, len(BRUNT_T1))
    BRUNT_T3 = BRUNT_T[2]
    BRUNT_T3 = np.resize(BRUNT_T3, len(BRUNT_T3))
    PRES_T = PRESSURE[i, :]
    valid1 = np.logical_not(ma.getmask(PRES_T))
    PRES_T1 = PRES_T.compressed()
    interpfunc = interpolate.interp1d(BRUNT_T3, BRUNT_T1, kind='linear', fill_value="extrapolate")
    xxx = interpfunc(PRES_T1)
    BRUNT[i, valid1] = xxx

### STEP 7: COMPUTATION OF MIXED LAYER DEPTH ###
# Compute Mixed Layer Depth using Holt and Tally method
MLD_TEMP_ARRAY = np.empty(len(JULD)) + nan
for i in range(len(JULD)):
    A = PRESSURE[i, :]
    B = TEMP[i, :]
    valid1 = np.logical_not(ma.getmask(A))
    A = A[valid1]
    B = B[valid1]
    valid1 = np.logical_not(ma.getmask(B))
    A = A[valid1]
    B = B[valid1]
    if len(A) > 100:
        h = HolteAndTalley(A, B)
        ##The temperature algorithms mixed layer depth
        MLD_TEMP_ARRAY[i] = h.tempMLD

### STEP 8: COMPUTATION OF KD(PAR) ###
# Firstly compute Kd and Chl in top 10 m
Kd_INITAL       = np.empty([len(JULD)]) + nan
CHL_SURF_INITAL = np.empty([len(JULD)]) + nan  # top 10 m
for INDEX_PROFILE in range(len(JULD)):
    # Fitting Kd
    Y = PAR[INDEX_PROFILE, :]
    X = DEPTH[INDEX_PROFILE, :]
    valid1 = np.logical_not(ma.getmask(Y))
    X = X[valid1]
    Y = Y[valid1]
    valid2 = np.logical_not(ma.getmask(X))
    X = X[valid2]
    Y = Y[valid2]
    idx = np.where(X <= 100)  # top 100 m for comoputing Kd
    X = X[idx]
    Y = Y[idx]
    idx = np.where(Y >= 0)
    X = X[idx]
    Y = Y[idx]
    Y1 = np.log(Y)
    X1 = X
    params1 = Parameters()
    params1.add('G1', value=4.0)  # Initial guess
    params1.add('KD', value=0.05)  # Initial guess
    if len(Y1) > 10:  # has to be at least 10 data point to fit otherwise call nan
        out = Minimizer(fcn2min_Kds, params1, fcn_args=(X1, Y1))
        result = out.minimize()
        G1_FIT = result.params['G1'].value
        KD_FIT = result.params['KD'].value
    else:
        G1_FIT = nan
        KD_FIT = nan
    Kd_INITAL[INDEX_PROFILE] = KD_FIT
    # Compute surface Chl
    Y4 = CHLA[INDEX_PROFILE, :]
    X4 = DEPTH[INDEX_PROFILE, :]
    valid1 = np.logical_not(ma.getmask(Y4))
    X4 = X4[valid1]
    Y4 = Y4[valid1]
    valid2 = np.logical_not(ma.getmask(X4))
    X4 = X4[valid2]
    Y4 = Y4[valid2]
    idx = np.where(X4 <= 10.0)
    if len(idx[0]) > 0:
        CHL_SURF_INITAL[INDEX_PROFILE] = np.mean(Y4[idx])
# Quantify relationship between Surf Chl and Kd
asq = np.where(Kd_INITAL > 0.02)  # 0.02 pure water
Chl_test = CHL_SURF_INITAL[asq]
Kd_test = Kd_INITAL[asq]
asq = np.where(Chl_test > 0.02)  # 0.02 pure water
Chl_test = Chl_test[asq]
Kd_test = Kd_test[asq]
params2 = Parameters()
params2.add('P1', value=0.05)  # Initial guess
params2.add('P2', value=0.1)  # Initial guess
X1 = Chl_test
Y1 = Kd_test
out = Minimizer(fcn2min_LINREGRESS, params2, fcn_args=(X1, Y1))
result = out.minimize()
P1_FIT = result.params['P1'].value
P2_FIT = result.params['P2'].value
SLOPE_ARRAY = result.params['P2'].value
INTER_ARRAY = result.params['P1'].value
SLOPE_ARRAY_SE = result.params['P2'].stderr
INTER_ARRAY_SE = result.params['P1'].stderr
Kd_mod_stats = (P1_FIT + P2_FIT * X1)
res = Kd_mod_stats - Y1
# Compute Kd model
Kd_mod = (P1_FIT + P2_FIT * CHL_SURF_INITAL)
SLOPE_ARRAY_ST = ("{0:.3f}".format(SLOPE_ARRAY))
SLOPE_ARRAY_SE_ST = ("{0:.3f}".format(SLOPE_ARRAY_SE))
INTER_ARRAY_ST = ("{0:.3f}".format(INTER_ARRAY))
INTER_ARRAY_SE_ST = ("{0:.3f}".format(INTER_ARRAY_SE))
# Compute Pearson
p_y = SLOPE_ARRAY * X1 + INTER_ARRAY
STATS_REG = stats.pearsonr(Y1, p_y)
R_CORR_ST = ("{0:.2f}".format(STATS_REG[0]))
P_CORR_ST = ("{0:.3f}".format(STATS_REG[1]))
# Compute confidence intervals of linear regression Kd and Chl-a top 10
y_err = Y1 - p_y
p_x = np.arange(0, 1, 0.1)
mean_x = np.mean(X1)  # mean of x
n = len(X1)
t = 2.31  # appropriate t value (where n=9, two tailed 95%)
s_err = np.sum(np.power(y_err, 2))  # sum of the squares of the residuals
confs = t * np.sqrt((s_err / (n - 2)) * (1.0 / n + (np.power((p_x - mean_x), 2) /
        ((np.sum(np.power(X1, 2))) - n * (np.power(mean_x, 2))))))
# now predict y based on test x-values
p_y = SLOPE_ARRAY * p_x + INTER_ARRAY
# get lower and upper confidence limits based on predicted y and confidence intervals
lower = p_y - abs(confs)
upper = p_y + abs(confs)
# SUPPLEMENTARY FIGURE 1
# Set-up the plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), \
      gridspec_kw={'hspace': 0.05})
fig.patch.set_facecolor('White')
ax1.plot(Chl_test, Kd_test, color='b', marker='o', alpha=.2, linestyle='None', ms=10)
ax1.plot(CHL_SURF_INITAL, Kd_mod, color='k', linestyle='-', ms=10,
         label='$K_d$=' + SLOPE_ARRAY_ST + '($\pm$' + SLOPE_ARRAY_SE_ST \
               + ')$B_{s10}$+' + INTER_ARRAY_ST \
               + '($\pm$' + INTER_ARRAY_SE_ST + ') \n $r$=' + R_CORR_ST \
               + ', $p$ = ' + P_CORR_ST + '')
ax1.set_ylabel('$K_d$ (m$^{-1}$)', fontsize=15)
ax1.set_title('(a)', fontsize=25)
ax1.yaxis.set_tick_params(labelsize=12)
ax1.set_ylim([0.04, 0.085])
ax1.set_xlabel('Surface chlorophyll-a ($B_{s10}$, mg m$^{-3}$)', fontsize=15, color='k')
ax1.set_xlim([0, 1])
ax1.plot(p_x, lower, 'k--')  # ,label='Lower confidence limit (95%)')
ax1.plot(p_x, upper, 'k--')  # ,label='Upper confidence limit (95%)')
ax1.legend(loc="lower right", fontsize=12)
ax2.hist(res, color='b', alpha=.5)
ax2.set_ylabel('Frequency', fontsize=15)
ax2.set_xlabel('$K_d$ model $-$ $K_d$ Data (m$^{-1}$)', fontsize=15)
ax2.set_title('(b)', fontsize=25)
ax2.set_xlim([-0.015, 0.015])
plt.show()
# Normality test for linear Kd
result = anderson(res)
print('Statistic: %.3f' % result.statistic)
p = 0
for i in range(len(result.critical_values)):
    sl, cv = result.significance_level[i], result.critical_values[i]
    if result.statistic < result.critical_values[i]:
        print('%.3f: %.3f, data looks normal (fail to reject H0)' % (sl, cv))
    else:
        print('%.3f: %.3f, data does not look normal (reject H0)' % (sl, cv))
# Fill Kd in with Kdmod where not available remove any Kd <= pure water
asq            = np.argwhere(np.isnan(Kd_INITAL)) #nans
Kd_INITAL[asq] = Kd_mod[asq]
asq            = np.where(Kd_INITAL <= 0.02) #0.02 pure water
Kd_INITAL[asq] = Kd_mod[asq]
Kd_ALL         = Kd_INITAL

### STEP 9: COMPUTATION OF Zp, integrated Chl-a and optical CHl-a and bbp) ###
CHL_SAT_ALL = np.empty([len(JULD)]) + nan
CHL_SAT_SD = np.empty([len(JULD)]) + nan
CHL_INT_1_ALL = np.empty([len(JULD)]) + nan
CHL_INT_2_ALL = np.empty([len(JULD)]) + nan
FIRST_OD_ALL = np.empty([len(JULD)]) + nan
EUPHOTIC_D_ALL = np.empty([len(JULD)]) + nan
BBP_SAT_ALL = np.empty([len(JULD)]) + nan
BBP_SAT_SD = np.empty([len(JULD)]) + nan
for INDEX_PROFILE in range(len(JULD)):
    # Compute 1st optical depth and euphotic depth from Kd
    FIRST_OD_ALL[INDEX_PROFILE] = 1.0 / Kd_ALL[INDEX_PROFILE]
    EUPHOTIC_D_ALL[INDEX_PROFILE] = 4.6 / Kd_ALL[INDEX_PROFILE]
    # Compute satellite(1st optical depth) Chl and column integrated
    Y4 = CHLA[INDEX_PROFILE, :]
    X4 = DEPTH[INDEX_PROFILE, :]
    valid1 = np.logical_not(ma.getmask(Y4))
    X4 = X4[valid1]
    Y4 = Y4[valid1]
    valid2 = np.logical_not(ma.getmask(X4))
    X4 = X4[valid2]
    Y4 = Y4[valid2]
    idx = np.where(X4 <= FIRST_OD_ALL[INDEX_PROFILE])
    if len(idx[0]) > 0:
        CHL_TEMP = Y4[idx]
        asd = np.where(CHL_TEMP >= 0.01)
        CHL_SAT_ALL[INDEX_PROFILE] = np.median(CHL_TEMP[asd])
        CHL_SAT_SD[INDEX_PROFILE] = np.std(CHL_TEMP[asd])
    idx = np.where(X4 <= EUPHOTIC_D_ALL[INDEX_PROFILE])
    if len(idx[0]) > 0:
        CHL_INT_1_ALL[INDEX_PROFILE] = spi.trapz(Y4[idx], X4[idx])
    idx = np.where(X4 <= EUPHOTIC_D_ALL[INDEX_PROFILE] * 1.5)
    if len(idx[0]) > 0:
        CHL_INT_2_ALL[INDEX_PROFILE] = spi.trapz(Y4[idx], X4[idx])
    # Compute satellite(1st optical depth) BBP
    Y4 = BBP[INDEX_PROFILE, :]
    X4 = DEPTH[INDEX_PROFILE, :]
    valid1 = np.logical_not(ma.getmask(Y4))
    X4 = X4[valid1]
    Y4 = Y4[valid1]
    valid2 = np.logical_not(ma.getmask(X4))
    X4 = X4[valid2]
    Y4 = Y4[valid2]
    idx = np.where(X4 <= FIRST_OD_ALL[INDEX_PROFILE])
    if len(idx[0]) > 0:
        BBP_TEMP = Y4[idx]
        asd = np.where(BBP_TEMP > 0.0)
        BBP_SAT_ALL[INDEX_PROFILE] = np.median(BBP_TEMP[asd])
        BBP_SAT_SD[INDEX_PROFILE] = np.std(BBP_TEMP[asd])

### STEP 10: PLOT MAP OF FLOAT DATA (FIGURE 1 OF PAPER) ###
# FIGURE 1 of paper (Map of floats)
XSIZE = 10  # Define the xsize of the figure window
YSIZE = 7  # Define the ysize of the figure window
Title_font_size = 7.5  # Define the font size of the titles
Label_font_size = 6  # Define the font size of the labels
cmap = mpl.cm.copper  # Defines the colour scale
cmap2 = mpl.cm.get_cmap('copper_r')  # Defines the colour scale
# initialise figure
fig = plt.figure(figsize=(XSIZE, YSIZE))
fig.patch.set_facecolor('White')
ax = fig.add_subplot(111, projection=ccrs.PlateCarree(central_longitude=35))
transform = ccrs.PlateCarree()._as_mpl_transform(ax)
# Loop through profiles and colour according to time
dates = [pd.to_datetime(d) for d in TIME_BGC]
for i in range(len(JULD)):
    ax.plot(LONGITUDE[i], LATITUDE[i], color=cmap(i / float(len(JULD))), marker='o', markersize=6,
            transform=ccrs.PlateCarree())
idex_label = [0, 50, 87, 122, 138, 158, 177]
for i in range(len(idex_label)):
    j = idex_label[i]
    ax.annotate(dates[j].strftime("%Y-%m-%d"),
                xy=(LONGITUDE[j] + 0.1, LATITUDE[j]),
                xycoords=transform, fontsize=12)
lbl_color = 'k'
# dd continents
ax.add_feature(
    cfeature.NaturalEarthFeature('physical', 'land', '10m', linewidth=0.5, edgecolor='gray', facecolor='whitesmoke'))
# set geographical extent of figure
ax.set_extent([32, 40, 22, 30], ccrs.PlateCarree())
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.375, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.bottom_labels = True
gl.left_labels = True
gl.right_labels = False
gl.xlines = True
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.ylabel_style = {'fontsize': 15, 'color': lbl_color}
gl.xlabel_style = {'fontsize': 15, 'color': lbl_color}
ax.text(-0.17, 0.55, 'Latitude', va='bottom', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax.transAxes, fontsize=15)
ax.text(0.5, -0.12, 'Longitude', va='bottom', ha='center',
        rotation='horizontal', rotation_mode='anchor',
        transform=ax.transAxes, fontsize=15)
# plot colorbar
divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="5%", pad=0.3, axes_class=plt.Axes)
fig.add_axes(ax_cb)
cb = plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap2), cax=ax_cb, ticks=[0, 1])
cb.ax.set_yticklabels(['Float-finish', 'Float-start'], fontsize=12)
plt.show()

### STEP 11: DEFINE NUMBER OF BOOTSTRAPS FOR RUNS ###
### NB: TO RECREATE PAPER RESULTS NEED TO BE 1000 ###
### USING 10 HERE AS RUN OF 1000 TAKES A LOT OF TIME
BOOTIT = 10#00

### STEP 12: COMPUTE RELATIONSHIP BETWEEN PARAMETERS OF COMMUNITY 1 for Zp<Zm WHEN THEY DOMINATE ##
# Create array of parameters
P1_ARRAY       = np.empty([len(JULD)])+nan
P1_ARRAY_LB    = np.empty([len(JULD)])+nan
P1_ARRAY_UB    = np.empty([len(JULD)])+nan
S1_ARRAY       = np.empty([len(JULD)])+nan
S1_ARRAY_LB    = np.empty([len(JULD)])+nan
S1_ARRAY_UB    = np.empty([len(JULD)])+nan
TAU1_ARRAY     = np.empty([len(JULD)])+nan
TAU1_ARRAY_LB  = np.empty([len(JULD)])+nan
TAU1_ARRAY_UB  = np.empty([len(JULD)])+nan
#Loop through each profile
for INDEX_PROFILE in range(len(JULD)):
    # Cases of mixed waters
    if EUPHOTIC_D_ALL[INDEX_PROFILE]/MLD_TEMP_ARRAY[INDEX_PROFILE] < 1.0:
        # Data for fitting
        CHL_DIM = CHLA[INDEX_PROFILE,:]/CHL_SAT_ALL[INDEX_PROFILE]
        OPT_DIM = DEPTH[INDEX_PROFILE,:]*Kd_ALL[INDEX_PROFILE]
        DEP_TEM = DEPTH[INDEX_PROFILE,:]
        MLD_OD  = MLD_TEMP_ARRAY[INDEX_PROFILE]*Kd_ALL[INDEX_PROFILE]
        #Create an array of parameters
        P1_TEMP      = np.empty([BOOTIT])+nan
        S1_TEMP      = np.empty([BOOTIT])+nan
        TAU1_TEMP    = np.empty([BOOTIT])+nan
        # Process data to be used for Chl fit
        X1        = CHL_DIM
        Y1        = OPT_DIM
        valid1    = np.logical_not(ma.getmask(X1))
        X1        = X1[valid1]
        Y1        = Y1[valid1]
        valid1    = np.logical_not(ma.getmask(Y1))
        X1        = X1[valid1]
        Y1        = Y1[valid1]
        if MLD_OD < 9.2:
            asd = np.where(Y1 < 9.2)
            X1        = X1[asd]
            Y1        = Y1[asd]
        for j in range(BOOTIT):
            RANDOM_IND = np.random.choice(len(Y1), len(Y1))
            X22 = X1[RANDOM_IND]
            Y22 = Y1[RANDOM_IND]
            params1  = Parameters()
            params1.add('P1', value=9., min = 4.6, max = 100.)
            params1.add('P2', value=MLD_OD)
            out      = Minimizer(fcn2min_1pop, params1, fcn_args=(X22, Y22))
            result   = out.minimize() #
            P1_FIT   = result.params['P1'].value
            P2_FIT   = result.params['P2'].value
            MLD_pop  = 1 - 1./(1+np.exp(-(P1_FIT/P2_FIT)*(Y22-P2_FIT)))
            r        = np.corrcoef(X22, MLD_pop)
            if r[1,0]**2 >= 0.90:
                P1_TEMP[j]   = result.params['P1'].value
                S1_TEMP[j]   = result.params['P1'].value/result.params['P2'].value
                TAU1_TEMP[j] = result.params['P2'].value
        # Compute median and percentiles of bootstrap
        UB_PERC = 97.5
        LB_PERC = 2.5
        P1_ARRAY[INDEX_PROFILE]        = np.nanmedian(P1_TEMP)
        P1_ARRAY_LB[INDEX_PROFILE]     = np.nanpercentile(P1_TEMP,LB_PERC)
        P1_ARRAY_UB[INDEX_PROFILE]     = np.nanpercentile(P1_TEMP,UB_PERC)
        S1_ARRAY[INDEX_PROFILE]        = np.nanmedian(S1_TEMP)
        S1_ARRAY_LB[INDEX_PROFILE]     = np.nanpercentile(S1_TEMP,LB_PERC)
        S1_ARRAY_UB[INDEX_PROFILE]     = np.nanpercentile(S1_TEMP,UB_PERC)
        TAU1_ARRAY[INDEX_PROFILE]      = np.nanmedian(TAU1_TEMP)
        TAU1_ARRAY_LB[INDEX_PROFILE]   = np.nanpercentile(TAU1_TEMP,LB_PERC)
        TAU1_ARRAY_UB[INDEX_PROFILE]   = np.nanpercentile(TAU1_TEMP,UB_PERC)
# Fit model Tau vrs Kd*Zm and get coefficients of linear fit and standard errors
Y1 = TAU1_ARRAY
X1 = MLD_TEMP_ARRAY * Kd_ALL
valid1 = np.where(X1 > 0)
X1 = X1[valid1]
Y1 = Y1[valid1]
valid1 = np.where(Y1 > 0)
x = X1[valid1]
y = Y1[valid1]
if len(y) > 3:
    params1 = Parameters()
    params1.add('P1', value=0.0)  # , vary=False)
    params1.add('P2', value=1.0)
    out = Minimizer(fcn2min_LINREGRESS, params1, fcn_args=(x, y))
    result = out.minimize()  # use powell minimisation method
    SLOPE_ARRAY = result.params['P2'].value
    INTER_ARRAY = result.params['P1'].value
    SLOPE_ARRAY_SE = result.params['P2'].stderr
    INTER_ARRAY_SE = result.params['P1'].stderr
SLOPE_ARRAY_ST = ("{0:.2f}".format(SLOPE_ARRAY))
SLOPE_ARRAY_SE_ST = ("{0:.2f}".format(SLOPE_ARRAY_SE))
INTER_ARRAY_ST = ("{0:.2f}".format(INTER_ARRAY))
INTER_ARRAY_SE_ST = ("{0:.2f}".format(INTER_ARRAY_SE))
# Pearson
p_y       = SLOPE_ARRAY * x + INTER_ARRAY
STATS_REG = stats.pearsonr(y, p_y)
R_CORR_ST = ("{0:.2f}".format(STATS_REG[0]))
P_CORR_ST = ("{0:.3f}".format(STATS_REG[1]))
####Plot linear regression with confidence intervals
y_err = y - p_y
p_x = np.arange(0, 20, 1)
mean_x = np.mean(x)
n = len(x)
t = 2.31  # appropriate t value (where n=9, two tailed 95%)
s_err = np.sum(np.power(y_err, 2))  # sum of the squares of the residuals
confs = t * np.sqrt((s_err / (n - 2)) * (1.0 / n + (np.power((p_x - mean_x), 2) /
        ((np.sum(np.power(x, 2))) - n * (np.power(mean_x, 2))))))
p_y = SLOPE_ARRAY * p_x + INTER_ARRAY
# get lower and upper confidence limits based on predicted y and confidence intervals
lower = p_y - abs(confs)
upper = p_y + abs(confs)
# # SUPPLEMENTARY FIGURE 2
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharex=True, \
                               gridspec_kw={'hspace': 0.05})
fig.patch.set_facecolor('White')
# Temperature subplot
ax1.plot(X1, Y1, color='b', marker='o', linestyle='None', alpha=.2, ms=10)
ax1.set_ylabel('$\\tau_1$', fontsize=15)
ax1.yaxis.set_tick_params(labelsize=12)
ax1.set_ylim([0, 20])
ax1.set_title('(a)', fontsize=25)
ax1.set_xlabel('Z$_mK_d$', fontsize=15, color='k')
# ax1.xaxis.set_tick_params(labelsize=12)
ax1.set_xlim([0, 20])
# plot line of best fit
ax1.plot(p_x, p_y, 'k-', label='$\\tau_1$=' + SLOPE_ARRAY_ST + '($\pm$' + SLOPE_ARRAY_SE_ST \
                               + ')$\\times$Z$_m\\times$K$_d$+' + INTER_ARRAY_ST \
                               + '($\pm$' + INTER_ARRAY_SE_ST + ') \n $r$=' + R_CORR_ST \
                               + ', $p$ = ' + P_CORR_ST + '')
ax1.plot(p_x, lower, 'k--', label='Lower confidence limit (95%)')
ax1.plot(p_x, upper, 'k--', label='Upper confidence limit (95%)')
ax1.legend(loc="upper left", fontsize=12)
p_y2 = p_y
####Fit model P1 and Tau get coefficients of linear fit and standard errors
Y1 = P1_ARRAY
X1 = TAU1_ARRAY
valid1 = np.where(X1 > 0)
X1 = X1[valid1]
Y1 = Y1[valid1]
valid1 = np.where(Y1 > 0)
x = X1[valid1]
y = np.log10(Y1[valid1])
params1 = Parameters()
params1.add('P1', value=np.log10(4.6), vary=False)
params1.add('P2', value=1.0)
out = Minimizer(fcn2min_LINREGRESS, params1, fcn_args=(x, y))
result = out.minimize()  # use powell minimisation method
SLOPE_ARRAY_P1 = result.params['P2'].value
INTER_ARRAY_P1 = result.params['P1'].value
SLOPE_ARRAY_P1_SE = result.params['P2'].stderr
INTER_ARRAY_P1_SE = result.params['P1'].stderr
SLOPE_ARRAY_ST = ("{0:.2f}".format(SLOPE_ARRAY_P1))
SLOPE_ARRAY_SE_ST = ("{0:.2f}".format(SLOPE_ARRAY_P1_SE))
INTER_ARRAY_ST = ("{0:.2f}".format(INTER_ARRAY_P1))
INTER_ARRAY_SE_ST = ("{0:.2f}".format(INTER_ARRAY_P1_SE))
p_y = SLOPE_ARRAY_P1 * x + INTER_ARRAY_P1
# Pearson
STATS_REG_P1 = stats.pearsonr(y, p_y)
STATS_REG_P1 = stats.pearsonr(y, p_y)
R_CORR_ST = ("{0:.2f}".format(STATS_REG_P1[0]))
P_CORR_ST = ("{0:.3f}".format(STATS_REG_P1[1]))
# calculate the y-error (residuals)
y_err = y - p_y
p_x = np.arange(0, 20, 1)
mean_x = np.mean(x)
n = len(x)
t = 2.31  # appropriate t value (where n=9, two tailed 95%)
s_err = np.sum(np.power(y_err, 2))  # sum of the squares of the residuals
confs = t * np.sqrt((s_err / (n - 2)) * (1.0 / n + (np.power((p_x - mean_x), 2) /
                                                    ((np.sum(np.power(x, 2))) - n * (np.power(mean_x, 2))))))
# now predict y based on test x-values
p_y = SLOPE_ARRAY_P1 * p_x + INTER_ARRAY_P1
# get lower and upper confidence limits based on predicted y and confidence intervals
lower = p_y - abs(confs)
upper = p_y + abs(confs)
ax2.plot(x, 10 ** (y), color='b', marker='o', linestyle='None', alpha=.2, ms=10)
ax2.plot(p_x, 10 ** (p_y), 'k-', label='$P_1$ = 10$^{' + SLOPE_ARRAY_ST + '(\pm' + SLOPE_ARRAY_SE_ST \
                                       + ')\\times \\tau_1 + ' + INTER_ARRAY_ST + '(\pm' + SLOPE_ARRAY_SE_ST \
                                       + ')}$ \n $r$=' + R_CORR_ST + ', $p$ = ' + P_CORR_ST + '')
ax2.plot(p_x, 10 ** (lower), 'k--')  # ,label='Lower confidence limit (95%)')
ax2.plot(p_x, 10 ** (upper), 'k--')  # ,label='Upper confidence limit (95%)')
ax2.set_ylabel('$P_1$', fontsize=15)
ax2.set_title('(b)', fontsize=25)
ax2.yaxis.set_tick_params(labelsize=12)
ax2.set_ylim([1, 1000])
ax2.set_yscale('log')
ax2.set_xlabel('$\\tau_1$', fontsize=15, color='k')
ax2.xaxis.set_tick_params(labelsize=12)
ax2.set_xlim([0, 20])
ax2.legend(loc="upper left", fontsize=12)
plt.show()

### STEP 13: FIT MODEL TO DATA ###
# Create array of parameters
P1_ARRAY       = np.empty([len(JULD)])+nan
P1_ARRAY_LB    = np.empty([len(JULD)])+nan
P1_ARRAY_UB    = np.empty([len(JULD)])+nan
S1_ARRAY       = np.empty([len(JULD)])+nan
S1_ARRAY_LB    = np.empty([len(JULD)])+nan
S1_ARRAY_UB    = np.empty([len(JULD)])+nan
TAU1_ARRAY     = np.empty([len(JULD)])+nan
TAU1_ARRAY_LB  = np.empty([len(JULD)])+nan
TAU1_ARRAY_UB  = np.empty([len(JULD)])+nan
BM2_ARRAY      = np.empty([len(JULD)])+nan
BM2_ARRAY_LB   = np.empty([len(JULD)])+nan
BM2_ARRAY_UB   = np.empty([len(JULD)])+nan
TAU2_ARRAY     = np.empty([len(JULD)])+nan
TAU2_ARRAY_LB  = np.empty([len(JULD)])+nan
TAU2_ARRAY_UB  = np.empty([len(JULD)])+nan
SIG2_ARRAY     = np.empty([len(JULD)])+nan
SIG2_ARRAY_LB  = np.empty([len(JULD)])+nan
SIG2_ARRAY_UB  = np.empty([len(JULD)])+nan
bbpS1_ARRAY    = np.empty([len(JULD)])+nan
bbpS1_ARRAY_LB = np.empty([len(JULD)])+nan
bbpS1_ARRAY_UB = np.empty([len(JULD)])+nan
bbpS2_ARRAY    = np.empty([len(JULD)])+nan
bbpS2_ARRAY_LB = np.empty([len(JULD)])+nan
bbpS2_ARRAY_UB = np.empty([len(JULD)])+nan
bbpk_ARRAY     = np.empty([len(JULD)])+nan
bbpk_ARRAY_LB  = np.empty([len(JULD)])+nan
bbpk_ARRAY_UB  = np.empty([len(JULD)])+nan
bbpST1_ARRAY    = np.empty([len(JULD)])+nan
bbpST1_ARRAY_LB = np.empty([len(JULD)])+nan
bbpST1_ARRAY_UB = np.empty([len(JULD)])+nan
bbpST2_ARRAY    = np.empty([len(JULD)])+nan
bbpST2_ARRAY_LB = np.empty([len(JULD)])+nan
bbpST2_ARRAY_UB = np.empty([len(JULD)])+nan
bbpSk_ARRAY     = np.empty([len(JULD)])+nan
bbpSk_ARRAY_LB  = np.empty([len(JULD)])+nan
bbpSk_ARRAY_UB  = np.empty([len(JULD)])+nan
for INDEX_PROFILE in range(len(JULD)):
    # Data for fitting
    BBP_DIM  = BBP[INDEX_PROFILE,:]/BBP_SAT_ALL[INDEX_PROFILE]
    CHL_DIM  = CHLA[INDEX_PROFILE,:]/CHL_SAT_ALL[INDEX_PROFILE]
    OPT_DIM  = DEPTH[INDEX_PROFILE,:]*Kd_ALL[INDEX_PROFILE]
    OPT_DIM2 = OPT_DIM
    DEP_TEM  = DEPTH[INDEX_PROFILE,:]
    MLD_OD   = MLD_TEMP_ARRAY[INDEX_PROFILE]*Kd_ALL[INDEX_PROFILE]
    # Only select data <9.2 Optical depths if MLD_opt < 9.2
    if MLD_OD < 9.2:
        valid1  = np.where(OPT_DIM < 9.2)
        CHL_DIM = CHL_DIM[valid1]
        OPT_DIM = OPT_DIM[valid1]
        DEP_TEM = DEP_TEM[valid1]
    # Top 500m for bbp
    valid1 = np.where(DEP_TEM < 500)
    BBP_DIM  = BBP_DIM[valid1]
    OPT_DIM2 = OPT_DIM2[valid1]
    #Create an array of parameters
    P1_TEMP      = np.empty([BOOTIT])+nan
    S1_TEMP      = np.empty([BOOTIT])+nan
    TAU1_TEMP    = np.empty([BOOTIT])+nan
    BM2_TEMP     = np.empty([BOOTIT])+nan
    TAU2_TEMP    = np.empty([BOOTIT])+nan
    SIG2_TEMP    = np.empty([BOOTIT])+nan
    bbpST1_TEMP  = np.empty([BOOTIT])+nan
    bbpST2_TEMP  = np.empty([BOOTIT])+nan
    bbpSk_TEMP   = np.empty([BOOTIT])+nan
    bbpS1_TEMP   = np.empty([BOOTIT])+nan
    bbpS2_TEMP   = np.empty([BOOTIT])+nan
    bbpk_TEMP    = np.empty([BOOTIT])+nan
    for j in range(BOOTIT):
        RANDOM_IND    = np.random.choice(len(OPT_DIM), len(OPT_DIM))
        CHL_DIM_TEMP  = CHL_DIM[RANDOM_IND]
        OPT_DIM_TEMP1 = OPT_DIM[RANDOM_IND]
        RANDOM_IND    = np.random.choice(len(OPT_DIM2), len(OPT_DIM2))
        BBP_DIM_TEMP  = BBP_DIM[RANDOM_IND]
        OPT_DIM_TEMP2 = OPT_DIM2[RANDOM_IND]
        # Process data to be used for Chl fit
        X1        = CHL_DIM_TEMP
        Y1        = OPT_DIM_TEMP1
        valid1    = np.logical_not(ma.getmask(X1))
        X1        = X1[valid1]
        Y1        = Y1[valid1]
        valid1    = np.logical_not(ma.getmask(Y1))
        X1        = X1[valid1]
        Y1        = Y1[valid1]
        if len(X1) > 6:
            # Fit 1st population
            params1  = Parameters()
            params1.add('P1', value=9., min = 4.6, max = 100.)
            params1.add('P2', value=MLD_OD)
            out      = Minimizer(fcn2min_1pop, params1, fcn_args=(X1, Y1))
            result   = out.minimize() #  use powell minimisation method ####method = 'powell'
            P1_FIT   = result.params['P1'].value
            P2_FIT   = result.params['P2'].value
            AIC_FIT1 = result.aic
            CHI_FIT1 = result.chisqr
            MLD_pop  = 1 - 1./(1+np.exp(-(P1_FIT/P2_FIT)*(Y1-P2_FIT)))
            r        = np.corrcoef(X1, MLD_pop)
            # Fit 2nd population
            params2 = Parameters()
            if r[1,0]**2 >= 0.90:
                P3_FIT = nan
                P4_FIT = nan
                P5_FIT = nan
            else:
                Tau1_temp = INTER_ARRAY+(MLD_OD*SLOPE_ARRAY)
                P1_temp   = 10**(SLOPE_ARRAY_P1 * Tau1_temp + INTER_ARRAY_P1)
                params2.add('P1', value=P1_temp, vary=False)#, min = 4.6, max = 8.)
                params2.add('P2', value=Tau1_temp, vary=False)#, min=1.0, max = 3.0)
                params2.add('P3', value=5.0, min = 0.0, max = 100.0)
                params2.add('P4', value=4.6, min = 0.0)
                params2.add('P5', value=1.0, min = 0.0)
                res     = Minimizer(fcn2min_2pop,  params2, fcn_args=(X1, Y1))
                result  = res.minimize()
                AIC_FIT2 = result.aic
                CHI_FIT2 = result.chisqr
                if AIC_FIT2 < AIC_FIT1:
                    P1_FIT = result.params['P1'].value
                    P2_FIT = result.params['P2'].value
                    P3_FIT = result.params['P3'].value
                    P4_FIT = result.params['P4'].value
                    P5_FIT = result.params['P5'].value
                else:
                    P3_FIT = nan
                    P4_FIT = nan
                    P5_FIT = nan
        else:
            P1_FIT = nan
            P2_FIT = nan
            P3_FIT = nan
            P4_FIT = nan
            P5_FIT = nan
        P1_TEMP[j]   = P1_FIT
        S1_TEMP[j]   = P1_FIT/P2_FIT
        TAU1_TEMP[j] = P2_FIT
        BM2_TEMP[j]  = P3_FIT
        TAU2_TEMP[j] = P4_FIT + P5_FIT * 3.000
        SIG2_TEMP[j] = P5_FIT
        # Process data to be used for bbp fit
        X1        = BBP_DIM_TEMP
        Y1        = OPT_DIM_TEMP2
        A1        = 1 - 1./(1+np.exp(-(P1_TEMP[j]/TAU1_TEMP[j])*(Y1-TAU1_TEMP[j])))
        A2        = BM2_TEMP[j]*np.exp(-((Y1 - TAU2_TEMP[j])/SIG2_TEMP[j])**2.)
        valid1    = np.logical_not(ma.getmask(X1))
        X1        = X1[valid1]
        Y1        = Y1[valid1]
        A1        = A1[valid1]
        A2        = A2[valid1]
        valid1    = np.logical_not(ma.getmask(Y1))
        X1        = X1[valid1]
        Y1        = Y1[valid1]
        A1        = A1[valid1]
        A2        = A2[valid1]
        if len(X1) > 6:
            if np.isnan(BM2_TEMP[j]):
                params1.add('P2', value=0.3, min = 0.01, max = 0.95)
                out     = Minimizer(fcn3min_bbp_norm_1pop, params1, fcn_args=(X1, A1))
                result  = out.minimize()
                bbpST1_TEMP[j]  = 1 - result.params['P2'].value
                bbpST2_TEMP[j]  = nan
                bbpSk_TEMP[j]   = result.params['P2'].value
                bbpS1_TEMP[j]   = bbpST1_TEMP[j]/(CHL_SAT_ALL[INDEX_PROFILE]/BBP_SAT_ALL[INDEX_PROFILE])
                bbpS2_TEMP[j]   = bbpST2_TEMP[j]/(CHL_SAT_ALL[INDEX_PROFILE]/BBP_SAT_ALL[INDEX_PROFILE])
                bbpk_TEMP[j]    = bbpSk_TEMP[j]*BBP_SAT_ALL[INDEX_PROFILE]
            else:
                params1.add('P1', value=0.2, min = 0.01)
                params1.add('P2', value=0.3, min = 0.01, max = 0.95)
                out     = Minimizer(fcn3min_bbp_norm, params1, fcn_args=(X1, A1, A2))
                result  = out.minimize() # use powell minimisation method  #method = 'powell'
                bbpST1_TEMP[j]  = 1 - result.params['P2'].value
                bbpST2_TEMP[j]  = result.params['P1'].value
                bbpSk_TEMP[j]   = result.params['P2'].value
                bbpS1_TEMP[j]   = bbpST1_TEMP[j]/(CHL_SAT_ALL[INDEX_PROFILE]/BBP_SAT_ALL[INDEX_PROFILE])
                bbpS2_TEMP[j]   = bbpST2_TEMP[j]/(CHL_SAT_ALL[INDEX_PROFILE]/BBP_SAT_ALL[INDEX_PROFILE])
                bbpk_TEMP[j]    = bbpSk_TEMP[j]*BBP_SAT_ALL[INDEX_PROFILE]
        else:
            bbpST1_TEMP[j]  = nan
            bbpST2_TEMP[j]  = nan
            bbpSk_TEMP[j]   = nan
            bbpS1_TEMP[j]   = nan
            bbpS2_TEMP[j]   = nan
            bbpk_TEMP[j]    = nan
    # Compute median and percentiles of bootstrap
    UB_PERC = 97.5
    LB_PERC = 2.5
    P1_ARRAY[INDEX_PROFILE]        = np.nanmedian(P1_TEMP)
    P1_ARRAY_LB[INDEX_PROFILE]     = np.nanpercentile(P1_TEMP,LB_PERC)
    P1_ARRAY_UB[INDEX_PROFILE]     = np.nanpercentile(P1_TEMP,UB_PERC)
    S1_ARRAY[INDEX_PROFILE]        = np.nanmedian(S1_TEMP)
    S1_ARRAY_LB[INDEX_PROFILE]     = np.nanpercentile(S1_TEMP,LB_PERC)
    S1_ARRAY_UB[INDEX_PROFILE]     = np.nanpercentile(S1_TEMP,UB_PERC)
    TAU1_ARRAY[INDEX_PROFILE]      = np.nanmedian(TAU1_TEMP)
    TAU1_ARRAY_LB[INDEX_PROFILE]   = np.nanpercentile(TAU1_TEMP,LB_PERC)
    TAU1_ARRAY_UB[INDEX_PROFILE]   = np.nanpercentile(TAU1_TEMP,UB_PERC)
    BM2_ARRAY[INDEX_PROFILE]       = np.nanmedian(BM2_TEMP)
    BM2_ARRAY_LB[INDEX_PROFILE]    = np.nanpercentile(BM2_TEMP,LB_PERC)
    BM2_ARRAY_UB[INDEX_PROFILE]    = np.nanpercentile(BM2_TEMP,UB_PERC)
    TAU2_ARRAY[INDEX_PROFILE]      = np.nanmedian(TAU2_TEMP)
    TAU2_ARRAY_LB[INDEX_PROFILE]   = np.nanpercentile(TAU2_TEMP,LB_PERC)
    TAU2_ARRAY_UB[INDEX_PROFILE]   = np.nanpercentile(TAU2_TEMP,UB_PERC)
    SIG2_ARRAY[INDEX_PROFILE]      = np.nanmedian(SIG2_TEMP)
    SIG2_ARRAY_LB[INDEX_PROFILE]   = np.nanpercentile(SIG2_TEMP,LB_PERC)
    SIG2_ARRAY_UB[INDEX_PROFILE]   = np.nanpercentile(SIG2_TEMP,UB_PERC)
    bbpST1_ARRAY[INDEX_PROFILE]    = np.nanmedian(bbpST1_TEMP)
    bbpST1_ARRAY_LB[INDEX_PROFILE] = np.nanpercentile(bbpST1_TEMP,LB_PERC)
    bbpST1_ARRAY_UB[INDEX_PROFILE] = np.nanpercentile(bbpST1_TEMP,UB_PERC)
    bbpST2_ARRAY[INDEX_PROFILE]    = np.nanmedian(bbpST2_TEMP)
    bbpST2_ARRAY_LB[INDEX_PROFILE] = np.nanpercentile(bbpST2_TEMP,LB_PERC)
    bbpST2_ARRAY_UB[INDEX_PROFILE] = np.nanpercentile(bbpST2_TEMP,UB_PERC)
    bbpSk_ARRAY[INDEX_PROFILE]     = np.nanmedian(bbpSk_TEMP)
    bbpSk_ARRAY_LB[INDEX_PROFILE]  = np.nanpercentile(bbpSk_TEMP,LB_PERC)
    bbpSk_ARRAY_UB[INDEX_PROFILE]  = np.nanpercentile(bbpSk_TEMP,UB_PERC)
    bbpS1_ARRAY[INDEX_PROFILE]     = np.nanmedian(bbpS1_TEMP)
    bbpS1_ARRAY_LB[INDEX_PROFILE]  = np.nanpercentile(bbpS1_TEMP,LB_PERC)
    bbpS1_ARRAY_UB[INDEX_PROFILE]  = np.nanpercentile(bbpS1_TEMP,UB_PERC)
    bbpS2_ARRAY[INDEX_PROFILE]     = np.nanmedian(bbpS2_TEMP)
    bbpS2_ARRAY_LB[INDEX_PROFILE]  = np.nanpercentile(bbpS2_TEMP,LB_PERC)
    bbpS2_ARRAY_UB[INDEX_PROFILE]  = np.nanpercentile(bbpS2_TEMP,UB_PERC)
    bbpk_ARRAY[INDEX_PROFILE]      = np.nanmedian(bbpk_TEMP)
    bbpk_ARRAY_LB[INDEX_PROFILE]   = np.nanpercentile(bbpk_TEMP,LB_PERC)
    bbpk_ARRAY_UB[INDEX_PROFILE]   = np.nanpercentile(bbpk_TEMP,UB_PERC)
    print(['Profile fit:', str(INDEX_PROFILE)])

### STEP 14: SUPPLEMENTARY FIGURE 6: TIMESERIES OF PARAMETERS ###
fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7,) = plt.subplots(7, sharex=True, figsize=(7, 7), \
      gridspec_kw={'hspace': 0.1})
fig.patch.set_facecolor('White')
# P1
ax1.errorbar(TIME_BGC, P1_ARRAY, yerr=(P1_ARRAY - P1_ARRAY_LB, P1_ARRAY_UB - P1_ARRAY), ms=2, linestyle='None',
             marker='o', color='k', alpha=.5)
ax1.set_ylabel('$P_1$', fontsize=Title_font_size, color='k')
ax1.yaxis.set_tick_params(labelsize=6)
ax1.set_ylim(0, 120)
# tau 1
ax2.errorbar(TIME_BGC, TAU1_ARRAY, yerr=(TAU1_ARRAY - TAU1_ARRAY_LB, TAU1_ARRAY_UB - TAU1_ARRAY), ms=2,
             linestyle='None',
             marker='o', color='b', alpha=.5)
ax2.set_ylabel('$\\tau_1$', fontsize=Title_font_size, color='b')
ax2.yaxis.set_tick_params(labelsize=6)
ax2.set_ylim(0, 20)
# B2m*
ax3.errorbar(TIME_BGC, BM2_ARRAY, yerr=(BM2_ARRAY - BM2_ARRAY_LB, BM2_ARRAY_UB - BM2_ARRAY), ms=2, linestyle='None',
             marker='o', color='r', alpha=.5)
ax3.set_ylabel('$B^{*}_{2,m}$', fontsize=Title_font_size, color='r')
ax3.yaxis.set_tick_params(labelsize=6)
# tau 2
ax4.errorbar(TIME_BGC, TAU2_ARRAY, yerr=(TAU2_ARRAY - TAU2_ARRAY_LB, TAU2_ARRAY_UB - TAU2_ARRAY), ms=2,
             linestyle='None',
             marker='o', color='g', alpha=.5)
ax4.set_ylabel('${\\tau}_2$', fontsize=Title_font_size, color='g')
ax4.yaxis.set_tick_params(labelsize=6)
ax4.set_ylim(0, 20)
# Sigma
ax5.errorbar(TIME_BGC, SIG2_ARRAY, yerr=(SIG2_ARRAY - SIG2_ARRAY_LB, SIG2_ARRAY_UB - SIG2_ARRAY), ms=2,
             linestyle='None',
             marker='o', color='c', alpha=.5)
ax5.set_ylabel('$\sigma$', fontsize=Title_font_size, color='c')
ax5.yaxis.set_tick_params(labelsize=6)
ax5.set_ylim(0, 5)
# chl-bbp1
ax6.errorbar(TIME_BGC, bbpS1_ARRAY, yerr=(bbpS1_ARRAY - bbpS1_ARRAY_LB, bbpS1_ARRAY_UB - bbpS1_ARRAY), ms=2,
             linestyle='None',
             marker='o', color='c', label='$b^B_{bp,1}$', alpha=.5)
ax6.errorbar(TIME_BGC, bbpS2_ARRAY, yerr=(bbpS2_ARRAY - bbpS2_ARRAY_LB, bbpS2_ARRAY_UB - bbpS2_ARRAY), ms=2,
             linestyle='None',
             marker='o', color='r', label='$b^B_{bp,2}$', alpha=.5)
ax6.set_ylabel('$b^B_{bp}$ (m$^2$ [mg $B$]$^{-1}$)', fontsize=Title_font_size, color='k')
ax6.yaxis.set_tick_params(labelsize=6)
ax6.set_yscale('log')
ax6.legend(loc="upper right", fontsize=6)
ax6.set_ylim(0.00001, 0.2)
# chl-bbp1#
ax7.errorbar(TIME_BGC, bbpk_ARRAY, yerr=(bbpk_ARRAY - bbpk_ARRAY_LB, bbpk_ARRAY_UB - bbpk_ARRAY), ms=2,
             linestyle='None',
             marker='o', color='k', alpha=.2)
ax7.set_ylabel('$b^k_{bp}$ (m$^{-1}$)', fontsize=Title_font_size, color='k')
ax7.yaxis.set_tick_params(labelsize=6)
ax7.xaxis.set_tick_params(labelsize=6)
plt.show()

### STEP 15: PRINT AVERAGE VALUES FOR BBPK, bbpB1 and bbpB2 and confidence ###
# bbpk
print(np.nanmean(bbpk_ARRAY))
print(np.nanmedian(bbpk_ARRAY))
# Comm 1
a = bbpS1_ARRAY
confidence = 0.95
ads = np.where(a > 0.)
a = a[ads]
n = len(a)
m, se = np.mean(a), scipy.stats.sem(a)
h = se * scipy.stats.t.ppf((1 + confidence) / 2., n - 1)
print(np.mean(a))
print(h)
# Comm 2
a = bbpS2_ARRAY
confidence = 0.95
ads = np.where(a > 0.)
a = a[ads]
n = len(a)
m, se = np.mean(a), scipy.stats.sem(a)
h = se * scipy.stats.t.ppf((1 + confidence) / 2., n - 1)
print(np.mean(a))
print(h)

### STEP 16: COMPUTE MEAN LIGHT (PAR) IN MIXED LAYER AND BETWEEN Zm and Zp ###
# Compute light in the mixed layer and at DCM
Light_ML = np.empty([len(JULD)]) + nan
Light_BML = np.empty([len(JULD)]) + nan
Light_DCM = np.empty([len(JULD)]) + nan
Light_MLD_EU = np.empty([len(JULD)]) + nan
for i in range(len(JULD)):
    X1 = DAY_PAR[i, :]
    Y1 = DEPTH[i, :]
    valid1 = np.logical_not(ma.getmask(Y1))
    X1 = X1[valid1]
    Y1 = Y1[valid1]
    valid1 = np.logical_not(ma.getmask(X1))
    X1 = X1[valid1]
    Y1 = Y1[valid1]
    asd = np.where(Y1 <= MLD_TEMP_ARRAY[i])
    Light_ML[i] = np.mean(X1[asd])
    asd = np.where((Y1 >= MLD_TEMP_ARRAY[i]) & (Y1 <= EUPHOTIC_D_ALL[i]))
    if len(asd[0]) > 0:
        Light_MLD_EU[i] = np.mean(X1[asd])
    else:
        Light_MLD_EU[i] = nan
    asd = np.argmin(np.abs(Y1 - MLD_TEMP_ARRAY[i]))
    Light_BML[i] = X1[asd]
    asd = np.argmin(np.abs(Y1 - TAU2_ARRAY[i] / (Kd_ALL[i])))
    Light_DCM[i] = X1[asd]

### STEP 17: SUPPLEMENTARY FIGURE 7 ###
fig, (ax1) = plt.subplots(1, 1, figsize=(8, 6), \
                          gridspec_kw={'hspace': 0.1})
fig.patch.set_facecolor('White')
ax1.scatter(Light_ML, bbpS1_ARRAY, linestyle='None',
            marker='o', color='c', alpha=.5, label='Community 1')
ax1.scatter(Light_MLD_EU, bbpS2_ARRAY, linestyle='None',
            marker='o', color='r', alpha=.5, label='Community 2')
ax1.set_ylabel('$b^B_{bp}$ (m$^2$ [mg $B$]$^{-1}$)', fontsize=12, color='k')
ax1.set_xlabel('PAR (mol quanta m$^{-2}$ d$^{-1}$)', fontsize=12, color='k')
ax1.yaxis.set_tick_params(labelsize=12)
ax1.set_xlim(0.05, 100)
ax1.set_ylim(0.000001, 0.2)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.legend(loc="upper left", fontsize=12)
plt.show()

### STEP 18: CORRELATIONS IN LIGHT AND BBP^B ###
# Community 1
XX = np.log10(Light_ML)
YY = np.log10(bbpS1_ARRAY)
asq = np.where(YY > -500.00)
XX = XX[asq]
YY = YY[asq]
asq = np.where(XX > -500.00)
XX = XX[asq]
YY = YY[asq]
STATS_REG_P1 = stats.pearsonr(XX, YY)
R_CORR_ST = ("{0:.2f}".format(STATS_REG_P1[0]))
P_CORR_ST = ("{0:.3f}".format(STATS_REG_P1[1]))
print([R_CORR_ST, P_CORR_ST])
STATS_REG_P1 = stats.spearmanr(XX, YY)
R_CORR_ST = ("{0:.2f}".format(STATS_REG_P1[0]))
P_CORR_ST = ("{0:.3f}".format(STATS_REG_P1[1]))
print([R_CORR_ST, P_CORR_ST])
# Community 2
XX = np.log10(Light_MLD_EU)
YY = np.log10(bbpS2_ARRAY)
asq = np.where(YY > -500.000)
XX = XX[asq]
YY = YY[asq]
asq = np.where(XX > -500.000)
XX = XX[asq]
YY = YY[asq]
STATS_REG_P1 = stats.pearsonr(XX, YY)
R_CORR_ST = ("{0:.2f}".format(STATS_REG_P1[0]))
P_CORR_ST = ("{0:.3f}".format(STATS_REG_P1[1]))
print([R_CORR_ST, P_CORR_ST])
STATS_REG_P1 = stats.spearmanr(XX, YY)
R_CORR_ST = ("{0:.2f}".format(STATS_REG_P1[0]))
P_CORR_ST = ("{0:.3f}".format(STATS_REG_P1[1]))
print([R_CORR_ST, P_CORR_ST])

### STEP 19a: INTERPOLATE ALL BGC-ARGO DATA ON SAME PRESSURE (DEPTH) AXIS FOR CONTOUR PLOTS###
### STEP 19b: PRODUCE MODELLED DATA ON SAME PRESSURE (DEPTH) AXIS ###
# Define arrays for interpolation
TEMP_INT = np.empty([len(TEMP[:, 0]), len(TEMP[0, :])]) + nan
PSAL_INT = np.empty([len(TEMP[:, 0]), len(TEMP[0, :])]) + nan
DENS_INT = np.empty([len(TEMP[:, 0]), len(TEMP[0, :])]) + nan
CHLA_INT = np.empty([len(TEMP[:, 0]), len(TEMP[0, :])]) + nan
BBP_INT = np.empty([len(TEMP[:, 0]), len(TEMP[0, :])]) + nan
DISS_OXY_INT = np.empty([len(TEMP[:, 0]), len(TEMP[0, :])]) + nan
PAR_INT = np.empty([len(TEMP[:, 0]), len(TEMP[0, :])]) + nan
CHLA_INT_SIM = np.empty([len(TEMP[:, 0]), len(TEMP[0, :])]) + nan
CHLA_INT_DCM = np.empty([len(TEMP[:, 0]), len(TEMP[0, :])]) + nan
CHLA_INT_MLD = np.empty([len(TEMP[:, 0]), len(TEMP[0, :])]) + nan
BVF_INT = np.empty([len(TEMP[:, 0]), len(TEMP[0, :])]) + nan
BBP_INT_SIM = np.empty([len(TEMP[:, 0]), len(TEMP[0, :])]) + nan
BBP_INT_DCM = np.empty([len(TEMP[:, 0]), len(TEMP[0, :])]) + nan
BBP_INT_MLD = np.empty([len(TEMP[:, 0]), len(TEMP[0, :])]) + nan
BBP_INT_BBK = np.empty([len(TEMP[:, 0]), len(TEMP[0, :])]) + nan
# Loop through each profile (float cycle)
for i in range(len(JULD)):
    # Temperature data interpolation
    a = TEMP[i, :]
    b = DEPTH[i, :]
    valid1 = np.logical_not(ma.getmask(b))
    valid2 = np.logical_not(ma.getmask(a))
    a = a[valid2]
    c = b[valid2]
    a = a.compressed()
    b = b.compressed()
    if len(a) > 500:
        interpfunc         = interpolate.interp1d(c,a, kind='linear',fill_value="extrapolate")
        xxx                = interpfunc(b)
        TEMP_INT[i,valid1] = xxx
    # Salinity data interpolation
    a = PSAL[i, :]
    b = DEPTH[i, :]
    valid1 = np.logical_not(ma.getmask(b))
    valid2 = np.logical_not(ma.getmask(a))
    a = a[valid2]
    c = b[valid2]
    a = a.compressed()
    b = b.compressed()
    if len(a) > 500:
        interpfunc = interpolate.interp1d(c, a, kind='linear', fill_value="extrapolate")
        xxx = interpfunc(b)
        PSAL_INT[i, valid1] = xxx
    # Density data interpolation
    a = DENSITY[i, :]
    b = DEPTH[i, :]
    valid1 = np.logical_not(ma.getmask(b))
    valid2 = np.logical_not(ma.getmask(a))
    a = a[valid2]
    c = b[valid2]
    a = a.compressed()
    b = b.compressed()
    if len(a) > 500:
        interpfunc = interpolate.interp1d(c, a, kind='linear', fill_value="extrapolate")
        xxx = interpfunc(b)
        DENS_INT[i, valid1] = xxx
    # Dissolved Oxygen data interpolation
    a = DISS_OXY[i, :]
    b = DEPTH[i, :]
    valid1 = np.logical_not(ma.getmask(b))
    valid2 = np.logical_not(ma.getmask(a))
    a = a[valid2]
    c = b[valid2]
    a = a.compressed()
    b = b.compressed()
    if len(a) > 10:
        interpfunc = interpolate.interp1d(c, a, kind='linear', fill_value="extrapolate")
        xxx = interpfunc(b)
        DISS_OXY_INT[i, valid1] = xxx
    # Chlorophyll-a data interpolation
    a = CHLA[i, :]
    b = DEPTH[i, :]
    valid1 = np.logical_not(ma.getmask(b))
    valid2 = np.logical_not(ma.getmask(a))
    a = a[valid2]
    c = b[valid2]
    a = a.compressed()
    b = b.compressed()
    if len(a) > 10:
        interpfunc = interpolate.interp1d(c, a, kind='linear', fill_value="extrapolate")
        xxx = interpfunc(b)
        CHLA_INT[i, valid1] = xxx
    # Backscattering data interpolation
    a = BBP[i, :]
    b = DEPTH[i, :]
    valid1 = np.logical_not(ma.getmask(b))
    valid2 = np.logical_not(ma.getmask(a))
    a = a[valid2]
    c = b[valid2]
    a = a.compressed()
    b = b.compressed()
    if len(a) > 10:
        interpfunc = interpolate.interp1d(c, a, kind='linear', fill_value="extrapolate")
        xxx = interpfunc(b)
        # BBP_INT[i,valid1] = xxx
        # Median filter the backscattering data for profiles with more than 11 data points
        if len(xxx) > 11:
            BBP_INT[i, valid1] = scipy.signal.medfilt(xxx, kernel_size=11)
        if len(xxx) < 11:
            BBP_INT[i, valid1] = xxx
    # PAR interpolation
    a = PAR[i, :]
    b = DEPTH[i, :]
    valid1 = np.logical_not(ma.getmask(b))
    valid2 = np.logical_not(ma.getmask(a))
    a = a[valid2]
    c = b[valid2]
    a = a.compressed()
    b = b.compressed()
    if len(a) > 10:
        interpfunc = interpolate.interp1d(c, a, kind='linear', fill_value="extrapolate")
        xxx = interpfunc(b)
        PAR_INT[i, valid1] = xxx
    # BVF Interpolation
    a = BRUNT[i, :]
    b = DEPTH[i, :]
    valid1 = np.logical_not(ma.getmask(b))
    valid2 = np.logical_not(ma.getmask(a))
    a = a[valid2]
    c = b[valid2]
    a = a.compressed()
    b = b.compressed()
    if len(a) > 500:
        interpfunc = interpolate.interp1d(c, a, kind='linear', fill_value="extrapolate")
        xxx = interpfunc(b)
        BVF_INT[i, valid1] = xxx
    # Modelled Chlorophyll-a data
    Y2 = DEPTH[i, :] * Kd_ALL[i]
    if np.isnan(BM2_ARRAY[i]):
        MLD_pop = (1 - 1. / (1 + np.exp(-(S1_ARRAY[i]) * (Y2 - TAU1_ARRAY[i]))))
        DCM_pop = MLD_pop * 0
        bbk_pop = MLD_pop * 0 + bbpSk_ARRAY[i]
        CHLA_INT_SIM[i, :] = (MLD_pop + DCM_pop) * CHL_SAT_ALL[i]
        CHLA_INT_MLD[i, :] = MLD_pop * CHL_SAT_ALL[i]
        CHLA_INT_DCM[i, :] = DCM_pop * CHL_SAT_ALL[i]
        BBP_INT_MLD[i, :] = MLD_pop * bbpST1_ARRAY[i] * BBP_SAT_ALL[i]
        BBP_INT_DCM[i, :] = MLD_pop * 0 * BBP_SAT_ALL[i]
        BBP_INT_BBK[i, :] = bbk_pop * BBP_SAT_ALL[i]
        BBP_INT_SIM[i, :] = (MLD_pop * bbpST1_ARRAY[i] + bbpSk_ARRAY[i]) * BBP_SAT_ALL[i]
    else:
        MLD_pop = (1 - 1. / (1 + np.exp(-(S1_ARRAY[i]) * (Y2 - TAU1_ARRAY[i]))))
        DCM_pop = BM2_ARRAY[i] * np.exp(-((Y2 - TAU2_ARRAY[i]) / SIG2_ARRAY[i]) ** 2.)
        bbk_pop = MLD_pop * 0 + bbpSk_ARRAY[i]
        CHLA_INT_SIM[i, :] = (MLD_pop + DCM_pop) * CHL_SAT_ALL[i]
        CHLA_INT_MLD[i, :] = MLD_pop * CHL_SAT_ALL[i]
        CHLA_INT_DCM[i, :] = DCM_pop * CHL_SAT_ALL[i]
        BBP_INT_MLD[i, :] = MLD_pop * bbpST1_ARRAY[i] * BBP_SAT_ALL[i]
        BBP_INT_DCM[i, :] = DCM_pop * bbpST2_ARRAY[i] * BBP_SAT_ALL[i]
        BBP_INT_BBK[i, :] = bbk_pop * BBP_SAT_ALL[i]
        BBP_INT_SIM[i, :] = (MLD_pop * bbpST1_ARRAY[i] + DCM_pop * bbpST2_ARRAY[i] + bbpSk_ARRAY[i]) * BBP_SAT_ALL[i]

### STEP 20: SUPPLEMENTARY FIGURE 5 ###
# Figure parameters that can be changed
TEMP_COL          = mpl.cm.magma  #Temp colour scale (see https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html)
PSAL_COL          = mpl.cm.winter #Salinity colour scale (see https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html)
DOXY_COL          = mpl.cm.copper #Diss OXY colour scale (see https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html)
CHLA_COL          = mpl.cm.viridis#Chl-a colour scale (see https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html)
BBP_COL           = mpl.cm.cividis#bbp colour scale (see https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html)
PAR_COL           = mpl.cm.inferno
DENS_COL          = mpl.cm.hot
BVF_COL           = mpl.cm.gist_heat
XSIZE             = 7            #Define the xsize of the figure window
YSIZE             = 8            #Define the ysize of the figure window
Title_font_size   = 6            #Define the font size of the titles
Label_font_size_x = 6            #Define the font size of the x-labels
Label_font_size_y = 6            #Define the font size of the y-labels
Cbar_title_size   = 7            #Define the font size of the Colourbar title
Cbar_label_size   = 4            #Define the font size of the Colourbar labels
Percentiles_upper = 99           #Upper percentiles used to constrain the colour scale
Percentiles_lower = 1            #Upper percentiles used to constrain the colour scale
#Define the figure window including 5 subplots orientated vertically
fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(8, sharex=True, figsize=(XSIZE,YSIZE), \
    gridspec_kw={'hspace': 0.4})
fig.patch.set_facecolor('White')
#SUBPLOT 1: TEMPERATURE TIME-SERIES
##Constrain data to be between 1 and 99 percentile (avoids outliers in data colour scaling)
IN_DATA        = TEMP_INT
PCT_1          = np.nanpercentile(IN_DATA, Percentiles_lower)
PCT_2          = np.nanpercentile(IN_DATA, Percentiles_upper)
print(PCT_1)
print(PCT_2)
valid          = (IN_DATA < PCT_1)
IN_DATA[valid] = PCT_1
valid          = (IN_DATA > PCT_2)
IN_DATA[valid] = PCT_2
##Define colour levels
levels         = np.arange(PCT_1, PCT_2+((PCT_2-PCT_1)/49.), (PCT_2-PCT_1)/50.)
im1            = ax1.contourf(TIME_BGC_MATRIX, DEPTH, IN_DATA, levels,cmap = TEMP_COL)
##Set axis info and titles
ax1.set_ylim([200,0])
ax1.set_ylabel('Depth (m)', fontsize= Title_font_size, color='k')
ax1.set_title('Temperature ($^o$C)', fontsize = Title_font_size, color='k')
ax1.yaxis.set_tick_params(labelsize= Label_font_size_y)
##Add colourbar
cbar1 = fig.colorbar(im1, ax=ax1)
cbar1.ax.locator_params()
cbar1.set_label("$^o$C", size  = Cbar_title_size)
cbar1.ax.tick_params(labelsize = Cbar_label_size)
#SUBPLOT 2: SALINITY TIME-SERIES
##Constrain data to be between 1 and 99 percentile (avoids outliers in data colour scaling)
IN_DATA        = PSAL_INT
PCT_1          = np.nanpercentile(IN_DATA, Percentiles_lower)
PCT_2          = np.nanpercentile(IN_DATA, Percentiles_upper)
valid          = (IN_DATA < PCT_1)
IN_DATA[valid] = PCT_1
valid          = (IN_DATA > PCT_2)
IN_DATA[valid] = PCT_2
##Define colour levels
levels         = np.arange(PCT_1, PCT_2+((PCT_2-PCT_1)/49.), (PCT_2-PCT_1)/50.)
im2 = ax2.contourf(TIME_BGC_MATRIX, DEPTH, IN_DATA, levels, cmap = PSAL_COL)
##Set axis info and titles
ax2.set_ylim([200,0])
ax2.set_ylabel('Depth (m)', fontsize= Title_font_size, color='k')
ax2.set_title('Salinity (PSU)', fontsize= Title_font_size, color='k')
ax2.yaxis.set_tick_params(labelsize= Label_font_size_y)
##Add colourbar
cbar2 = fig.colorbar(im2, ax=ax2)
cbar2.ax.locator_params()
cbar2.set_label("PSU", size    = Cbar_title_size)
cbar2.ax.tick_params(labelsize = Cbar_label_size)
#SUBPLOT 3: DENSITY TIME-SERIES
##Constrain data to be between 1 and 99 percentile (avoids outliers in data colour scaling)
IN_DATA        = DENS_INT
PCT_1          = np.nanpercentile(IN_DATA, Percentiles_lower)
PCT_2          = np.nanpercentile(IN_DATA, Percentiles_upper)
valid          = (IN_DATA < PCT_1)
IN_DATA[valid] = PCT_1
valid          = (IN_DATA > PCT_2)
IN_DATA[valid] = PCT_2
##Define colour levels
levels         = np.arange(PCT_1, PCT_2+((PCT_2-PCT_1)/49.), (PCT_2-PCT_1)/50.)
im3 = ax3.contourf(TIME_BGC_MATRIX, DEPTH, IN_DATA, levels, cmap = DENS_COL)
##Set axis info and titles
ax3.set_ylim([200,0])
ax3.set_ylabel('Depth (m)', fontsize= Title_font_size, color='k')
ax3.set_title('Density (kg m$^{-3}$)', fontsize= Title_font_size, color='k')
ax3.yaxis.set_tick_params(labelsize= Label_font_size_y)
##Add colourbar
cbar3 = fig.colorbar(im3, ax=ax3)
cbar3.ax.locator_params()
cbar3.set_label("kg m$^{-3}$", size    = Cbar_title_size)
cbar3.ax.tick_params(labelsize = Cbar_label_size)
#SUBPLOT 4: PAR TIME-SERIES
##Constrain data to be between 1 and 99 percentile (avoids outliers in data colour scaling)
IN_DATA        = np.log(PAR_INT)
PCT_1          = np.nanpercentile(IN_DATA, Percentiles_lower)
PCT_2          = np.nanpercentile(IN_DATA, Percentiles_upper)
valid          = (IN_DATA < PCT_1)
IN_DATA[valid] = PCT_1
valid          = (IN_DATA > PCT_2)
IN_DATA[valid] = PCT_2
##Define colour levels
levels         = np.arange(PCT_1, PCT_2+((PCT_2-PCT_1)/49.), (PCT_2-PCT_1)/50.)
im4 = ax4.contourf(TIME_BGC_MATRIX, DEPTH, IN_DATA, levels, cmap = PAR_COL)
##Set axis info and titles
ax4.set_ylim([200,0])
ax4.set_ylabel('Depth (m)', fontsize= Title_font_size, color='k')
ax4.set_title('log(PAR)', fontsize= Title_font_size, color='k')
ax4.yaxis.set_tick_params(labelsize= Label_font_size_y)
##Add colourbar
cbar4 = fig.colorbar(im4, ax=ax4)
cbar4.ax.locator_params()
cbar4.set_label("log(PAR)", size    = Cbar_title_size)
cbar4.ax.tick_params(labelsize = Cbar_label_size)
#SUBPLOT 5: DISSOLVED OXYGEN TIME-SERIES
#Constrain data to be between 1 and 99 percentile (avoids outliers in data colour scaling)
IN_DATA        = DISS_OXY_INT
PCT_1          = np.nanpercentile(IN_DATA, Percentiles_lower)
PCT_2          = np.nanpercentile(IN_DATA, Percentiles_upper)
valid          = (IN_DATA < PCT_1)
IN_DATA[valid] = PCT_1
valid          = (IN_DATA > PCT_2)
IN_DATA[valid] = PCT_2
##Define colour levels
levels         = np.arange(PCT_1, PCT_2+((PCT_2-PCT_1)/49.), (PCT_2-PCT_1)/50.)
im5 = ax5.contourf(TIME_BGC_MATRIX, DEPTH, IN_DATA, levels,cmap = DOXY_COL)
##Set axis info and titles
ax5.set_ylim([200,0])
ax5.set_ylabel('Depth (m)', fontsize= Title_font_size, color='k')
ax5.set_title('Dissolved Oxygen (micro mol kg$^{-3}$)', fontsize= Title_font_size, color='k')
ax5.yaxis.set_tick_params(labelsize= Label_font_size_y)
##Add colourbar
cbar5 = fig.colorbar(im5, ax=ax5)
cbar5.ax.locator_params()
cbar5.set_label("micro mol kg$^{-3}$", size= Cbar_title_size)
cbar5.ax.tick_params(labelsize= Cbar_label_size)
#SUBPLOT 6: CHLOROPHYLL-A TIME-SERIES
##Constrain data to be between 1 and 99 percentile (avoids outliers in data colour scaling)
IN_DATA        = CHLA_INT
PCT_1          = np.nanpercentile(IN_DATA, Percentiles_lower)
PCT_2          = np.nanpercentile(IN_DATA, Percentiles_upper)
valid          = (IN_DATA < PCT_1)
IN_DATA[valid] = PCT_1
valid          = (IN_DATA > PCT_2)
IN_DATA[valid] = PCT_2
##Define colour levels
levels         = np.arange(PCT_1, PCT_2+((PCT_2-PCT_1)/49.), (PCT_2-PCT_1)/50.)
im6 = ax6.contourf(TIME_BGC_MATRIX, DEPTH, IN_DATA, levels, cmap = CHLA_COL)
##Set axis info and titles
ax6.set_ylim([200,0])
ax6.set_ylabel('Depth (m)', fontsize= Title_font_size, color='k')
ax6.set_title('Chl-a (mg m$^{-3}$)', fontsize= Title_font_size, color='k')
ax6.yaxis.set_tick_params(labelsize= Label_font_size_y)
##Add colourbar
cbar6 = fig.colorbar(im6, ax=ax6)
cbar6.ax.locator_params()
cbar6.set_label("mg m$^{-3}$", size = Cbar_title_size)
cbar6.ax.tick_params(labelsize = Cbar_label_size)
#SUBPLOT 7: BBP TIME-SERIES
##Constrain data to be between 1 and 99 percentile (avoids outliers in data colour scaling)
IN_DATA        = BBP_INT
PCT_1          = np.nanpercentile(IN_DATA, Percentiles_lower)
PCT_2          = np.nanpercentile(IN_DATA, Percentiles_upper)
valid          = (IN_DATA < PCT_1)
IN_DATA[valid] = PCT_1
valid          = (IN_DATA > PCT_2)
IN_DATA[valid] = PCT_2
##Define colour levels
levels         = np.arange(PCT_1, PCT_2+((PCT_2-PCT_1)/49.), (PCT_2-PCT_1)/50.)
im7 = ax7.contourf(TIME_BGC_MATRIX, DEPTH, IN_DATA, levels,cmap = BBP_COL)
##Set axis info and titles
ax7.set_ylim([200,0])
ax7.set_ylabel('Depth (m)', fontsize= Title_font_size, color='k')
ax7.set_title('$b_{bp}$ (m$^{-1}$)', fontsize= Title_font_size, color='k')
ax7.yaxis.set_tick_params(labelsize= Label_font_size_y)
ax7.xaxis.set_tick_params(labelsize= Label_font_size_x)
##Add colourbar
cbar7 = fig.colorbar(im7, ax=ax7)
cbar7.ax.locator_params()
cbar7.set_label("m$^{-1}$", size = Cbar_title_size)
cbar7.ax.tick_params(labelsize = Cbar_label_size)
#SUBPLOT 8: BVFTIME-SERIES
#Constrain data to be between 1 and 99 percentile (avoids outliers in data colour scaling)
IN_DATA        = BVF_INT
PCT_1          = np.nanpercentile(IN_DATA, Percentiles_lower)
PCT_2          = np.nanpercentile(IN_DATA, Percentiles_upper)
valid          = (IN_DATA < PCT_1)
IN_DATA[valid] = PCT_1
valid          = (IN_DATA > PCT_2)
IN_DATA[valid] = PCT_2
##Define colour levels
levels         = np.arange(PCT_1, PCT_2+((PCT_2-PCT_1)/49.), (PCT_2-PCT_1)/50.)
im8 = ax8.contourf(TIME_BGC_MATRIX, DEPTH, IN_DATA, levels,cmap = BVF_COL)
##Set axis info and titles
ax8.set_ylim([200,0])
ax8.set_ylabel('Depth (m)', fontsize= Title_font_size, color='k')
ax8.set_title('Brunt-Väisälä Frequency', fontsize= Title_font_size, color='k')
ax8.yaxis.set_tick_params(labelsize= Label_font_size_y)
ax8.xaxis.set_tick_params(labelsize= Label_font_size_y)
ax8.set_xlabel('Time', fontsize=8, color='k')
##Add colourbar
cbar8 = fig.colorbar(im8, ax=ax8)
cbar8.ax.locator_params()
cbar8.set_label(" ", size= Cbar_title_size)
cbar8.ax.tick_params(labelsize= Cbar_label_size)
plt.show()

### STEP 21: FIGURE 3 ###
#Figure parameters that can be changed
CHLA_COL          = mpl.cm.viridis#Chl-a colour scale (see https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html)
DIFF_COL          = mpl.cm.bwr    #bbp colour scale (see https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html)
XSIZE             = 7            #Define the xsize of the figure window
YSIZE             = 7            #Define the ysize of the figure window
Title_font_size   = 7            #Define the font size of the titles
Label_font_size_x = 5            #Define the font size of the x-labels
Label_font_size_y = 6            #Define the font size of the y-labels
Cbar_title_size   = 7            #Define the font size of the Colourbar title
Cbar_label_size   = 6            #Define the font size of the Colourbar labels
#Define the figure window including 5 subplots orientated vertically
fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, sharex=True, figsize=(XSIZE,YSIZE), \
    gridspec_kw={'hspace': 0.3})
fig.patch.set_facecolor('White')
#SUBPLOT 1: CHLOROPHYLL-A TIME-SERIES
##Set any CHl-a data below 0.001 to 0.001 (minimim detectable limit)
IN_DATA        = CHLA_INT
PCT_1          = 0.0
PCT_2          = 0.3
valid          = (IN_DATA < PCT_1)
IN_DATA[valid] = PCT_1
valid          = (IN_DATA > PCT_2)
IN_DATA[valid] = PCT_2
##Define colour levels
levels       = np.linspace(PCT_1, PCT_2, 50)
im1 = ax1.contourf(TIME_BGC_MATRIX, PRESSURE, IN_DATA, levels, cmap = CHLA_COL)
##Set axis info and titles
ax1.set_ylim([200,0])
ax1.set_ylabel('Depth (m)', fontsize= Title_font_size, color='k')
ax1.set_title('(a) Data', fontsize= Title_font_size, color='k')
ax1.yaxis.set_tick_params(labelsize= Label_font_size_y)
##Add colourbar
cbar1 = fig.colorbar(im1, ax=ax1)
cbar1.ax.locator_params()
cbar1.set_label("mg m$^{-3}$", size = Cbar_title_size)
cbar1.ax.tick_params(labelsize = Cbar_label_size)
cbar1.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.3f'))
#SUBPLOT 2: CHLOROPHYLL-A TIME-SERIES
##Constrain data to be between 1 and 99 percentile (avoids outliers in data colour scaling)
IN_DATA        = CHLA_INT_SIM
PCT_1          = 0.0
PCT_2          = 0.3
valid          = (IN_DATA < PCT_1)
IN_DATA[valid] = PCT_1
valid          = (IN_DATA > PCT_2)
IN_DATA[valid] = PCT_2
##Define colour levels
levels       = np.linspace(PCT_1, PCT_2, 50)
im2          = ax2.contourf(TIME_BGC_MATRIX, PRESSURE, IN_DATA, levels, cmap = CHLA_COL)
##Set axis info and titles
ax2.set_ylim([200,0])
ax2.set_ylabel('Depth (m)', fontsize= Title_font_size, color='k')
ax2.set_title('(b) Model', fontsize= Title_font_size, color='k')
ax2.yaxis.set_tick_params(labelsize= Label_font_size_y)
##Add colourbar
cbar2 = fig.colorbar(im2, ax=ax2)
cbar2.ax.locator_params()
cbar2.set_label("mg m$^{-3}$", size = Cbar_title_size)
cbar2.ax.tick_params(labelsize = Cbar_label_size)
cbar2.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.3f'))
#SUBPLOT 3: Difference
##Constrain data to be between 1 and 99 percentile (avoids outliers in data colour scaling)
IN_DATA        = CHLA_INT_SIM - CHLA_INT
PCT_1          = -0.3
PCT_2          = 0.3
valid          = (IN_DATA < PCT_1)
IN_DATA[valid] = PCT_1
valid          = (IN_DATA > PCT_2)
IN_DATA[valid] = PCT_2
##Define colour levels
levels       = np.linspace(PCT_1, PCT_2, 50)
im3          = ax3.contourf(TIME_BGC_MATRIX, PRESSURE, IN_DATA, levels, cmap = DIFF_COL)
##Set axis info and titles
ax3.set_ylim([200,0])
ax3.set_ylabel('Depth (m)', fontsize= Title_font_size, color='k')
ax3.set_title('(c) Model - Data', fontsize= Title_font_size, color='k')
ax3.yaxis.set_tick_params(labelsize= Label_font_size_y)
##Add colourbar
cbar3 = fig.colorbar(im3, ax=ax3)
cbar3.ax.locator_params()
cbar3.set_label("mg m$^{-3}$", size = Cbar_title_size)
cbar3.ax.tick_params(labelsize = Cbar_label_size)
cbar3.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.3f'))
#SUBPLOT 4: CHLOROPHYLL-A TIME-SERIES Surface
##Constrain data to be between 1 and 99 percentile (avoids outliers in data colour scaling)
IN_DATA        = CHLA_INT_MLD
PCT_1          = 0.0
PCT_2          = 0.3
valid          = (IN_DATA < PCT_1)
IN_DATA[valid] = PCT_1
valid          = (IN_DATA > PCT_2)
IN_DATA[valid] = PCT_2
##Define colour levels
levels       = np.linspace(PCT_1, PCT_2, 50)
im4          = ax4.contourf(TIME_BGC_MATRIX, PRESSURE, IN_DATA, levels, cmap = CHLA_COL)
##Set axis info and titles
ax4.set_ylim([200,0])
ax4.set_ylabel('Depth (m)', fontsize= Title_font_size, color='k')
ax4.set_title('(d) Community 1', fontsize= Title_font_size, color='k')
ax4.yaxis.set_tick_params(labelsize= Label_font_size_y)
##Add colourbar
cbar4 = fig.colorbar(im4, ax=ax4)
cbar4.ax.locator_params()
cbar4.set_label("mg m$^{-3}$", size = Cbar_title_size)
cbar4.ax.tick_params(labelsize = Cbar_label_size)
cbar4.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.3f'))
#SUBPLOT 5: CHLOROPHYLL-A TIME-SERIES Sub-Surface
##Constrain data to be between 1 and 99 percentile (avoids outliers in data colour scaling)
IN_DATA        = CHLA_INT_DCM
PCT_1          = 0.0
PCT_2          = 0.3
valid          = (IN_DATA < PCT_1)
IN_DATA[valid] = PCT_1
valid          = (IN_DATA > PCT_2)
IN_DATA[valid] = PCT_2
##Define colour levels
levels       = np.linspace(PCT_1, PCT_2, 50)
im5          = ax5.contourf(TIME_BGC_MATRIX, PRESSURE, IN_DATA, levels, cmap = CHLA_COL)
##Set axis info and titles
ax5.set_ylim([200,0])
ax5.set_ylabel('Depth (m)', fontsize= Title_font_size, color='k')
ax5.set_title('(e) Community 2', fontsize= Title_font_size, color='k')
ax5.yaxis.set_tick_params(labelsize= Label_font_size_y)
ax5.xaxis.set_tick_params(labelsize= Label_font_size_y)
ax5.set_xlabel('Time', fontsize=Title_font_size, color='k')
##Add colourbar
cbar5 = fig.colorbar(im5, ax=ax5)
cbar5.ax.locator_params()
cbar5.set_label("mg m$^{-3}$", size = Cbar_title_size)
cbar5.ax.tick_params(labelsize = Cbar_label_size)
cbar5.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.3f'))
plt.show()

### STEP 22: FIGURE 5 ###
#Figure parameters that can be changed
BBP_COL           = mpl.cm.cividis#bbp colour scale (see https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html)
DIFF_COL          = mpl.cm.bwr    #bbp colour scale (see https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html)
XSIZE             = 7            #Define the xsize of the figure window
YSIZE             = 8            #Define the ysize of the figure window
Title_font_size   = 7            #Define the font size of the titles
Label_font_size_x = 5            #Define the font size of the x-labels
Label_font_size_y = 6            #Define the font size of the y-labels
Cbar_title_size   = 7            #Define the font size of the Colourbar title
Cbar_label_size   = 6            #Define the font size of the Colourbar labels
#Define the figure window including 5 subplots orientated vertically
fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, sharex=True, figsize=(XSIZE,YSIZE), \
    gridspec_kw={'hspace': 0.25})
fig.patch.set_facecolor('White')
#SUBPLOT 1: BBP TIME-SERIES
##Set any BBP data below 0.001 to 0.001 (minimum detectable limit)
IN_DATA        = BBP_INT
PCT_1          = 0.0
PCT_2          = 0.001
valid          = (IN_DATA < PCT_1)
IN_DATA[valid] = PCT_1
valid          = (IN_DATA > PCT_2)
IN_DATA[valid] = PCT_2
##Define colour levels
levels       = np.linspace(PCT_1, PCT_2, 50)
im1 = ax1.contourf(TIME_BGC_MATRIX, PRESSURE, IN_DATA, levels, cmap = BBP_COL)
##Set axis info and titles
ax1.set_ylim([200,0])
ax1.set_ylabel('Depth (m)', fontsize= Title_font_size, color='k')
ax1.set_title('(a) Data', fontsize= Title_font_size, color='k')
ax1.yaxis.set_tick_params(labelsize= Label_font_size_y)
##Add colourbar
cbar1 = fig.colorbar(im1, ax=ax1)
cbar1.ax.locator_params()
cbar1.set_label("m$^{-1}$", size = Cbar_title_size)
cbar1.ax.tick_params(labelsize = Cbar_label_size)
cbar1.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.4f'))
#SUBPLOT 2: BBP TIME-SERIES
##Constrain data to be between 1 and 99 percentile (avoids outliers in data colour scaling)
IN_DATA        = BBP_INT_SIM
PCT_1          = 0.0
PCT_2          = 0.001
valid          = (IN_DATA < PCT_1)
IN_DATA[valid] = PCT_1
valid          = (IN_DATA > PCT_2)
IN_DATA[valid] = PCT_2
##Define colour levels
levels       = np.linspace(PCT_1, PCT_2, 50)
im2          = ax2.contourf(TIME_BGC_MATRIX, PRESSURE, IN_DATA, levels, cmap = BBP_COL)
##Set axis info and titles
ax2.set_ylim([200,0])
ax2.set_ylabel('Depth (m)', fontsize= Title_font_size, color='k')
ax2.set_title('(b) Model', fontsize= Title_font_size, color='k')
ax2.yaxis.set_tick_params(labelsize= Label_font_size_y)
##Add colourbar
cbar2 = fig.colorbar(im2, ax=ax2)
cbar2.ax.locator_params()
cbar2.set_label("m$^{-1}$", size = Cbar_title_size)
cbar2.ax.tick_params(labelsize = Cbar_label_size)
cbar2.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.4f'))
#SUBPLOT 3: Difference
##Constrain data to be between 1 and 99 percentile (avoids outliers in data colour scaling)
IN_DATA        = BBP_INT_SIM - BBP_INT
PCT_1          = -0.001
PCT_2          = 0.001
valid          = (IN_DATA < PCT_1)
IN_DATA[valid] = PCT_1
valid          = (IN_DATA > PCT_2)
IN_DATA[valid] = PCT_2
##Define colour levels
levels       = np.linspace(PCT_1, PCT_2, 50)
im3          = ax3.contourf(TIME_BGC_MATRIX, PRESSURE, IN_DATA, levels, cmap = DIFF_COL)
##Set axis info and titles
ax3.set_ylim([200,0])
ax3.set_ylabel('Depth (m)', fontsize= Title_font_size, color='k')
ax3.set_title('(c) Model - Data', fontsize= Title_font_size, color='k')
ax3.yaxis.set_tick_params(labelsize= Label_font_size_y)
##Add colourbar
cbar3 = fig.colorbar(im3, ax=ax3)
cbar3.ax.locator_params()
cbar3.set_label("m$^{-1}$", size = Cbar_title_size)
cbar3.ax.tick_params(labelsize = Cbar_label_size)
cbar3.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.4f'))
#SUBPLOT 4: BBP TIME-SERIES Surface
##Constrain data to be between 1 and 99 percentile (avoids outliers in data colour scaling)
IN_DATA        = BBP_INT_MLD
PCT_1          = 0.0
PCT_2          = 0.001
valid          = (IN_DATA < PCT_1)
IN_DATA[valid] = PCT_1
valid          = (IN_DATA > PCT_2)
IN_DATA[valid] = PCT_2
##Define colour levels
levels       = np.linspace(PCT_1, PCT_2, 50)
im4          = ax4.contourf(TIME_BGC_MATRIX, PRESSURE, IN_DATA, levels, cmap = BBP_COL)
##Set axis info and titles
ax4.set_ylim([200,0])
ax4.set_ylabel('Depth (m)', fontsize= Title_font_size, color='k')
ax4.set_title('(d) Community 1', fontsize= Title_font_size, color='k')
ax4.yaxis.set_tick_params(labelsize= Label_font_size_y)
##Add colourbar
cbar4 = fig.colorbar(im4, ax=ax4)
cbar4.ax.locator_params()
cbar4.set_label("m$^{-1}$", size = Cbar_title_size)
cbar4.ax.tick_params(labelsize = Cbar_label_size)
cbar4.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.4f'))
#SUBPLOT 5: BBP TIME-SERIES Sub-Surface
##Constrain data to be between 1 and 99 percentile (avoids outliers in data colour scaling)
IN_DATA        = BBP_INT_DCM
PCT_1          = 0.0
PCT_2          = 0.001
valid          = (IN_DATA < PCT_1)
IN_DATA[valid] = PCT_1
valid          = (IN_DATA > PCT_2)
IN_DATA[valid] = PCT_2
##Define colour levels
levels       = np.linspace(PCT_1, PCT_2, 50)
im5          = ax5.contourf(TIME_BGC_MATRIX, PRESSURE, IN_DATA, levels, cmap = BBP_COL)
##Set axis info and titles
ax5.set_ylim([200,0])
ax5.set_ylabel('Depth (m)', fontsize= Title_font_size, color='k')
ax5.set_title('(e) Community 2', fontsize= Title_font_size, color='k')
ax5.yaxis.set_tick_params(labelsize= Label_font_size_y)
#ax5.set_xlabel('Time', fontsize=15, color='k')
##Add colourbar
cbar5 = fig.colorbar(im5, ax=ax5)
cbar5.ax.locator_params()
cbar5.set_label("m$^{-1}$", size = Cbar_title_size)
cbar5.ax.tick_params(labelsize = Cbar_label_size)
cbar5.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.4f'))
#SUBPLOT 6: BBPk TIME-SERIES
##Constrain data to be between 1 and 99 percentile (avoids outliers in data colour scaling)
IN_DATA        = BBP_INT_BBK
PCT_1          = 0.0
PCT_2          = 0.001
valid          = (IN_DATA < PCT_1)
IN_DATA[valid] = PCT_1
valid          = (IN_DATA > PCT_2)
IN_DATA[valid] = PCT_2
##Define colour levels
levels       = np.linspace(PCT_1, PCT_2, 50)
im6          = ax6.contourf(TIME_BGC_MATRIX, PRESSURE, IN_DATA, levels, cmap = BBP_COL)
##Set axis info and titles
ax6.set_ylim([200,0])
ax6.set_ylabel('Depth (m)', fontsize= Title_font_size, color='k')
ax6.set_title('(f) $b^k_{bp}$', fontsize= Title_font_size, color='k')
ax6.yaxis.set_tick_params(labelsize= Label_font_size_y)
ax6.xaxis.set_tick_params(labelsize= Label_font_size_y)
ax6.set_xlabel('Time', fontsize=Title_font_size, color='k')
##Add colourbar
cbar6 = fig.colorbar(im6, ax=ax6)
cbar6.ax.locator_params()
cbar6.set_label("m$^{-1}$", size = Cbar_title_size)
cbar6.ax.tick_params(labelsize = Cbar_label_size)
cbar6.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.4f'))
plt.show()

### STEP 23: INTEGRATE CHL-A DATA ###
#Define variables
CHL_INT_DATA    = np.empty([len(JULD)])+nan
CHL_INT_MOD_ALL = np.empty([len(JULD)])+nan
CHL_INT_MOD_AS1 = np.empty([len(JULD)])+nan
CHL_INT_MOD_AS2 = np.empty([len(JULD)])+nan
BVF_INDEX_AV    = np.empty([len(JULD)])+nan
for INDEX_PROFILE in range(len(JULD)):
    DEPTH_PROF = DEPTH[INDEX_PROFILE,:]
    idx       = np.where(DEPTH_PROF <= EUPHOTIC_D_ALL[INDEX_PROFILE]*1.5)
    CHL_INT_DATA[INDEX_PROFILE]     = spi.trapz(CHLA_INT[INDEX_PROFILE,idx],DEPTH[INDEX_PROFILE,idx])
    CHL_INT_MOD_ALL[INDEX_PROFILE]  = spi.trapz(CHLA_INT_SIM[INDEX_PROFILE,idx],DEPTH[INDEX_PROFILE,idx])
    CHL_INT_MOD_AS1[INDEX_PROFILE]  = spi.trapz(CHLA_INT_MLD[INDEX_PROFILE,idx],DEPTH[INDEX_PROFILE,idx])
    CHL_INT_MOD_AS2[INDEX_PROFILE]  = spi.trapz(CHLA_INT_DCM[INDEX_PROFILE,idx],DEPTH[INDEX_PROFILE,idx])
    BVF_INDEX_AV[INDEX_PROFILE]     = np.nanmean(BVF_INT[INDEX_PROFILE,idx])

### STEP 24: COMPUTE PHENOLOGY METRICS ###
# SURFACE COMMUNITY
asd = np.where(CHL_INT_MOD_AS1 >= 0.0)
TIME_BGC_AS1 = TIME_BGC[asd]
JULD_AS1 = JULD[asd]
CHL_INT_MOD_AS1b = CHL_INT_MOD_AS1[asd]
# Define New Time
JULD_INT_AS1 = np.linspace(24014., 24378., 365)
REFERENCE_DATE = 2433282.50000  # JULIAN DATE for reference 1950-01-01 00:00:00
TIME_BGC_INT_AS1 = np.empty(len(JULD_INT_AS1), dtype='datetime64[s]')
for i in range(len(JULD_INT_AS1)):
    TIME_BGC_INT_AS1[i] = julian.from_jd(JULD_INT_AS1[i] + REFERENCE_DATE, fmt='jd')
# Interpolate data to new time
a = JULD_AS1
b = CHL_INT_MOD_AS1b
valid = np.logical_not(ma.getmask(a))
a = a[valid]
b = b[valid]
a = a.compressed()
interpfunc = interpolate.interp1d(a, b, kind='linear', fill_value="extrapolate")
CHL_INT_MOD_AS1_INT = interpfunc(JULD_INT_AS1)
# Smooth the integrated communities
CHL_INT_MOD_AS1_INT = savgol_filter(CHL_INT_MOD_AS1_INT, 15, 2, mode='nearest')
# Compute median + 5%, subtract for anomaly and compute cumulative sum + gradient
PERCENT_INC = 0.05
CHL_INT_MOD_AS1_MED = np.median(CHL_INT_MOD_AS1_INT) + (np.median(CHL_INT_MOD_AS1_INT) * PERCENT_INC)
CHL_INT_MOD_AS1_NOM = CHL_INT_MOD_AS1_INT - CHL_INT_MOD_AS1_MED
CHL_INT_MOD_AS1_CUM = np.cumsum(CHL_INT_MOD_AS1_NOM)
CHL_INT_MOD_AS1_GRAD = np.gradient(CHL_INT_MOD_AS1_CUM)
# Window of above and below zero
WINDOW = 15
# Initiation index
asd = np.where(CHL_INT_MOD_AS1_GRAD > 0)
IND_A = np.min(asd)
IND_B = np.max(asd)
CHL_INT_MOD_AS1_GRAD_1 = CHL_INT_MOD_AS1_GRAD[IND_A:IND_B]
JULD_INT_AS1_2 = JULD_INT_AS1[IND_A:IND_B]
TIME_BGC_INT_AS1_2 = TIME_BGC_INT_AS1[IND_A:IND_B]
INDEX_asd = np.empty([len(CHL_INT_MOD_AS1_GRAD_1)]) * nan
for k in range(len(INDEX_asd) - WINDOW):
    ads_16 = CHL_INT_MOD_AS1_GRAD_1[k:k + (WINDOW - 1)]
    if np.min(ads_16) > 0:
        INDEX_asd[k] = 1
asd2 = np.where(INDEX_asd > 0)
INITIATION_JULD_AS1 = JULD_INT_AS1_2[np.min(asd2)]
INITIATION_TIME_AS1 = TIME_BGC_INT_AS1_2[np.min(asd2)]
# Termination index
asd = np.argmax(CHL_INT_MOD_AS1_GRAD)
CHL_INT_MOD_AS1_GRAD_2 = CHL_INT_MOD_AS1_GRAD[asd:]
JULD_INT_AS1_2 = JULD_INT_AS1[asd:]
TIME_INT_BGC_AS1_2 = TIME_BGC_INT_AS1[asd:]
asd = np.where(CHL_INT_MOD_AS1_GRAD_2 < 0)
IND_A = np.min(asd)
IND_B = np.max(asd)
CHL_INT_MOD_AS1_GRAD_3 = CHL_INT_MOD_AS1_GRAD_2[IND_A:IND_B]
JULD_INT_AS1_3 = JULD_INT_AS1_2[IND_A:IND_B]
TIME_BGC_INT_AS1_3 = TIME_INT_BGC_AS1_2[IND_A:IND_B]
INDEX_asd = np.empty([len(CHL_INT_MOD_AS1_GRAD_3)]) * nan
for k in range(len(INDEX_asd) - WINDOW):
    ads_16 = CHL_INT_MOD_AS1_GRAD_3[k:k + (WINDOW - 1)]
    if np.max(ads_16) < 0:
        INDEX_asd[k] = 1
asd2 = np.where(INDEX_asd > 0)
TERMINATION_JULD_AS1 = JULD_INT_AS1_3[np.min(asd2)]
TERMINATION_TIME_AS1 = TIME_BGC_INT_AS1_3[np.min(asd2)]
# SUB-SURFACE COMMUNITY
asd = np.where(CHL_INT_MOD_AS2 >= 0.0)
TIME_BGC_AS2 = TIME_BGC[asd]
JULD_AS2 = JULD[asd]
CHL_INT_MOD_AS2b = CHL_INT_MOD_AS2[asd]
# Define New Time
JULD_INT_AS2 = np.linspace(24124., 24478., 365)
REFERENCE_DATE = 2433282.50000  # JULIAN DATE for reference 1950-01-01 00:00:00
TIME_BGC_INT_AS2 = np.empty(len(JULD_INT_AS2), dtype='datetime64[s]')
for i in range(len(JULD_INT_AS2)):
    TIME_BGC_INT_AS2[i] = julian.from_jd(JULD_INT_AS2[i] + REFERENCE_DATE, fmt='jd')
# Interpolate data to new time
a = JULD_AS2
b = CHL_INT_MOD_AS2b
valid = np.logical_not(ma.getmask(a))
a = a[valid]
b = b[valid]
a = a.compressed()
interpfunc = interpolate.interp1d(a, b, kind='linear', fill_value="extrapolate")
CHL_INT_MOD_AS2_INT = interpfunc(JULD_INT_AS2)
# Smooth the integrated communities
CHL_INT_MOD_AS2_INT = savgol_filter(CHL_INT_MOD_AS2_INT, 15, 2, mode='nearest')
# Compute median + 5%, subtract for anomaly and compute cumulative sum + gradient
PERCENT_INC = 0.05
CHL_INT_MOD_AS2_MED = np.median(CHL_INT_MOD_AS2_INT) + (np.median(CHL_INT_MOD_AS2_INT) * PERCENT_INC)
CHL_INT_MOD_AS2_NOM = CHL_INT_MOD_AS2_INT - CHL_INT_MOD_AS2_MED
CHL_INT_MOD_AS2_CUM = np.cumsum(CHL_INT_MOD_AS2_NOM)
CHL_INT_MOD_AS2_GRAD = np.gradient(CHL_INT_MOD_AS2_CUM)
# Initiation index
asd = np.where(CHL_INT_MOD_AS2_GRAD > 0)
IND_A = np.min(asd)
IND_B = np.max(asd)
CHL_INT_MOD_AS2_GRAD_1 = CHL_INT_MOD_AS2_GRAD[IND_A:IND_B]
JULD_INT_AS2_2 = JULD_INT_AS2[IND_A:IND_B]
TIME_BGC_INT_AS2_2 = TIME_BGC_INT_AS2[IND_A:IND_B]
INDEX_asd = np.empty([len(CHL_INT_MOD_AS2_GRAD_1)]) * nan
for k in range(len(INDEX_asd) - WINDOW):
    ads_16 = CHL_INT_MOD_AS2_GRAD_1[k:k + (WINDOW - 1)]
    if np.min(ads_16) > 0:
        INDEX_asd[k] = 1
asd2 = np.where(INDEX_asd > 0)
INITIATION_JULD_AS2 = JULD_INT_AS2_2[np.min(asd2)]
INITIATION_TIME_AS2 = TIME_BGC_INT_AS2_2[np.min(asd2)]
# Termination index
asd = np.argmax(CHL_INT_MOD_AS2_GRAD)
CHL_INT_MOD_AS2_GRAD_2 = CHL_INT_MOD_AS2_GRAD[asd:]
JULD_INT_AS2_2 = JULD_INT_AS2[asd:]
TIME_INT_BGC_AS2_2 = TIME_BGC_INT_AS2[asd:]
asd = np.where(CHL_INT_MOD_AS2_GRAD_2 < 0)
IND_A = np.min(asd)
IND_B = np.max(asd)
CHL_INT_MOD_AS2_GRAD_3 = CHL_INT_MOD_AS2_GRAD_2[IND_A:IND_B]
JULD_INT_AS2_3 = JULD_INT_AS2_2[IND_A:IND_B]
TIME_BGC_INT_AS2_3 = TIME_INT_BGC_AS2_2[IND_A:IND_B]
INDEX_asd = np.empty([len(CHL_INT_MOD_AS2_GRAD_3)]) * nan
for k in range(len(INDEX_asd) - WINDOW):
    ads_16 = CHL_INT_MOD_AS2_GRAD_3[k:k + (WINDOW - 1)]
    # print(ads_16)
    if np.max(ads_16) < 0:
        INDEX_asd[k] = 1
asd2 = np.where(INDEX_asd > 0)
TERMINATION_JULD_AS2 = JULD_INT_AS2_3[np.min(asd2)]
TERMINATION_TIME_AS2 = TIME_BGC_INT_AS2_3[np.min(asd2)]

### STEP 24: FIGURE 4 OF THE PAPER ###
#Smooth the integrated communities
CHL_INT_MOD_AS1_SM = savgol_filter(CHL_INT_MOD_AS1, 11, 2, mode='nearest')
CHL_INT_MOD_AS2_SM = savgol_filter(CHL_INT_MOD_AS2, 11, 2, mode='nearest')
#Define the figure window including 5 subplots orientated vertically
XSIZE             = 6            #Define the xsize of the figure window
YSIZE             = 7            #Define the ysize of the figure window
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, figsize=(XSIZE,YSIZE), \
    gridspec_kw={'hspace': 0.25})
fig.patch.set_facecolor('White')
ax1.plot(TIME_BGC,CHL_INT_DATA, label = 'Data', color = 'k')
ax1.plot(TIME_BGC,CHL_INT_MOD_ALL, label = 'Model', color = 'g')
ax1.set_ylabel('Chl-a (mg m$^{-2}$)', fontsize= 7, color='k')
ax1.set_title('(a)', fontsize= 12, color='k')
ax1.legend(loc="upper left", fontsize=7)
ax1.yaxis.set_tick_params(labelsize= 6)
ax2.fill_between([INITIATION_TIME_AS1, TERMINATION_TIME_AS1],[40,40],facecolor='red', alpha=0.1)
ax2.fill_between([INITIATION_TIME_AS2, TERMINATION_TIME_AS2],[40,40],facecolor='blue', alpha=0.1)
ax2.plot(TIME_BGC,CHL_INT_MOD_AS1, label = 'Community 1', color = 'pink')
ax2.plot(TIME_BGC,CHL_INT_MOD_AS1_SM, label = 'Community 1 Smooth', color = 'r')
ax2.plot(TIME_BGC,CHL_INT_MOD_AS2, label = 'Community 2', color = 'skyblue')
ax2.plot(TIME_BGC,CHL_INT_MOD_AS2_SM, label = 'Community 2 Smooth', color = 'b')
ax2.set_ylabel('Chl-a (mg m$^{-2}$)', fontsize= 7, color='k')
ax2.set_title('(b)', fontsize= 12, color ='k')
ax2.yaxis.set_tick_params(labelsize= 6)
ax2.legend(loc="upper right", fontsize=4)
ax3.plot(TIME_BGC,BVF_INDEX_AV, label = 'Brunt-Väisälä Frequency (BVF)', color = 'mediumpurple')
ax3.set_ylabel('Mean BVF', fontsize= 7, color='k')
ax3.legend(loc="upper right", fontsize=7)
ax3.set_title('(c)', fontsize= 12, color='k')
ax3.set_ylim([0, 0.0002])
ax3.set_yticks([0,0.00005,0.00010,0.00015,0.00020])
ax3.yaxis.set_tick_params(labelsize= 6)
ax4.plot(TIME_BGC,MLD_TEMP_ARRAY,label = 'Mixed-layer depth ($Z_m$)', color = 'orange')
ax4.set_ylabel('$Z_m$ (m)', fontsize= 7, color='k')
ax4.set_xlabel('Time', fontsize= 7, color='k')
ax4.legend(loc="upper right", fontsize=7)
ax4.set_title('(d)', fontsize= 12, color='k')
ax4.yaxis.set_tick_params(labelsize= 6)
ax4.xaxis.set_tick_params(labelsize= 6)
plt.show()

### STEP 25: CORRELATIONS OF INTEGRATED CHL-A AND STRATIFICATION ###
# Community 1 Pearson
XX = BVF_INDEX_AV
YY = CHL_INT_MOD_AS1
asq = np.where(XX > -0.0001)
XX = XX[asq]
YY = YY[asq]
asq = np.where(YY > 0.000)
XX = XX[asq]
YY = YY[asq]
STATS_REG_P1  = stats.pearsonr(XX,YY)
R_CORR_ST         = ("{0:.2f}".format(STATS_REG_P1[0]))
P_CORR_ST         = ("{0:.3f}".format(STATS_REG_P1[1]))
print([R_CORR_ST,P_CORR_ST])
# Community 2 Pearson
XX = BVF_INDEX_AV
YY = CHL_INT_MOD_AS2
asq = np.where(XX > -0.0001)
XX = XX[asq]
YY = YY[asq]
asq = np.where(YY > 0.000)
XX = XX[asq]
YY = YY[asq]
STATS_REG_P1  = stats.pearsonr(XX,YY)
R_CORR_ST         = ("{0:.2f}".format(STATS_REG_P1[0]))
P_CORR_ST         = ("{0:.3f}".format(STATS_REG_P1[1]))
print([R_CORR_ST,P_CORR_ST])

sys.exit()
################# FINISH CODE ##########################