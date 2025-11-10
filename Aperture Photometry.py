#!/usr/bin/env python
# coding: utf-8

# In[1]:


#START OF MAIN PROGRAM
#collect all file names in array
from astropy.io import fits
from astropy.wcs import WCS
from datetime import datetime
import matplotlib.dates as date
import glob
import os
os.scandir(r"C:\Users\Joe\Documents\Uni\Year 3\PX3350 Physics Project\Test Data\WASP-52b-lights (Solved)")
Files = glob.glob('*[b]*.fits')

#2-d array containing all reference star and transit star positions
ref = [[348.3346328,8.7883528],
       [348.3716107,8.7527600],
       [348.3909231,8.6942919],
       [348.3832744,8.6821373],
       [348.5499822,8.8488922],
       [348.5303259,8.6581843],
       [348.4948000,8.7609959]]


#define aperture function

def aperture(ra,dec,image,ap_rad,ap_bkg):
    #open image
    file = fits.open(image)
    #pixel values from image
    img_data = file[0].data
    #convert ra and dec to pixel positions for given image
    w = WCS(file[0].header)
    x_pos,y_pos = w.world_to_pixel_values(ra,dec)
    #to nearest pixel(whole number)
    x=int(x_pos)
    y=int(y_pos)
    
    #define sum and count for aperture and background
    ap_sum = 0
    N_ap = 0
    bkg_sum = 0
    N_bkg = 0  
    
    #loop over background radius (includes aperture and background)
    for i in range(-ap_bkg,ap_bkg+1):
        for j in range(-ap_bkg,ap_bkg+1):
            #if pixel position is inside aperture but not inside background, sum the pixel values and pixel number
            if (i**2+j**2<=ap_rad**2):
                ap_sum+=img_data[y+i,x+j]
                N_ap+=1
            #if pixel position is inside backgoudn but not inside aperture, sum and count pixels
            if ((i**2+j**2<=ap_bkg**2)&(i**2+j**2>ap_rad**2)):
                bkg_sum+=img_data[y+i,x+j]
                N_bkg+=1
             
    #average value for background pixel
    avg_bkg = bkg_sum/N_bkg
    #flux = aperture sum minus background for each pixel in aperture
    flux = ap_sum-(avg_bkg*N_ap)
    file.close()
    return flux
        
#ap_sum: sum of pixel magnitudes in central aperture 
#N_ap: number of pixels in central aperture
#N_bkg: number of pixels in background 
#bkg_sum: sum of pixels magnitudes in background 


# In[2]:


import matplotlib.pyplot as plt
import numpy as np

#time array point for each image
time = np.arange(0,len(Files),1)

#create array(7,63) of zeros to fill
#7 columns for each star
#63 rows for each image
flux = np.zeros((len(ref),len(Files)))

#for i over each star
#for j over all image files
for i in range(0,len(ref)):
    for j in range(0,len(Files)):
        #fill flux array using aperture function
        flux[i][j] = aperture(ref[i][0],ref[i][1],Files[j],20,40) #defined size of aperture and annulus   


# In[3]:


from scipy.optimize import curve_fit
from datetime import datetime
import matplotlib.dates as date
from astropy.time import Time


ref_frame = 16 #select the reference frame

multi_factor = np.zeros((len(ref)-1,len(Files)))

#loop over all reference stars 
for i in range(0,len(ref)-1):
    multi_factor[i] = (flux[i]/flux[i,ref_frame])#fill multi factor array

#zeros array to fill (length 63)
multi_factor_avg = np.zeros(len(Files))


for j in range(0,len(Files)):
    for i in range(0,len(ref)-1):
        multi_factor_avg[j] += (multi_factor[i,j]/(len(ref)-1))

print(multi_factor[0,8])

time = np.zeros(len(Files)).astype(datetime)
for j in range(0,len(Files)):
    file = fits.open(Files[j])
    hdr = file[0].header
    time[j] = datetime.fromisoformat(hdr['DATE-AVG']) 
    
date = date.date2num(time)
#print(time)
#t = Time(time, format='isot', scale='utc')
#print(t.jd)

#zeros array to fill with normalised flux (same size as flux array)       
normalised = np.zeros((len(ref),len(Files)))

#apply this average to the target star
normalised[len(ref)-1] = (flux[len(ref)-1]/flux[len(ref)-1,ref_frame])/multi_factor_avg
normalised_flux = normalised[len(ref)-1]

#plot target star flux 
plt.figure(figsize=(9, 5))
plt.title("Transit Light Curve")
plt.xlabel("Time [dd hh:mm]")
plt.ylabel("Relative Flux")
plt.ylim(0.96,1.02)
plt.plot(date,normalised_flux,'.', markersize=10)
plt.grid()


# In[9]:


############################# ERROR BARS ########################################
def std(x,mu,N):
    return np.sqrt(((np.sum(x-mu))**2/N))


#loop over all reference stars 
for i in range(0,len(ref)-1):
    multi_factor[i] = (flux[i]/flux[i,ref_frame])#fill multi factor array

#zeros array to fill (length 63)
standard_dev = np.zeros(len(Files))

#i = current image
#j = current star
for i in range(len(Files)):
    for j in range(len(ref)-1):
         standard_dev[i] = std(multi_factor[j,i],multi_factor_avg[i],len(ref)-1)

standard_err = standard_dev/np.sqrt(6)
std = np.std(multi_factor,axis=0)
std_err = std/np.sqrt(6)

#plot target star flux 
plt.figure(figsize=(9, 5))
plt.title("Transit Light Curve")
plt.xlabel("Time [dd hh:mm]")
plt.ylabel("Relative Flux")
plt.ylim(0.96,1.02)

plt.errorbar(time,normalised_flux,yerr=std_err,fmt ='.',markersize=10)
#plt.plot(date,normalised_flux,'.', markersize=10)
plt.grid()
print(std_err)
#standard_error = standard_dev/np.sqrt()


# In[5]:


############################## DETRENDING #########################################
def func(x,m,c):
    return m*x+c
    
#plot normalised flux 
t = np.arange(63)
plt.plot(t,normalised_flux,'.')

#plot only the edges of the transit 
time21 = np.arange(21)
time13 =np.arange(50,63)
time = np.concatenate((time21,time13))
flux_to_fit = np.concatenate((normalised_flux[t<21],normalised_flux[t>49]))
plt.plot(time,flux_to_fit,'.')


#run curve fit on the edges of the transit
m_guess = 0.012/34
c_guess = 0.996
constants = curve_fit(func, time, flux_to_fit)
m_fit = constants[0][0]
c_fit = constants[0][1]


print(m_fit,c_fit)


# In[6]:


#plot the curve fit of the edges of the transit
fit_flux = func(t,m_fit,c_fit)

plt.figure(figsize=(9, 5))
plt.plot(t,fit_flux,label='detrend fit')
plt.plot(t,normalised_flux,'.',label='data')
#plot target star flux 

plt.title("Transit Light Curve")
plt.xlabel("Image")
plt.ylabel("Relative Flux")
plt.ylim(0.96,1.02)
plt.grid()
plt.legend(loc='lower left')


# In[11]:


#apply the trend to the whole transit a
detrend_normalised_flux = (normalised_flux - fit_flux) + 1

time = np.zeros(len(Files)).astype(datetime)
for j in range(0,len(Files)):
    file = fits.open(Files[j])
    hdr = file[0].header
    time[j] = datetime.fromisoformat(hdr['DATE-AVG'])  

#plot target star flux 
plt.figure(figsize=(9, 5))
plt.title("Transit Light Curve")
plt.xlabel("Time [dd hh:mm]")
plt.ylabel("Relative Flux")
plt.ylim(0.96,1.02)
plt.errorbar(time,detrend_normalised_flux,yerr=std_err,fmt ='.',markersize=10)
plt.grid()


# In[141]:


################################### FITTTING ######################################
import batman
import numpy as np
import matplotlib.pyplot as plt

params = batman.TransitParams()       #object to store transit parameters
params.t0 = 0.004                        #time of inferior conjunction
params.per = 1.75                       #orbital period
params.rp = 0.16                       #planet radius (in units of stellar radii)
params.a =  15                       #semi-major axis (in units of stellar radii)
params.inc = 88                     #orbital inclination (in degrees)
params.ecc = 0.                       #eccentricity
params.w = 90.                        #longitude of periastron (in degrees)
params.limb_dark = "nonlinear"        #limb darkening model
params.u = [0.5, 0.1, 0.1, -0.1]      #limb darkening coefficients [u1, u2, u3, u4]

t = np.linspace(-0.03, 0.03, 63)  #times at which to calculate light curve
m = batman.TransitModel(params, t)    #initializes model


# In[142]:


new_flux = m.light_curve(params)    #recalculates light curve
plt.figure(figsize=(12, 5))     
plt.plot(t,new_flux)
plt.errorbar(t,detrend_normalised_flux,yerr=std_err,fmt ='.',markersize=10)
#plt.plot(t,detrend_normalised_flux,'.', markersize=10)
plt.title("Transit Light Curve")
plt.xlabel("Time [dd hh:mm]")
plt.ylabel("Relative Flux")
plt.ylim(0.96,1.02)
plt.grid()


# In[8]:


from astropy.timeseries import BoxLeastSquares
import batman
from scipy.optimize import minimize
#import corner
import numpy as np
import matplotlib.pyplot as plt
#import emcee
import astropy.units as u


# In[106]:


#max_power = np.argmax(periodogram.power)
period_bls = 1.75
epoch_bls = 0
depth_bls = 0.0269
#duration_bls = periodogram.duration[max_power]
t = np.linspace(-0.03, 0.03, 63)  #times at which to calculate light curve
depth_value = np.zeros(511)
scaled_seperation_value = np.zeros(511)
inclination_value = np.zeros(511)


# In[109]:


def log_likelihood(parameters, times, observed_flux):
    global count

    depth_value[count] = parameters[2]
    scaled_seperation_value[count] = parameters[3]
    inclination_value[count] = parameters[4]
    # Generate the model flux using the batman package and the provided parameters.
    model_flux = f_batman(times, *parameters)
    
    # Calculate the variance of the observed flux (squared error).
   
    # Compute the log-likelihood using the Gaussian log-likelihood formula.
    log_likelihood_value = -0.5 * np.sum(
        (observed_flux - model_flux) ** 2 
    )
    
    count +=1
    return log_likelihood_value


def f_batman(
    times,
    time_of_conjunction,
    orbital_period,
    planet_radius,
    semi_major_axis,
    orbital_inclination,
    baseline_flux=0.0,
    eccentricity=0,
    longitude_of_periastron=90,
    limb_darkening_coefficients=[0.3, 0.28],
    limb_darkening_model="quadratic",
):
   
    # Initialize the parameters for the batman transit model
    params = batman.TransitParams()
    params.t0 = time_of_conjunction
    params.per = orbital_period
    params.rp = planet_radius
    params.a = semi_major_axis
    params.inc = orbital_inclination
    params.ecc = eccentricity
    params.w = longitude_of_periastron
    params.u = limb_darkening_coefficients
    params.limb_dark = limb_darkening_model

    # Initialize the batman transit model and calculate the light curve
    transit_model = batman.TransitModel(params, times)
    flux_model = transit_model.light_curve(params)

    # Add the baseline flux to the model output and return
    return np.array(flux_model) + baseline_flux


# In[110]:


# Constants for conversion
SOLAR_RADIUS_IN_AU = 215.0  # Solar radius in astronomical units

# Calculation of the semi-major axis in solar radii assuming a 1 Msol, 1 Rsol host
semi_major_axis = ((period_bls / 365.0) ** (2.0 / 3.0)) * SOLAR_RADIUS_IN_AU
semi_major_axis =7.5
# Initial guess for the transit model parameters
initial_guess = np.array(
    [
        epoch_bls,  # Time of central transit
        period_bls,  # Orbital period
        depth_bls**0.5,  # Approximation of planet radius
        semi_major_axis,  # Semi-major axis
        85,  # Orbital inclination (degrees)
        0.0,  # Baseline flux offset
    ]
)

count = 0
# Define the negative log likelihood function for optimization
negative_log_likelihood = lambda *args: -log_likelihood(*args)

# Perform optimization to find the best-fit parameters
solution = minimize(
    negative_log_likelihood,
    initial_guess,
    args=(t, detrend_normalised_flux), 
)

# Parameter labels for output
parameter_labels = [
    "t0 (Epoch)",
    "per (Period)",
    "rp (Planet Radius)",
    "a (Semi-major Axis)",
    "inc (Inclination)",
    "baseline (Flux Offset)",
]


# Print the header of the table
header = f"{'Parameter':<25} {'Initial Value':<15} {'Fitted Value':<15}"
print(header)
print("-" * len(header))
#solution.x[1]=1.75
# Iterate over each parameter and print the values in a table-like format
for idx in range(len(solution.x)):
    parameter_label = parameter_labels[idx]
    initial_value = initial_guess[idx]
    fitted_value = solution.x[idx]

    row = f"{parameter_label:<25} {initial_value:<15.4f} {fitted_value:<15.4f}"
    print(row)


# In[128]:


std_dev_depth = np.std(depth_value)
std_dev_scaled_seperation = np.std(scaled_seperation_value)
std_dev_inclination = np.std(inclination_value)


print(std_dev_depth)
print(std_dev_scaled_seperation)
print(std_dev_inclination)


# In[143]:


# Constants
NUMBER_OF_POINTS = 10000
PHASE_FOLD_THRESHOLD = 0.5

# Create a time array for the transit model
model_times = np.linspace(np.min(t), np.max(t), NUMBER_OF_POINTS)

# Compute the initial model flux and the fitted model flux
initial_model_flux = f_batman(model_times, *initial_guess)
fitted_model_flux = f_batman(model_times, *solution.x)

# Phase-fold the observed data time-array
phase_observed = (t - solution.x[0]) % solution.x[1] / solution.x[1]
phase_observed[phase_observed > PHASE_FOLD_THRESHOLD] -= 1

# Phase-fold the BLS model time-array and sort
phase_bls_model = (model_times - epoch_bls) % period_bls / period_bls
phase_bls_model[phase_bls_model > PHASE_FOLD_THRESHOLD] -= 1
sorted_bls_model_idx = np.argsort(phase_bls_model)

# Phase-fold the fitted model time-array and sort
phase_fitted_model = (model_times - solution.x[0]) % solution.x[1] / solution.x[1]
phase_fitted_model[phase_fitted_model > PHASE_FOLD_THRESHOLD] -= 1
sorted_fitted_model_idx = np.argsort(phase_fitted_model)

# Plotting
# Time series plot of data and models
plt.figure(figsize=(10, 5))
plt.plot(t, detrend_normalised_flux, "o", color="xkcd:dusty red", label="Data (Time Series)", alpha=0.4)
#plt.plot(model_times, initial_model_flux, "k", label="Initial BLS Model", alpha=0.6)
plt.plot(model_times, fitted_model_flux, "--k", label="Fitted Model")
plt.xlabel("Time [days]")
plt.ylabel("Normalized Flux")
plt.legend()
plt.grid

# Phase-folded plot of data and fitted model
plt.figure(figsize=(12, 5))
#plt.errorbar(date,normalised_flux,yerr=std_err,fmt ='.',markersize=10)
plt.errorbar(
    phase_observed * solution.x[1],
    detrend_normalised_flux,
    yerr=std_err,
    fmt = "o",
    color="xkcd:dusty red",
    label="Data",
    markersize=5,
    alpha=0.9,
)
plt.plot(
    phase_fitted_model[sorted_fitted_model_idx] * solution.x[1],
    fitted_model_flux[sorted_fitted_model_idx],
    "--k",
    label="Fitted Model",
)
plt.xlabel("Time from Central Transit [days]")
plt.ylabel("Normalized Flux")
plt.xlim(-0.04, 0.03)
plt.legend()
plt.grid()


# In[72]:


##########################################################################################################################

#########################################################################################################################


# In[10]:


from astropy.timeseries import BoxLeastSquares
import batman
from scipy.optimize import minimize
from scipy.stats import norm
#from IPython.display import display, Math
#import corner
import numpy as np
import matplotlib.pyplot as plt
#import emcee
import astropy.units as u
#from multiprocessing import Pool
t = np.linspace(-0.03, 0.03, 63)  #times at which to calculate light curve


# In[13]:


from scipy import stats
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt


n = 20


print(detrend_normalised_flux)

def gaussian(params):
    mean = params[0]   
    sd = params[1]

    # Calculate negative log likelihood
    nll = -np.sum(stats.norm.logpdf(detrend_normalised_flux, loc=mean, scale=sd))

    return nll


initParams = [1, 1]

results = minimize(gaussian, initParams, method='Nelder-Mead')
print(results.x)


# In[14]:


from scipy import stats
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

known_log_pdfs = stats.norm.logpdf(detrend_normalised_flux)
print(known_log_pdfs)
def gaussian(params):
    mean = params[0]   
    sd = params[1]

    estimated_log_pdfs = stats.norm.logpdf(detrend_normalised_flux, loc=mean, scale=sd)
    
    mse = ((estimated_log_pdfs - known_log_pdfs) ** 2).mean()

    return mse

initParams = [1, 1]

results = minimize(gaussian, initParams, method='Nelder-Mead')
print(results.x)


# In[123]:


#max_power = np.argmax(periodogram.power)
period_bls = 1.75
epoch_bls = 0
depth_bls = 0.0269
duration_bls = 1.8

def log_likelihood(parameters, times, observed_flux, observed_flux_error):
    """
    Calculates the logarithm of the likelihood for a given model and observed data.

    This function computes the log-likelihood for a planetary transit model, assuming
    Gaussian noise in the observed data. It uses the 'batman' package to generate
    the transit model.

    Parameters:
    parameters : array_like
        An array of parameters for the transit model. These typically include the
        time of central transit, orbital period, planet radius, semi-major axis,
        orbital inclination, eccentricity, argument of periastron, and baseline flux.
    times : array_like
        The time points at which observations were made.
    observed_flux : array_like
        The observed flux values corresponding to the given times.
    observed_flux_error : array_like
        The standard error of each observed flux measurement.

    Returns:
    float
        The log-likelihood of the observed data given the model parameters.
    """
    # Generate the model flux using the batman package and the provided parameters.
    model_flux = f_batman(times, *parameters)

    # Calculate the variance of the observed flux (squared error).
    
    observed_flux_error[16]=0.0015
    
    variance = observed_flux_error**2
    
    # Compute the log-likelihood using the Gaussian log-likelihood formula.
   
    log_likelihood_value = (-0.5 * np.sum((observed_flux - model_flux) ** 2 )/ variance + np.log(2.0 * np.pi * variance))

    return log_likelihood_value


def f_batman(
    times,
    time_of_conjunction,
    orbital_period,
    planet_radius,
    semi_major_axis,
    orbital_inclination,
    baseline_flux=0.0,
    eccentricity=0,
    longitude_of_periastron=90,
    limb_darkening_coefficients=[0.34, 0.28],
    limb_darkening_model="quadratic",
):
    """
    Computes the transit light curve model using the batman package.

    This function generates a transit light curve for an exoplanet based on the specified
    transit parameters and times. It uses the batman package to create a model of how
    the flux from a star decreases as a planet transits in front of it.

    Parameters:
    times : array_like
        The time points at which to calculate the light curve.
    time_of_conjunction : float
        The time of central transit or inferior conjunction.
    orbital_period : float
        The orbital period of the exoplanet.
    planet_radius : float
        The radius of the planet in units of the host star's radius.
    semi_major_axis : float
        The semi-major axis of the planet's orbit, in stellar radii.
    orbital_inclination : float
        The orbital inclination in degrees.
    baseline_flux : float, optional
        The baseline flux level to add to the model output. Default is 0.0.
    eccentricity : float, optional
        The eccentricity of the planet's orbit. Default is 0.
    longitude_of_periastron : float, optional
        The longitude of periastron in degrees. Default is 90.
    limb_darkening_coefficients : list, optional
        The limb darkening coefficients [u1, u2]. Default is [0.34, 0.28].
    limb_darkening_model : str, optional
        The model of limb darkening to use. Default is "quadratic".

    Returns:
    numpy.ndarray
        An array of flux values corresponding to each time point.
    """
    # Initialize the parameters for the batman transit model
    params = batman.TransitParams()
    params.t0 = time_of_conjunction
    params.per = orbital_period
    params.rp = planet_radius
    params.a = semi_major_axis
    params.inc = orbital_inclination
    params.ecc = eccentricity
    params.w = longitude_of_periastron
    params.u = limb_darkening_coefficients
    params.limb_dark = limb_darkening_model

    # Initialize the batman transit model and calculate the light curve
    transit_model = batman.TransitModel(params, times)
    flux_model = transit_model.light_curve(params)

    # Add the baseline flux to the model output and return
    return np.array(flux_model) + baseline_flux


# In[124]:


# Constants for conversion
SOLAR_RADIUS_IN_AU = 215.0  # Solar radius in astronomical units

# Calculation of the semi-major axis in solar radii assuming a 1 Msol, 1 Rsol host
semi_major_axis = ((period_bls / 365.0) ** (2.0 / 3.0)) * SOLAR_RADIUS_IN_AU

# Initial guess for the transit model parameters
initial_guess = np.array(
    [
        epoch_bls,  # Time of central transit
        period_bls,  # Orbital period
        depth_bls**0.5,  # Approximation of planet radius
        semi_major_axis,  # Semi-major axis
        88.5,  # Orbital inclination (degrees)
        0.0,  # Baseline flux offset
    ]
)



# Define the negative log likelihood function for optimization
negative_log_likelihood = lambda *args: -log_likelihood(*args)


# Perform optimization to find the best-fit parameters
solution = minimize(
    negative_log_likelihood,
    initial_guess,
    args=(t, detrend_normalised_flux, std_err),
)

# Parameter labels for output
parameter_labels = [
    "t0 (Epoch)",
    "per (Period)",
    "rp (Planet Radius)",
    "a (Semi-major Axis)",
    "inc (Inclination)",
    "baseline (Flux Offset)",
]


# Print the header of the table
header = f"{'Parameter':<25} {'Initial Value':<15} {'Fitted Value':<15}"
print(header)
print("-" * len(header))

# Iterate over each parameter and print the values in a table-like format
for idx in range(len(solution.x)):
    parameter_label = parameter_labels[idx]
    initial_value = initial_guess[idx]
    fitted_value = solution.x[idx]

    row = f"{parameter_label:<25} {initial_value:<15.4f} {fitted_value:<15.4f}"
    print(row)


# In[125]:


# Constants
NUMBER_OF_POINTS = 10000
PHASE_FOLD_THRESHOLD = 0.5

# Create a time array for the transit model
model_times = np.linspace(np.min(t), np.max(t), NUMBER_OF_POINTS)

# Compute the initial model flux and the fitted model flux
initial_model_flux = f_batman(model_times, *initial_guess)
fitted_model_flux = f_batman(model_times, *solution.x)

# Phase-fold the observed data time-array
phase_observed = (t - solution.x[0]) % solution.x[1] / solution.x[1]
phase_observed[phase_observed > PHASE_FOLD_THRESHOLD] -= 1

# Phase-fold the BLS model time-array and sort
phase_bls_model = (model_times - epoch_bls) % period_bls / period_bls
phase_bls_model[phase_bls_model > PHASE_FOLD_THRESHOLD] -= 1
sorted_bls_model_idx = np.argsort(phase_bls_model)

# Phase-fold the fitted model time-array and sort
phase_fitted_model = (model_times - solution.x[0]) % solution.x[1] / solution.x[1]
phase_fitted_model[phase_fitted_model > PHASE_FOLD_THRESHOLD] -= 1
sorted_fitted_model_idx = np.argsort(phase_fitted_model)

# Plotting
# Time series plot of data and models
plt.figure(figsize=(10, 5))
plt.plot(t, detrend_normalised_flux, "o", color="xkcd:dusty red", label="Data (Time Series)", alpha=0.4)
plt.plot(model_times, initial_model_flux, "k", label="Initial BLS Model", alpha=0.6)
plt.plot(model_times, fitted_model_flux, "--k", label="Fitted Model")
plt.xlabel("Time [days]")
plt.ylabel("Normalized Flux")
plt.legend()

# Phase-folded plot of data and fitted model
plt.figure(figsize=(10, 5))
plt.plot(
    phase_observed * solution.x[1],
    detrend_normalised_flux,
    ".",
    color="xkcd:dusty red",
    label="Phase Folded Data",
    alpha=0.1,
)
plt.plot(
    phase_fitted_model[sorted_fitted_model_idx] * solution.x[1],
    fitted_model_flux[sorted_fitted_model_idx],
    "--k",
    label="Fitted Model",
)
plt.xlabel("Time from Central Transit [days]")
plt.ylabel("Normalized Flux")
plt.xlim(-1, 1)
plt.legend()


# In[ ]:





# In[ ]:




