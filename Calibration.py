# Welcome!
# If you're trying to read through the code in this repository, it is recommended
# to read in the following order: 
# 1. Instrument_Creator.py 
# 2. Calibration.py
# 3. DataLoader.py
# 4. Plotting.py


# The calibration procedure is KEY for understanding the prismatic concept. Please dig in 
# to this code and reach out to the author Adit Desai if you have questions.


##Here are the necessary import statements for the code

import numpy as np
import os
import scipy.signal
from lmfit.models import *
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler
from datetime import datetime as dt

## The bulk of the calibration procedure requires only the instrument object created using Instrument_Creator.py
## and the name of the folder where the calibration data is located. Note that the path to the folder should have been
## specified in the pathBase paramter when defining the Instrument() object.

## the other parameters are in case you want to plot the fitted Gaussians. set plot=True and then specify
## the x and y axis-bounds and whether you want to save the figure. plotVals should be passed as a list of
## the energies you want plotted. By default, all energies are plotted. 
def calibration(instrument, folder, plot=False, xlim = None, ylim = None, plotVals = "all", saveFig = False):
    print("sup")
    ## rawDataDict will be a layered dictionary, with the format {pixel1:{Ei1:Intensity, Ei2:Intensity,...}, ...}
    ## essentially each pixel will have a map with each energy that represents the raw data measured from the
    ## the calibration experiment.
    rawDataDict = {}
    ## totalNeutronDict has the format {Ef1:neutrons1(float), Ef2:neutrons2, ...} and keeps track of the total
    ## intensity measured. This will eventually be used for scaleDict
    #totalNeutronDict = {}
    ## Scale Dict will have the format {Ef1:float1, Ef2:float2, .... } and tells how much the intensity in 
    ## rawDataDict should be scaled by per energy. This is done to account for discrepancies in total intensity
    ## measured for different Efs. Essentially, in practice, from a white beam each Ef should be measured equally
    ## but instrument parameters will affect the distribution. This corrects for that.
    #scaleDict = {}
    pixelEdges = np.linspace(-0.45, 0.45, 1025)
    pixels = pixelEdges[:-1] + (pixelEdges[1]-pixelEdges[0])/2
    ## there are 1024 pixels used in the 0.9 meter active length ReuterStokes detector (based off CAMEA paper)
    ## thus 1025 pixel edges are defined. The center of the bin is then defined in pixels, based on
    ## the pixelEdges and the halfway point between the edges.

    ## This setting up the first layer of the dictionary as commented above
    for pixel in pixels:
        rawDataDict[pixel] = {}

    if plot == True:
        ## This section just sets up the colors for the plot. Feel free to ignore
        default_cycler = cycler(color=['#13306dff', '#347537', '#a65c85',  '#4565b9','#d65a5f', '#f56464',
                                '#eb8055', '#f9b64aff', '#efe350'])
        plt.rc('axes', prop_cycle=default_cycler)
        ## this next section takes care of some technicalities.
    for ei in instrument.energyList():
        if instrument.type == "full":
            ## the original McStas simulations had slightly different folder naming conventions (bad practice I know)
            ## and this section takes care of that.
            fileFormat = f"Ei-{ei}"
        elif instrument.type == "toy model":
            fileFormat = f"Ei-{ei}TwoTh15psi0"
        try:
            os.chdir(f"{instrument.pathBase}/{folder}/{fileFormat}")
        except:
            print(f"Could not find the file {fileFormat}!")
            continue
        dataPos = []
        dataIntensities = []
        for detRow in range(1,3):
            for det in range(1,8):
                if detRow == 2 and det == 7:
                    continue
                try:
                    ## Below is the naming convention for the data files. The "_1" indicates the channel of detectors
                    ## However all calibration used only one channel; that is all that is needed
                    fileOpener = open(f"ReuterStokes{detRow}_{det}_1.psd")
                except:
                    continue
                fileData = fileOpener.readlines()
                ## This next section reads the data from the ReuterStokes data and appends the location
                ## the neutron landed on the detctor (ypos) and the measured intensity at the point
                for line in fileData[instrument.startpoint():]:
                    splitLine = line.split()
                    ypos = float(splitLine[2])
                    intensity = float(splitLine[0])
                    if intensity > 0:
                        dataPos.append(ypos)
                        dataIntensities.append(intensity)
        ## Below is where I histogram all the data into 1024 pixels. This is only valid for the McStas simulations
        ## where postion is known absolutely.                
        ## unusedbins is equivalent to binEdges.
        hist, unusedbins = np.histogram(dataPos, bins=1024, weights = dataIntensities, range = (-0.45, 0.45))
        ## The next portion will plot the raw, histogrammed signal measured for each Ei 
        if plot == True:
            if plotVals == "all" or plotVals == "All" or ei in plotVals:
                ## The +0.45 term is essentially to make detector position (0,0.9) rather than (-0.45, 0.45)
                ## for plotting purposes
                plt.scatter(pixels+0.45, hist, marker = "x", s= 15)
        ## The line below finds the peaks across the detector for a given Ef
        ## There can be multiple peaks. The distance will keep peaks that are too close together from
         
        i_pk, _ = scipy.signal.find_peaks(hist, prominence = np.max(hist)/10, distance = 50)
        
        ## the next portion makes use of the lmfit module to assist in fitting
        ## lmfit is useful because it can handle overlapping and multiple peak fitting on the same axis
        ## Essentially each Ef will have multiple peaks as it can be scattered by multiple peaks, thus
        ## lmfit will keep track of all the peaks
        gaussModel = GaussianModel()
        ## Setting up initial guess of the Gaussian fit of the data 
        pars = gaussModel.guess(data=hist, x = pixels)
        modelList = []
        ## This makes sure for each peak i, found in scipy.signal.find_peaks, there is an associated fit 
        for i in range(len(i_pk)):
            peak_index = i_pk[i]
            gauss = GaussianModel(prefix=f'g{i+1}_')
            pars.update(gauss.make_params())
            pars[f'g{i+1}_center'].set(pixels[peak_index])
            pars[f'g{i+1}_sigma'].set(0.01)
            pars[f'g{i+1}_amplitude'].set(hist[peak_index])
            modelList.append(gauss)
        
        ## The two lines below just prepare the module by summing all the data
        modelArray = np.array(modelList)
        model = np.sum(modelArray)
        
        ## Now Lmfit will perform the fit to the raw signal.
        out = model.fit(hist, pars, x=pixels)
        ## The below section controls how the plotting, it will essentially plot the Gaussian fit
        if plot == True:
            if plotVals == "all" or plotVals == "All" or ei in plotVals:

                ## The +0.45 term is essentially to make detector position (0,0.9) rather than (-0.45, 0.45)
                ## out.best_fit is the best fit Gaussian for the associated Ei
                plt.plot(pixels+0.45, out.best_fit, label = '{0:.2f} meV'.format(ei))
        ## mySum keeps track of the total number of neutrons for each Ei for totalNeutronDict
        #mySum = 0
        for index in range(0, len(pixels)):
            ## While we used a Gaussian fit for the neutrons, parts of the regions need to be cut off
            ## Due to the baffle regions. BaffleChecker makes sure that the baffleRegions described in
            ## Instrument_Creator are mapped to 0 intensity
            baffleChecker = False
            for baffleRegion in instrument.baffleRegions():
                if pixels[index] > baffleRegion[0] and pixels[index] < baffleRegion[1]:
                    baffleChecker = True
            if baffleChecker == False:
                ## Then rawDataDict will get the Gaussian fit of the data not in the baffle region
                ## Note that only neutrons landing in the non-baffle regions are included in the sum
                rawDataDict[pixels[index]][ei] = out.best_fit[index]
        #        mySum += out.best_fit[index]
            else:
                rawDataDict[pixels[index]][ei] = 0
        ## totalNeutronDict then keeps the sum of the fitted Gaussian's data           
        #totalNeutronDict[ei] = mySum
    ## Now to define the scaling correcting for the previously described inhomegenous measurement of Ef
    ## the maxNperE is defined from the neutron sums of all Efs. Note that this exclusively used for calibration
    ## this will not extrapolate the experimental data in any fashion.
    #maxNperE = max(totalNeutronDict.values())

    ## After the data is collected, a pandas dataframe is created
    calibrationDF = pd.DataFrame.from_dict(rawDataDict)
    ## After collecting the data, we need to scale it to correct for the previously described inhomegenous measurement of Ef
    ## The first step in finding the scaling is find the net sum of all neutrons measured at the detector
    ## Then, the maximally measured energy is used to set the scale factor
    ## All other energies that are intrinsically scaled less are scaled by N(E_max)/N(E_i)
    calibrationSums = np.array(calibrationDF.sum(axis=1))
    calibrationScales = np.max(calibrationSums)/calibrationSums
    ## Note that the scale factor is exclusively used for calibration
    ## this will not extrapolate the experimental data in any fashion.
    ## Now I multiply the raw data by the scale factors
    calibrationDF = calibrationDF.multiply(calibrationScales, axis=0)
    ## Finally, the calibrationdf will make sure that each pixel has a distribution of 
    ## intensities that sum to 1. Essentially each pixel will have a probability
    ## distribution of energies; during an experiment if a neutron lands in a
    ## pixel, it is split into fractional neutrons with the different Efs corresponding
    ## to the probability distribution of the pixel.
    ## To do this, I first collect the sum of all the neutrons measured at a pixel
    ## Then I divide the pixels by the sum (calibrationSums) such that a probability 
    ## is associated with each Ef.   
    calibrationSums = calibrationDF.sum()
    ## The line below is just a way to avoid having a 0/0 scenario in the event
    ## that a Sum term is 0 (for example in a baffle region)
    calibrationSums[calibrationSums == 0] = 1
    ## Here is the division line.
    calibrationDF = calibrationDF/calibrationSums
    
    ## The rest of the code is for plotting the Gaussian fits
    if plot == True:
        ## control the x and y-axes
        if xlim != None:
            plt.xlim(xlim[0], xlim[1])
        if ylim != None:
            plt.ylim(ylim[0], ylim[1])
        plt.xlabel("Detector Position (m)")
        plt.ylabel("Intensity (a.u.)")
        plt.title(f"{instrument.stations} Stations Mosaic {instrument.mosaic} Calibration")
        plt.legend()
        if saveFig == True:
            ## As a default, the saveFig True will save the output as a pdf in the directory where your 
            ## calibration folder and (simulated) experimental data are located
            plt.savefig(f"{instrument.pathBase}/{instrument.stations}_Stations_Mosaic_{instrument.mosaic}_Calibration_{dt.now().strftime('%Y_%m_%d_%H_%M_%S')}.pdf", format = "pdf")
        plt.show()
    ## I round the calibration dataframe to get rid of extremely small terms in fits
    ## It will ease computational intensity later
    calibrationDF = calibrationDF.round(decimals=5)
    print("Calibration Successful!")
    return calibrationDF

## The pixelHistogram will let you look at an individual pixel and study the distribution
## Running it assumes you already have a calibration complete, which you would then pass as an input
## pixelNum is the nth pixel (an int) you want passed. I may one day add an option that lets you pass
## a yposition on the director (float) and then plot the nearest pixel to the yposition. 
def pixelHistogram(instrument, calibration, pixelNum, xlim = None, ylim = None, saveFig = False):
    ##pixelVals accesses the centers of the pixels
    pixelVals = calibration.columns.values

    pixel = pixelVals[pixelNum]
    ## Below I plot the histogram of the energies weighted by the
    ## Quite simple code, feel free to look at the documentation of plt.hist if you don't understand what is going on
    plt.hist(calibration.index, bins = calibration.index, weights = calibration[pixel])
    plt.xlabel("$E_f$")
    plt.ylabel("Intensity Density")
    plt.title(f"{instrument.stations} Stations Mosaic {instrument.mosaic} Calibration Pixel {pixelNum}")
    if xlim !=None:
        plt.xlim(xlim[0], xlim[1])
    if ylim != None:
        plt.ylim(ylim[0], ylim[1])
    if saveFig == True:
        plt.savefig(f"{instrument.pathBase}/{instrument.stations}_Stations_Mosaic_{instrument.mosaic}_Calibration_Pixel_{pixelNum}_{dt.now().strftime('%Y_%m_%d_%H_%M_%S')}.pdf", format = "pdf")
    plt.show()






            

    