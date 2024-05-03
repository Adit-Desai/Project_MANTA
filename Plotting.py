# Welcome!
# If you're trying to read through the code in this repository, it is recommended
# to read in the following order: 
# 1. Instrument_Creator.py 
# 2. Calibration.py
# 3. DataLoader.py
# 4. Plotting.py

## Necessary import statements
import matplotlib.pyplot as plt
import numpy as np
import scipy
from lmfit.models import * 
from datetime import datetime as dt

## Now you can take a look at the various plotting features included in this library
## It's currently designed for simulations and for the simple cubic sample used, and such
## only plots Q and E, but does this defined off the neutron's Q, rather than the sample

## When calling plotting functions, it assumes the variable names you call are in the exact form as listed
## in the pandas dataframe. Please print out a portion of the DataFrame if you are uncertain
## about how the variables are named

## cut2D requires the instrument and the dataframe prepared in the DataLoader.py
## It also will require the x-axis and y-axis variable, xVar and yVar you wish to plot (likely Qx, Qy, E).
## The color of the plot will always be the intensity of the signal.
## This will be plotted such that the third variable you are not putting on your axis, integrationVar,
## is held at a constant value intergrationVal with an acceptance of integrationWidth.
## the optional variables xlim, ylim, and colorBarLim, control the limits of the axes
## you can also set the saveFile = True to save the file as a pdf in the directory
## where the data is located.
def cut2D(instrument, dataframe, xVar, yVar, integrationVar, integrationVal, integrationWidth, 
          xlim = None, ylim = None, colorBarLim = None, saveFile= False):
    ## First we access the relevant data within the integration Volume

    ## Include the try except clause in case there was a mistake in 
    ## xVar or yVar
    try:
        data = dataframe[dataframe[integrationVar] > integrationVal-integrationWidth]
    except:
        print(f"{integrationVar} was not recognized as a variable within the dataframe!")
        return None
    data = data[data[integrationVar] < integrationVal + integrationWidth]
    ## Next we sort the values such that most intense is plotted last
    ## this is essential for the scatterplot method used
    data = data.sort_values(by="Intensity")
    
    ##Now we plot using plt.scatter, with vmin and vmax controlling the colorbar intensityh
    try:
        if colorBarLim != None:
            ##c  controls the color such that it corresponds to the intensity of the event.
            plt.scatter(data[xVar], data[yVar], c=data["Intensity"], s=8, vmin=colorBarLim[0], 
                        vmax = colorBarLim[1])
        else:
            plt.scatter(data[xVar], data[yVar], c=data["Intensity"], s=8)
    except:
        print(f"{xVar} and or {yVar} were not recognized as variables within the dataframe!")
    plt.colorbar(label = "Intensity (a.u.)")
    ## The rest just controls the axes labels and makes it so they use LaTeX font
    ## if applicable
    if xlim !=None:
        plt.xlim(xlim[0], xlim[1])
    if ylim != None:
        plt.ylim(ylim[0], ylim[1])
    if xVar == "Qy":
        plt.xlabel("$Q_y~(\AA^{-1})$")
    elif xVar == "Qx":
        plt.xlabel("$Q_x~(\AA^{-1})$")
    elif xVar == "E":
        plt.xlabel("$E$ (meV)")
    else:
        plt.xlabel(f"{xVar}")
    ##
    if yVar == "Qy":
        plt.ylabel("$Q_y~(\AA^{-1})$")
    elif yVar == "Qx":
        plt.ylabel("$Q_x~(\AA^{-1})$")
    elif yVar == "E":
        plt.ylabel("$E$ (meV)")
    else:
        plt.ylabel(f"{yVar}")
    ## This just sets the title to the following default
    ## Can be customized according to a users discretion
    plt.title(f"{instrument.stations} Stations Mosaic {instrument.mosaic} {integrationVar} = {integrationVal} $\pm$ {integrationWidth}")
    ## Saves the file if desired by the user
    if saveFile == True:
        plt.savefig(f"{instrument.pathBase}/{instrument.stations}_Stations_Mosaic_{instrument.mosaic}_{integrationVar}_{integrationVal}_{Intensity}_{dt.now().strftime('%Y_%m_%d_%H_%M_%S')}.pdf", format = pdf)
    plt.show()


## The next function creates a 1D plot with one variable where the y-axis is the intensity
## The data is automatically fit to a Gaussian, which lends itself to resolution calculations
## Because this is a scatter plot that is fit to a Gaussian, with a very large number of marginally different
## Q points, the x-axis data is always histogrammed.
## Thus the binSize needs to be specified, which controls the size of histogram binning
## Then you always have to specify the two variables you are integrating over, their
## fixed value, and their integration volume
## The threshold value is an option to help prevent overfitting of the Gaussian and prevents
## small, noisy peaks from being fitted
## the binRange variable not only controls the x-axis plotting region, 
## but also controls the region for which Gaussians are fit to. Note the region shown in the plot
## of binRange plotted is slightly larger than the bounds of binRange just for visual purposes.
## This is useful in particular if you only want 
## A singular Gaussian fit in a specific range to find the fit of a specific peak
## It is particularly useful for the resolution function
## binRange was chosen as the name rather than xlim to help differentiate when using the resolution()
## function.
## ylim is not as special and exclusively controls the y-axis range
## showPlot lets you turn off the plotting function, particularly useful for the resolution()
## function below.
def cut1D(instrument, dataframe, xVar, binSize, integrationVar1, integrationVal1, integrationWidth1, integrationVar2, 
          integrationVal2, integrationWidth2, threshold = None, binRange = None, ylim = None,
          showPlot = True, saveFile=False):
    ## Here we extract the integration region that's valid.
    ## Try except blocks are included in case integrationVar1 or integrationVar2 are incorrectly named.
    try:
        data = dataframe[dataframe[integrationVar1] > integrationVal1-integrationWidth1]
    except:
        print(f"{integrationVar1} was not recognized as a variable within the dataframe!")
        return None

    data = data[data[integrationVar1] < integrationVal1 + integrationWidth1]
    try:
        data = data = data[data[integrationVar2] > integrationVal2-integrationWidth2]
    except:
        print(f"{integrationVar2} was not recognized as a variable within the dataframe!")
        return None
    
    data = data[data[integrationVar2] < integrationVal2 + integrationWidth2]
    ## Here we control the minimum regions to be plotted based off the binRange
    ## This is needed for the gaussian fit
    if binRange == None:
        minVal, maxVal = data[xVar].min(), data[xVar].max()

    else:
        minVal, maxVal = binRange[0], binRange[1]
    ## Now the data is histogrammed based off the intensity. This is to help the identification of a clear Gaussian 
    ## peak. the binsizes are controlled by binSize, and the binRange specified.
    histData, binEdges = np.histogram(data[xVar], weights=data["Intensity"], bins = np.arange(minVal, maxVal, binSize))
    binCenters = binEdges[:-1]  + (binEdges[1] - binEdges[0])
    ## Now the index of the peak is found, which is essential for the Gaussian fitting
    ## The prominence term controls the minimum height it will look for for fitting
    ## The distance variable sets the minimum distance between peaks, and is there to help prevent
    ## overfitting. This is the the most nitpicky of the variables, and is highly dependent on 
    ## where you are fitting. Currently I have it set such that the distance scales automatically
    ## with the number of bins, but the exact scale factor is tricky. Optimization may be needed
    ## here.
    i_pk, _ = scipy.signal.find_peaks(histData, distance = len(binCenters)//6, prominence = threshold)
    
    ## Now we prepare the gaussian fitting package using LmFit.
    gaussModel = GaussianModel()
    ## this is the natural sequence for looking at multiple Gaussians.
    pars = gaussModel.guess(data=histData, x = binCenters)
    modelList = []
    ## setting up the initial guesses which it will refine from,.
    ## the procedure is identical to that described in Calibration.py
    ## For more details refer to there.
    for i in range(len(i_pk)):
        peak_index = i_pk[i]
        gauss = GaussianModel(prefix=f'g{i+1}_')
        pars.update(gauss.make_params())
        pars[f'g{i+1}_center'].set(binCenters[peak_index])
        pars[f'g{i+1}_sigma'].set(0.1)
        pars[f'g{i+1}_amplitude'].set(histData[peak_index])
        modelList.append(gauss)
    modelArray = np.array(modelList)
    model = np.sum(modelArray)
    ## this is the actual fitting procedure
    try:
        out = model.fit(histData, pars, x=binCenters)
    except:
        ## If the fit fails, which can happen particularly if your binRange doesn't capture any peaks 
        ## (slope=0) and so itll ask to check on that.
        print(f"Fit Failed for {integrationVar1}={integrationVal1} {integrationVar2}={integrationVal2}")
        if binRange != None:
            print("Please check if your bin range is large enough!")
            return None
        ## now this controls the plotting
    if showPlot == True:
        ## First the histogrammed raw data is plotted, and then the gaussian fit is overplotted
        plt.scatter(binCenters, histData, marker = "x")
        plt.plot(binCenters, out.best_fit)

        plt.xlabel(f"{xVar}")
        plt.ylabel("Intensity (a.u.)")
        plt.title(f"{instrument.stations} Stations Mosaic {instrument.mosaic}  {xVar} vs. Intensity {integrationVar1} = {integrationVal1} $\pm$ {integrationWidth1} {integrationVar2} = {integrationVal2} $\pm$ {integrationWidth2}")
        ## This just gives the x-axis plotting size a little bit of extra space so the edges aren't defined by
        ## binRange.
        spacing = (maxVal-minVal)*0.05
        plt.xlim(minVal-spacing, maxVal+spacing)
        if ylim != None:
            plt.ylim(ylim[0], ylim[1])
        ## Figures saved as pdfs by default if saveFig set to true.
        if saveFile == True:
            plt.savefig(f"{instrument.pathBase}/{instrument.stations}_Stations_Mosaic_{instrument.mosaic}_{xVar}_v_{Intensity}_{dt.now().strftime('%Y_%m_%d_%H_%M_%S')}.pdf", format = "pdf")
        plt.show()
    return out.best_values

## The resolution function depends heavily on the cut1D function
## Essentially, it sweeps over xVar and performs a resolution calculation for
## resVar at each point in xVar. xVar is integrated over a specific range specified
## by xStepSize at each point in the sweep. resVar will then be the x-axis of the 
## cut1D plot, and the outputted Gaussian fit's Full-width-half-maximum (FWHM)
## will be the y-axis of the output of resolution.

## Again the resolution function requires the instrument, the dataframe
## xVar, xStepSize, resVar, and binSize were all described above.
## Whichever variable is not being swept over is passed as integrationVar,
## and is held constant at integrationVal with integration volume
## integrationWidth. The threshold and binRange are optional and passed directly to the cut1D()
## function. If you would like to see each individual cut produced by cut1D(),
## then set showCuts =True. If you'd like to save all the cut1d() plots, then set saveCuts=True,
## which is passed directly to cut1D(). saveFile is the parameter
## that controls whether the resolution plot (xVar vs resVar) itself is saved.
## ylim controls the y-axis scale, but xlim will also control the number of points
## cut1D is calculated at. Essentially the points sweeped are in range(xlim[0], xlim[1], xStepSize)
## The actual plotted x-axis range is slightly larger than the specified range.
def resolution(instrument, dataframe, xVar, xStepSize, resVar,
                binSize, integrationVar, integrationVal, integrationWidth,   
                threshold = None, binRange = None, xlim = None, ylim = None,
                showCuts = False, saveCuts=False, saveFile=False):
    
    ## The below lists will be appended to and plotted
    xVarList = []
    resList = []

    ## This controls the range which the resolution is calculated
    ## If not specified, it'll just go to the min and max values within the dataframe
    if xlim == None:
        xMin, xMax = dataframe[xVar].min(), dataframe[xVar].max()
    else:
        xMin, xMax = xlim[0], xlim[1] 
    for num in np.arange(xMin, xMax, xStepSize):
        try:
            ## Now the outputted fit from cut1D, which will output a 
            ## dictionary with the best fit parameters
            ## is below
            bestFit = cut1D(instrument=instrument, dataframe=dataframe, 
                            xVar = resVar, binSize = binSize,
                            integrationVar1=integrationVar, integrationVal1 = integrationVal, 
                            integrationWidth1 = integrationWidth, integrationVar2 = xVar,
                            integrationVal2= num, integrationWidth2 = xStepSize/2, 
                            threshold=threshold, binRange = binRange,
                            showPlot=showCuts, saveFile=saveCuts)
        except:
            ## Some values may not work, so the points it fails at are printed. However, in some cases
            ## this is quite normal so the loop will continue instead of breaking.
            print(f"Resolution calculation of {resVar} failed at  {xVar} = {num}, {integrationVar} = {integrationVal}, ")
            continue
        if bestFit== None:
            continue    
        try:
            ## As the fit will automatically try fitting multiple Gaussians
            ## the secondPeak term will warn you of this. It will prevent
            ## any point in which multiple peaks are found from being plotted.
            ## this is when the binRange and threshold parameters are extremely useful
            ## I recommend plotting to see what causes this.
            secondPeak = bestFit["g2_sigma"]
            print("Warning! Multiple Gaussian peaks were found at the same value of "\
                f"{xVar} = {num}. Please change your xRange or the threshold"\
                    " variable to make sure there is only one peak used for "\
                        "calculating the resolution!")
            continue
        except:
            pass
        ## If the fit is successful, the x point and the y point are appended to the list

        xVarList.append(num)
        resList.append(bestFit['g1_sigma']*2.355)
    ## The fit is then plotted with a line and a scatterplot
    plt.plot(xVarList, resList)
    plt.scatter(xVarList, resList, marker = "x")
    ## Some if statements to allow x,y axes to be plotted with LaTex label
    ## Not super important.
    if xVar == "Qx":
        plt.xlabel(f"$Q_x~\AA^{-1}$")
    elif xVar == "Qy":
        plt.xlabel(f"$Q_y~\AA^{-1}$")
    elif xVar == "E":
        plt.xlabel(f"$E$ (meV)")
    else:
        plt.xlabel(f"{xVar}")
    if resVar == "Qx":
        plt.ylabel(f"$\Delta Q_x~\AA^{-1}$")
    elif resVar == "Qy":
        plt.ylabel(f"$\Delta Q_y~\AA^{-1}$")
    elif resVar == "E":
        plt.ylabel(f"$\Delta E$ (meV)")
    else:
        plt.ylabel(f"$\Delta${xVar}")
    if xlim != None:
        ## the spacing term adds a little breathing room from the plot
        ## so that the first and last x-y points are not the end of the plot
        spacing = (xlim[1]-xlim[0])*0.05
        plt.xlim(xlim[0]-spacing, xlim[1]+spacing)
    if ylim != None:
        plt.ylim(ylim[0], ylim[1])
    plt.title(f"{instrument.stations} Stations Mosaic {instrument.mosaic} {xVar} vs. {resVar} Resolution")
    if saveFile == True:
        plt.savefig(f"{instrument.pathBase}/{instrument.stations}_Stations_Mosaic_{instrument.mosaic}_{xVar}_v_{resVar}_Res_{dt.now().strftime('%Y_%m_%d_%H_%M_%S')}.pdf", format = "pdf")
    plt.show()    
    ## note that the resolution x and y values are returned, and thus can be stored if desired.
    ## It is also used in resolutionComp()
    return (xVarList, resList)


## The function below is in the event you want to compare different resolutions.
## You could easily do the function of this plot yourself using the output of resolution()
## But it is included for convenience. Essentially it takes in each instrument the user is comparing
## as a list [instrument1, instrument2,...] and the resolutions outputted from the resolution() function
## as a second variable resxy. So resxy = [(xres1, yres1,),.... ] where xres1, yres1 are the precalculated
## resolutions for instrument 1 as calculated from resoltuion(). This means that the resolution() function
## is already run manually by the user.
## Next, the xlim, ylim control axes (no special meaning), and saveFig, as always, will allow the user to
## save the figure as a pdf. Finally, the maskPoints is passed as a list of x-axis values
## you want removed. This is particularly useful for regions that are noisy, such as near
## q=0. Note that the closest value to the x-axis values you maskPoints is removed,
## so make sure you have an idea of what points you want to remove beforehand.
def resolutionComp(instrumentList, resxy, xVar, resVar, xlim = None, ylim = None, maskPoints = None, saveFig = False):
    for idx in range(len(instrumentList)):
        ## Essentially the x,y points for a given instrument are extracted.
        resx, resy = resxy[idx][0], resxy[idx][1]

        if maskPoints!= None:
            ## masked points are removed by finding the closest x-value index
            ## and removing it
            for point in maskPoints:
                ## index is finding the nearest x-point satisfying the 
                ## x value passed. Then the x, y value at that point 
                ## are removed.
                index = np.argmin(np.abs(np.array(resx) - point))
                resx.remove(resx[index])
                resy.remove(resy[index])
            ## Then the rest are plotted.
        plt.scatter(resx, resy, marker = "x")
        ## The instrument is passed to allow for the labeling in the legend.
        plt.plot(resx, resy, label = f'{instrumentList[idx].stations} Stations Mosaic {instrumentList[idx].mosaic}')
    ## again some conveniences included so that the x,y axis labels can have LaTex formatting.
    if xVar == "Qx":
        plt.xlabel(f"$Q_x~\AA^{-1}$")
    elif xVar == "Qy":
        plt.xlabel(f"$Q_y~\AA^{-1}$")
    elif xVar == "E":
        plt.xlabel(f"$E$ (meV)")
    else:
        plt.xlabel(f"{xVar}")
    if resVar == "Qx":
        plt.ylabel(f"$\Delta Q_x~\AA^{-1}$")
    elif resVar == "Qy":
        plt.ylabel(f"$\Delta Q_y~\AA^{-1}$")
    elif resVar == "E":
        plt.ylabel(f"$\Delta E$ (meV)")
    else:
        plt.ylabel(f"$\Delta${xVar}")
    plt.title(f"{xVar} vs. {resVar} Resolution")
    ## xlim ylim only control axes bounds.
    if xlim !=None:
        plt.xlim(xlim[0], xlim[1])
    if ylim != None:
        plt.ylim(ylim[0], ylim[1])
    if saveFig == True:
        plt.savefig(f"{instrument.pathBase}/{xVar}_v_{resVar}_Resolution_{dt.now().strftime('%Y_%m_%d_%H_%M_%S')}.pdf", format = "pdf")
    plt.legend()
    plt.show()
