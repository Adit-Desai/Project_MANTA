# Welcome!
# If you're trying to read through the code in this repository, it is recommended
# to read in the following order: 
# 1. Instrument_Creator.py 
# 2. Calibration.py
# 3. DataLoader.py
# 4. Plotting.py

## A general note: This is the most computationally intensive portion of the code.
## I've done some work to optimize it, but this is outside my region of expertise.
## Optimizations here would be the MOST valuable for future use to accompany a
## potential MANTA-like instrument.

##Various import statements used throughout the module
import numpy as np
import os
import pandas as pd
from tqdm import tqdm

## These functions are simple ways to calculate qx and qy from the experimental parameters
## based on the sample angle and twotheta (scattering angle)
## It would also be reasonable to modify this from qx, qy to reciprocal lattice units
## given a cif file. Or any basic structure. However this is not included yet as the
## simulations were conducted on a simple cubic structure.     
def qx_calculator(ki, kf, twoth, sampleAng):
    twothrad = np.deg2rad(twoth)
    sampleAngRad = np.deg2rad(sampleAng)
    Qx = ki*np.cos(-1*sampleAngRad) - kf*np.cos(-1*sampleAngRad + twothrad)
    return Qx
def qy_calculator(ki, kf, twoth, sampleAng):
    twothrad = np.deg2rad(twoth)
    sampleAngRad = np.deg2rad(sampleAng)
    Qy = ki*np.sin(-1*sampleAngRad) - kf*np.sin(-1*sampleAngRad + twothrad)
    return Qy

## This is the main function users will call on that accesses all their data
## dataLoader() requires the instrument object, the calibration from Calibration.py
## and the name of the folder where the data is located. Note the folder containing the data
## must be located in the same directory as the calibration data.
def dataLoader(instrument, calibration, folder):
    ## Access all datafiles there, any unwanted files currently have to be removed manually.
    allFiles = [f for f in os.listdir(f"{instrument.pathBase}/{folder}")]
    ## events will be the temporary list that will contain the Ei, Efs, intensities, sample angle,
    ## and twotheta of the simulation/data. It's the precursor to the pandas dataframe that will
    ## be created. It's just faster to work with lists/numpy arrays before building the pandas
    ## dataframe
    events = []
    
    ## Now I turn the calibration pandas dataframe into an array
    ## It is faster to turn it into essentially a matrix than to work
    ## with the dataframe.
    calibrationArr = np.array(calibration)
    
    
    ## Essentially toy models only have one angular channel while the full
    ## instrument has 8, this just reflects that for reading the data
    if instrument.type == "toy model":
        channelNum = 1
    elif instrument.type == "full":
        channelNum = 8
    ## This section sets up the tqdm progress bar (a convenience)
    ## So users can track how long their data will take to load
    progress = tqdm(allFiles)
    progress.set_description("FilesRead/TotalFiles")


    for file in progress:
        ## Go to the directory
        try:
            os.chdir(f"{instrument.pathBase}/{folder}/{file}")            
        except:
            print(f"{file} is not a folder!")
            continue
        ## This section is just to read the angle twotheta Ei, and sample Angle
        ## It is designed for McStas ReuterStokes.psd files. which also produce
        ## psdtube.dat files. ReuterStokes.psd files contain the actual data
        ## but aren't created unless a neutron lands in the specific tube.
        ## However psdtube.dat files, despite not containing the actual data
        ## are always created, which is why they're useful for extracting
        ## experimental parameters.

        ## I use this section to set up a true/false statement by having
        ## the experimental parameters initially undefined.
        twothBase = "undefined"
        Ei = "undefined"
        sampleAng = "undefined"
        ##Opening and reading the data
        fileOpener = open("psd_tube1_1a.dat", "r")
        fileData = fileOpener.readlines()
        ##This basically has it so it extracts the experimental parameters based on the structure of
        ## the psd_tube.dat files
        ## As soon as the three parameters, Ei, sampleAng, and twothBase are defined
        ## the for loop breaks. 
        for line in fileData:
            splitLine = line.split()
            if splitLine[1] == "Param:":
                paramSplit = splitLine[2].split("=")
                if paramSplit[0] == "Ei":
                    Ei = float(paramSplit[1])
                if paramSplit[0] == "TwoTh":
                    ## twothBase is the rotation of the entire angular detection system
                    ## So there are several angles within the CAMEA/MANTA subsystem that 
                    ## are increased by the constant term twoThBase
                    twothBase = float(paramSplit[1])
                if paramSplit[0] == "psi":
                    sampleAng = float(paramSplit[1])
                if twothBase != "undefined" and Ei != "undefined" and sampleAng != "undefined":
                    break
        ## This is just a catch in case something went wrong defining the parameters, instance
        ## has not yet occurred but could be useful for future debugging.
        if twothBase == "undefined" or Ei == "undefined" or sampleAng == "undefined":
            print(f"Something went wrong defining Ei, Two Theta, and the Sample Angle for File {file}")
            continue
        ## There are 8 angular channnels of detectors (controlled by channel and channelNum)
        ## there are 2 rows of detectors, the bottom row has 7 detectors and the top has 6 detectors
        ## the 6 on the top row are placed in between the 7, so each has a slightly different twoth
        for detRow in range(1, 3):
            if detRow == 1:
                detNum = 7
                ## there are the different angles within the bottom row of the tube relative
                ## to the center of the angular channel
                detAngList = np.array([-3.33, -2.22, -1.11, 0., 1.11, 2.22, 3.33])
                ## I think defining detAngList within the for loop is slightly less efficient
                ## but after tests it hardly seems to matter and it improves readability if its
                ## in the for loop. Future readers of code feel free to change that.
            elif detRow == 2:
                detNum = 6
                ## these are the different angles relative to the top row of the tube
                ## relative to the center of the angular channel
                detAngList = np.array([-2.775, -1.665, -0.555, 0.555, 1.665, 2.775])
            for channel in range(1, channelNum+1):
                for det in range(detNum):
                    try:
                        ##Opening the data based on the detector
                        fileOpener = open(f"ReuterStokes{str(detRow)}_{str(det+1)}_{channel}.psd")
                        fileData = fileOpener.readlines()
                    except FileNotFoundError:
                        continue
                    ## below is the true twotheta of the tube based on the twoThBase
                    ## each angular channel will be rotated by 7.5 degrees from the past one
                    ## and each tube will be somewhere in the middle of the 7.5. degree span
                    ## based on detAngList

                    ## For example, if twoThBase = 12, and we are looking at 5th detector tube 
                    ## on the bottom row of the third channel, twoTh is calculated as
                    ## twoth = 12 - 1.11 + 2*7.5 = 25.89 degrees
                    twoth = twothBase - detAngList[det] + ((channel-1) * 7.5)
                    ## next thing is creating a matrix of all the relevant parameters that will be needed
                    ## to calculate Q and E
                    ## This is the portion I think could be optimized further
                    fileEvents = np.zeros((5, len(calibration.index)))
                    fileEvents[0] += Ei
                    ## the indices are the different Efs used from the calibration

                    fileEvents[1] += np.array(calibration.index)
                    fileEvents[2] += twoth
                    fileEvents[3] += sampleAng
                    ## by the end of this, I have created a 5xN(Ef) matrix where each 
                    ## row corresponds to a specific parameter
                    ## The last row to be filled in is the intensity for each Ef

                    ## To do this, I make a list of all events and their corresponding intensities
                    ## This is similar to the way the calibration() in Calibration.py was run
                    ## Note that the list resets each new file as each one will have slightly
                    ## different experimental parameters.
                    fileyPos= []
                    fileIntensities = []
                    for line in fileData[instrument.startpoint():]:
                        splitLine = line.split()
                        ## the format of each event in the psd tube format is 
                        ## Intensity, x, y, z, ....
                        ## where the ... includes variables we don't care about
                        ## We need Intensity, and the y position, which is the vertical
                        ## location where the neutron lands on the detector
                        intensity = float(splitLine[0])
                        ## I ignore all 0 and "negative" intensity events (a weird quirk that shows up)
                        ## occasionally.
                        if intensity <= 0:
                            continue
                        else:
                            ## Now I keep track of the yposition and the corresponding intensity
                            ## By appending them to the relevant lists
                            ypos = float(splitLine[2])
                            fileyPos.append(ypos)
                            fileIntensities.append(intensity)
                    ## Now I histogram the data with the same bins as in the calibration.
                    ## The positions are passed as the "x" data and they are weighted by the 
                    ## intensities.
                    histogrammedData, binEdges = np.histogram(fileyPos, bins = 1024, weights = fileIntensities, range = (-0.45, 0.45))
                    
                    ## Now I matrix multiply. Essentially it multiplies an 1024 x N(Ef)
                    ## matrix by a 1024 column vector, turning it into an N(Ef) column
                    ## vector that has the intensities for each of the energies.
                    ## This is based off the prismatic weighting from the calibration.
                    updatedintensities = np.matmul(calibrationArr, histogrammedData)
                    ## I then add it to the fileEvents array and then append the final
                    ## fileEvents array to the events list.
                    fileEvents[4] += updatedintensities
                    ## I mentioned before that each row in fileEvents corresponds to a different parameter
                    ## It was easier to do this for the collection of data, but it will be easier to
                    ## perform the future calculations of Q using pandas, thus we transpose the matrix
                    ## to form the following row format: [Ei, Ef, twoth, sampleAng, Intensity]
                    events.append(fileEvents.transpose())
                      
    # Now that we have all the data, let's prepare it for the pandas dataframe                
    events = np.array(events)
    ## This step basically flattens the np array so that instead of having several matrices
    ## with the shape (N(Ef)*5) appended to each other, we get a single matrix with 5 columns
    ## with all the unique events being a different row
    events = np.reshape(events, newshape = (len(events)* len(events[0]), 5))
    ## now create the pandas dataframe
    data = pd.DataFrame(events, columns = ["Ei", "Ef", "Two Theta", "Sample Angle", "Intensity"])
    ## Now create the new columns used in plotting,
    data["E"] = (data["Ei"] - data["Ef"])
    data["ki"] = 2*np.pi/np.sqrt(81.8047/data["Ei"])
    data["kf"] = 2*np.pi/np.sqrt(81.8047/data["Ef"])
    ## Uses ki, kf, twoTheta, and sampleAng to calculate Qx and Qy
    ## Qx and Qy are calculated using the functions defined above in the code
    data["Qx"] = qx_calculator(data["ki"], data["kf"], data["Two Theta"], data["Sample Angle"])
    data["Qy"] = qy_calculator(data["ki"], data["kf"], data["Two Theta"], data["Sample Angle"])
    data["Intensity"] = data["Intensity"] * data["ki"]/data["kf"]
    return data
