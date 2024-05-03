# Welcome!
# If you're trying to read through the code in this repository, it is recommended
# to read in the following order: 
# 1. Instrument_Creator.py 
# 2. Calibration.py
# 3. DataLoader.py
# 4. Plotting.py

## The Instrument File contains useful information about the different designs
## and is used to assist in all the other files in some fashion, though it is 
## most valuable in the calibration file. It will essentially create an instrument object
## which is passed through the other files and used to help it make decisions.


##Here are the necessary import statements for this file
import numpy as np
import os
import sys

##Defining the basic Instrument class
class Instrument:
    ## the stations number (int), mosaic (int) in arc minutes are required inputs.
    ## The CAMEA design (and the one used for the main MANTA simulations)
    ## have 8 stations and an analyzer mosaic of 60'.

    def __init__(self, stations, mosaic, pathBase = "", type = "full"):
        ## The pathBase specifies the location where the calibration data and the experimental data is located.
        ## Note that the calibration data and experimental data must be in separate folders located in the pathBase
        ## directory. If pathBase is not specified, then the instrument will assume the data is located in the same
        ## directory as the repository.
        if pathBase == "":
            self.pathBase = os.path.dirname(os.path.abspath(sys.argv[0]))
        else:
            self.pathBase = pathBase
        self.stations = stations
        self.mosaic = mosaic
        ## Here you can specify if the instrument is full or toy model. By default, the instrument will assume you look at the full
        if type == "full" or type == "FULL" or type == "MANTA":
            self.type = "full"
        elif type == "toy model" or type == "Toy Model" or type == "Toy_Model" or type == "toy_model":
            self.type = "toy model"
        else:
            print("Instrument type not recognized! Please specify as 'full' or 'toy model'.")
    def stationList(self):
        ## these are the Bragg energies of each of the analyzers for each design in meV
        if self.stations == 8:
            return [3.21, 3.38, 3.58, 3.8, 4.05, 4.33, 4.64, 5.01]
        elif self.stations == 5:
            return [3.2, 3.6, 4.0, 4.4, 4.8]
        elif self.stations == 10:
            return [3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0]
        else:
            print("Something went wrong defining the Bragg energies of the analyzers! Please pick 5, 8, or 10 as your number of stations.")
    def energyList(self):
        ## energyList will return all the energies used for prismatic analysis
        ## Currently its setup for only one specific type of calibration per instrument
        ## but modifying this section would allow customization of calibration procedure
        myList = []
        for baseEnergy in self.stationList():
            ## I append different energies to myList, all delE != 0 energies
            ## are analyzed through the prismatic analysis procedure. Refer to
            ## the paper for more details. 
            ## As a starting point the Bragg energy is used along with
            ## additonal prismatic energies using delE, where 
            ## -0.08 meV <= delE <= 0.12 meV. 
            for delE in np.arange(-0.08, 0.13, 0.04):
                ##the round() function is just to deal with float instability
                energy = round(baseEnergy + delE, 3)
                if energy in myList:
                    continue
                else:
                    myList.append(energy)
        ## Occasionally more energies are needed for correct prismatic analysis, particularly for larger mosaic.
        ## Prismatic analysis will essentially "stretch" neutron events from a point to a line
        ## on the energy axis. If there are not enough energies, then artifacts that look like 
        ## the end of the line of prismatic energies seem unusually intense, then additional calibration
        ## energies are needed.
        if self.stations == 5 and self.mosaic == 60:
            myList += [3.08, 3.48, 3.88, 4.28, 4.56, 4.68]
        elif self.stations == 5 and self.mosaic == 120:
            myList = []
            ##for 5 station mosaic 120, enough addtional energies are needed where the basis of the calibration range
            ## (previously -0.08 meV <= delE <= 0.12, meV) effectively enlargens to -0.20 meV <= delE <= 0.20 meV  
            for baseEnergy in self.stationList():
                for delE in np.arange(-0.2, 0.21, 0.04):
                    energy = round(baseEnergy + delE, 3)
                    if energy in myList or energy == 3.0 or energy == 3.04 or energy == 3.08: 
                        continue
                    else:
                        myList.append(energy)
                    myList += [5.04]
        elif self.stations == 8 and self.mosaic == 60:
            myList += [4.52, 4.89]
        elif self.stations == 8 and self.mosaic == 120:
            myList += [4.84, 4.88, 4.49, 4.21, 5.17, 5.21, 5.25, 4.53, 4.8]
        elif self.stations == 10 and self.mosaic == 120:
            myList += [5.16, 5.20, 5.24]
        myList.sort()
        return myList
    ## The baffleRegions are used for the calibration. Essentially, the neutron analyzers are focused such that 
    ## each analyzer is focused to a specific "strip" of neutrons centered at
    ## a discrete position on the detector with a finite bandwidth. In a real experiment,
    ## baffles will be placed such that there is no crosstalk between analyzers, and essentially regions 
    ## strips will have no intensity. The purpose of the baffleRegions() method is to mimic that effect.
    ## It also helps deal with stray neutron events that were scattered into regions between the strips,
    ## where there is significantly larger uncertainty.

    ## The baffle Regions will vary based on number of stations and mosaic
    ## So currently each of them are predefined here. If a new design is used
    ## you'll have to specify custom baffleRegions by looking at the calibration plots.
    def baffleRegions(self):
        if self.stations == 5 and self.mosaic == 30:
            baffleRegions = [(0.115-0.45, 0.265-0.45), (0.315-0.45, 0.45-0.45),
            (0.5-0.45, 0.62-0.45), (0.675-0.45, 0.79-0.45)]
        elif self.stations == 5 and self.mosaic == 60:
            baffleRegions = [(0.12-0.45, 0.255-0.45), (0.32-0.45, 0.44-0.45),
            (0.505-0.45, 0.615-0.45), (0.675-0.45, 0.785-0.45)]
        elif self.stations == 5 and self.mosaic == 120:
            baffleRegions = [(0.145-0.45, 0.23-0.45), (0.34-0.45, 0.425-0.45),
            (0.525-0.45, 0.605-0.45), (0.705-0.45, 0.77-0.45)]
        elif self.stations == 8 and self.mosaic == 30:
            baffleRegions = [(0.115-0.45, 0.165-0.45), (0.22-0.45, 0.27-0.45),
            (0.325-0.45, 0.37-0.45), (0.435-0.45, 0.48-0.45), 
            (0.53-0.45, 0.585-0.45),(0.635-0.45, 0.69 -0.45),
            (0.745-0.45, 0.805-0.45)]
        elif self.stations == 8 and self.mosaic == 60:
            baffleRegions = [(0.115-0.45, 0.165-0.45), (0.22-0.45, 0.27-0.45),
            (0.325-0.45, 0.37-0.45), (0.435-0.45, 0.47-0.45), 
            (0.535-0.45, 0.58-0.45),(0.64-0.45, 0.685 -0.45),
            (0.75-0.45, 0.8-0.45)]
        elif self.stations == 8 and self.mosaic == 120:
            baffleRegions = [(0.115-0.45, 0.17-0.45), (0.22-0.45, 0.275-0.45),
            (0.33-0.45, 0.375-0.45), (0.435-0.45, 0.48-0.45), 
            (0.535-0.45, 0.59-0.45),(0.64-0.45, 0.69 -0.45),
            (0.75-0.45, 0.81-0.45)]
        elif self.stations == 10 and self.mosaic == 30:
            baffleRegions = [(0.07-0.45, 0.125-0.45), (0.175-0.45, 0.22-0.45),
            (0.275-0.45, 0.315-0.45), (0.36-0.45, 0.405-0.45), 
            (0.45-0.45, 0.49-0.45),(0.535-0.45, 0.575 -0.45),
            (0.62-0.45, 0.655-0.45), (0.70-0.45, 0.735-0.45),
                (0.785-0.45, 0.815-0.45)]
        elif self.stations == 10 and self.mosaic == 60:
            baffleRegions = [(0.075-0.45, 0.125-0.45), (0.175-0.45, 0.22-0.45),
            (0.275-0.45, 0.305-0.45), (0.37-0.45, 0.405-0.45), 
            (0.46-0.45, 0.49-0.45),(0.55-0.45, 0.575 -0.45),
            (0.625-0.45, 0.655-0.45), (0.71-0.45, 0.74-0.45),
                (0.79-0.45, 0.815-0.45)]
        elif self.stations == 10 and self.mosaic == 120:
            baffleRegions = [(0.065-0.45, 0.13-0.45), (0.175-0.45, 0.235-0.45),
            (0.27-0.45, 0.33-0.45), (0.365-0.45, 0.42-0.45), 
            (0.45-0.45, 0.505-0.45),(0.535-0.45, 0.59 -0.45),
            (0.62-0.45, 0.67-0.45), (0.705-0.45, 0.76-0.45),
                (0.78-0.45, 0.835-0.45)]
        return baffleRegions
    ## The startpoint is a minutia detail that essentially says when in McStas Reuter_Stokes_.psd 
    ## data files (used in Calibration.py and DataLoader.py) the data entries begin. It changes based 
    ## on how you write the McStas file. 
    def startpoint(self):
        if self.stations == 10:
            index = 49
        elif self.stations == 5 or self.stations == 8:
            index = 46
        return index
