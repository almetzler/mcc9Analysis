import ROOT, math, sys, os
import uproot
from ROOT import TH1, TAxis, gROOT, TCanvas

import pandas as pd
import numpy as np
import json

#Custom
from PlotUtils import PlotUtils

####################################################################################################

#EventTypes = {"1 Pi0 0 PiP" : 0, "1 Pi0 N PiP" : 1, "CC 0 Pi0" : 2, "CC N Pi0" : 3, "CC Other" : 4, "NC 0 Pi0" : 5, "NC N Pi0" : 6, "NC Other" : 7, "OOFV" : 8, "Non numu" : 9}
IntNames = ("1 Pi0 0 PiP", "1 Pi0 N PiP", "CC 0 Pi0", "CC N Pi0", "CC Other", "NC 0 Pi0", "NC N Pi0", "NC Other", "OOFV", "Non numu")

Plotter = PlotUtils(IntNames, "mc_EventType")


def getFlashWgt(x):
  returnVal = expDegSix(x)
  if(returnVal < 0):
    return 1.0
  else:   
    return returnVal

def expDegSix(x):
   p = [ 1.99384917e-01, -1.05022405e-10,  6.34916308e-05, -3.44705817e-03,
  7.32059307e-02, -5.91696006e-01,  2.14667463e+00, -1.02545380e+00,
  3.65915734e-01]
   return(np.exp(-p[0]*x)*(p[1]*pow(x, 6) + p[2]*pow(x,5) + p[3]*pow(x,4) + p[4]*pow(x,3) + p[5]*pow(x,2) + p[6]*x + p[7]) + p[8])  
  
def LoadShowerWeights(InputFrame):
  weights = InputFrame['wgt']
  showerL = InputFrame['_fShowerEnergy2'].str.len()#This is a dummy variable to get the number of the showers
  InputFrame.insert(InputFrame.shape[1], "ShowerWeights", [np.full(N, wgt) for wgt, N in zip(weights.to_numpy(dtype=object), showerL.to_numpy(dtype=object) ) ] ) 

def LoadTrackWeights(InputFrame):
  weights = InputFrame['wgt']
  showerL = InputFrame['_fTrackPID'].str.len()#This is a dummy variable to get the number of tracks
  InputFrame.insert(InputFrame.shape[1], "TrackWeights", [np.full(N, wgt) for wgt, N in zip(weights.to_numpy(dtype=object), showerL.to_numpy(dtype=object) ) ] ) 

def DefineEventType(isCC, isFV, nuPDG, nPi0, nPiP, nMuon):
    if(not isFV):
       return Plotter.eventDict["OOFV"]
    if(nuPDG != 14):
       return Plotter.eventDict["Non numu"]
    
    if(isCC == 0):
      if(nPi0 == 1 and nPiP == 0 and nMuon == 1):
        return Plotter.eventDict["1 Pi0 0 PiP"]
      elif(nPi0 == 1 and nPiP > 0):
        return Plotter.eventDict["1 Pi0 N PiP"]
      elif(nPi0 == 0):
        return Plotter.eventDict["CC 0 Pi0"]
      elif(nPi0 > 1): 
        return Plotter.eventDict["CC N Pi0"]
      else:
        return Plotter.eventDict["CC Other"]
    
    else:
      if(nPi0 == 0):
        return Plotter.eventDict["NC 0 Pi0"]
      elif(nPi0 > 1): 
        return Plotter.eventDict["NC N Pi0"]
      else:
        return Plotter.eventDict["NC Other"]                  

def LoadEventTypes(inputFrame):
    inputFrame.insert(inputFrame.shape[1], "mc_EventType", [DefineEventType(i, j, k, x, y, z) for i,  j, k, x, y, z in zip(inputFrame['_fNuCCNC'], inputFrame['_fNuInFV'], inputFrame['_fNuPDG'], inputFrame['_fNpi0'], inputFrame['_fNpiplus'],  inputFrame['_fNmuon'])] )

InputFiles  = ["/uboone/data/users/joelam/CCPi0Ntuples/Run1Overlay.root", "/uboone/data/users/joelam/CCPi0Ntuples/RunC1Ext.root", "/uboone/data/users/joelam/CCPi0Ntuples/Run1Data.root", "/uboone/data/users/joelam/CCPi0Ntuples/Run1Dirt.root"]
with open("SelectionVariables.txt") as file:
  selectionVariables = [line.strip() for line in file if not line.startswith('#')]

with open("PlottingVariables.txt") as file:
  plottingVariables = [line.strip() for line in file if not line.startswith('#')]

selectionVariables.extend(plottingVariables)

ExtScale         = 1.02767
extTriggersC1    = 33630174.0
dataPOT          =  1.547e+20
bnbSpills        = 34361582.0
signalMassLow    = 60.0
signalMassHigh   = 180.0

#Python library to read in ROOT ntuples files.
bnbEvents        = uproot.open(InputFiles[1])["efficiency/eventtree"]
eventsExt    = pd.DataFrame(bnbEvents.arrays(selectionVariables) )

bnbEvents        = uproot.open(InputFiles[2])["efficiency/eventtree"]
eventsOnBeam    = pd.DataFrame(bnbEvents.arrays(selectionVariables) )

selectionVariables.append("_fCVWeight")

bnbEvents        = uproot.open(InputFiles[0])["efficiency/eventtree"]
eventsOverlay    = pd.DataFrame(bnbEvents.arrays(selectionVariables) )

selectionVariables.append("_fNuEnergy")

bnbEvents        = uproot.open(InputFiles[3])["efficiency/eventtree"]
eventsDirt    = pd.DataFrame(bnbEvents.arrays(selectionVariables) )

overlayPOT    = uproot.open(InputFiles[0])["efficiency/tree"]
dirtPOT       = uproot.open(InputFiles[3])["efficiency/tree"]

sumOverlayPOT = (pd.Series(overlayPOT.array("pot"))).sum()
sumDirtPOT    = (pd.Series(dirtPOT.array("pot"))).sum()

overlayWeights = np.full(eventsOverlay.shape[0], dataPOT / sumOverlayPOT )
dirtWeights    = np.full(eventsDirt.shape[0], dataPOT / sumDirtPOT)

eventsOverlay.insert(eventsOverlay.shape[1], "pot_wgt", overlayWeights )
eventsDirt.insert(eventsDirt.shape[1], "pot_wgt", dirtWeights )

eventsOverlay.insert(eventsOverlay.shape[1], "flash_wgt", [getFlashWgt(x) for x in eventsOverlay['_fFlashChi2']])

eventsOverlay.eval('wgt = pot_wgt*_fCVWeight*flash_wgt', inplace=True)
eventsDirt.eval('wgt = pot_wgt*_fCVWeight', inplace=True)   

extWeights              = np.full(eventsExt.shape[0],  (bnbSpills / extTriggersC1) )
eventsExt.insert(eventsExt.shape[1], "wgt", extWeights)

eventsWithShowers      = eventsOverlay.query('_fNLeadCandidates > 0 and _fNSubleadCandidates > 0 and _fHasCandidateNeutrino == 1 and _fHasCandidateMuon == 1')
dirtEventsWithShowers  = eventsDirt.query('_fNLeadCandidates > 0 and _fNSubleadCandidates > 0 and _fHasCandidateNeutrino == 1 and _fHasCandidateMuon == 1')
extEventsWithShowers   = eventsExt.query('_fNLeadCandidates > 0 and _fNSubleadCandidates > 0 and _fHasCandidateNeutrino == 1 and _fHasCandidateMuon == 1')
dataEventsWithShowers  = eventsOnBeam.query('_fNLeadCandidates > 0 and _fNSubleadCandidates > 0 and _fHasCandidateNeutrino == 1 and _fHasCandidateMuon == 1')

LoadEventTypes(eventsWithShowers)

LoadShowerWeights(eventsWithShowers)
LoadTrackWeights(eventsWithShowers)
LoadShowerWeights(dirtEventsWithShowers)
LoadTrackWeights(dirtEventsWithShowers)
LoadShowerWeights(extEventsWithShowers)
LoadTrackWeights(extEventsWithShowers)

Plotter.setWeightName('wgt')
singalEventsOverlay = eventsWithShowers.query('_fNChargedPiCandidates == 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates == 1 and _fNOtherFancyPairs == 0')
singalEventsDirt    = dirtEventsWithShowers.query('_fNChargedPiCandidates == 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates == 1 and _fNOtherFancyPairs == 0')
singalEventsExt     = extEventsWithShowers.query('_fNChargedPiCandidates == 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates == 1 and _fNOtherFancyPairs == 0')
signalEventsData    = dataEventsWithShowers.query('_fNChargedPiCandidates == 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates == 1 and _fNOtherFancyPairs == 0')

InvariantMassRange = (0, 500)
axesLabels = ["Pi0 Invariant Mass (Two MIP Sideband)", "Invaraint Mass (MeV/c2)", "Number of Events"]
limits = {"xlimits" : (200, ), "Titles" : axesLabels }

#2 Mip+
twoMIPEventsOverlay    = eventsWithShowers.query('_fNChargedPiCandidates > 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates == 1 and _fNOtherFancyPairs == 0')
twoMIPEventsDirt       = dirtEventsWithShowers.query('_fNChargedPiCandidates > 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates == 1 and _fNOtherFancyPairs == 0')
twoMIPEventsExt        = extEventsWithShowers.query('_fNChargedPiCandidates > 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates == 1 and _fNOtherFancyPairs == 0')
twoMIPEventsData       = dataEventsWithShowers.query('_fNChargedPiCandidates > 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates == 1 and _fNOtherFancyPairs == 0')

'''
data types that have given me a memory error :/
- json
- hdf5

other data types
- csv with chunk size
  - didnt serialize/deserialize well
  - to_csv("Data/eventsExt.csv", chunksize=1, float_format='%g', na_rep=np.nan)

pickle worked :)
'''        

eventsExtMIP = twoMIPEventsExt.to_pickle("Data/eventsExtMIP.pkl")
eventsExtSig = singalEventsExt.to_pickle("Data/eventsExtSig.pkl")
# eventsExt.to_csv("Data/eventsExt.csv", chunksize=1, sep='|')
# eventsExt.to_hdf("Data/eventsExt.h5", key = 'df', mode='w', encoding = 'UTF-8')
# print(eventsExt.shape)
print('eventExt done')

eventsOverlayMIP = twoMIPEventsOverlay.to_pickle("Data/eventsOverlayMIP.pkl")
eventsOverlaySig = singalEventsOverlay.to_pickle("Data/eventsOverlaySig.pkl")
# eventsOverlay.to_csv("Data/eventsOverlay.csv", chunksize=1, sep='|')
# eventsOverlay.to_hdf("Data/eventsOverlay.h5", key = 'df', mode='w', encoding = 'UTF-8')
# print(eventsOverlay.shape)
print('eventsOverlay done')

eventsOnBeamMIP = twoMIPEventsData.to_pickle("Data/eventsOnBeamMIP.pkl")
eventsOnBeamSig = signalEventsData.to_pickle("Data/eventsOnBeamSig.pkl")
# eventsOnBeam.to_csv("Data/eventsOnBeam.csv", chunksize=1, sep='|')
# eventsOnBeam.to_hdf("Data/eventsOnBeam.h5", key = 'df', mode='w', encoding = 'UTF-8')
# print(eventsOnBeam.shape)
print('eventsOnBeam done')

eventsDirtMIP = twoMIPEventsDirt.to_pickle("Data/eventsDirtMIP.pkl")
eventsDirtSig = singalEventsDirt.to_pickle("Data/eventsDirtSig.pkl")
# eventsDirt.to_csv("Data/eventsDirt.csv", chunksize=1, sep='|')
# eventsDirt.to_hdf("Data/eventsDirt.h5", key = 'df', mode='w', encoding = 'UTF-8')
# print(eventsDirt.shape)
print('eventsDirt done')