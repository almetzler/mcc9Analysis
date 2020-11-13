import ROOT, math, sys, os
import uproot
from ROOT import TH1, TAxis, gROOT, TCanvas

import pandas as pd
import numpy as np
import json

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

'''
with open(filepath, 'w') as f:
    for chunk in json.JSONEncoder().iterencode(object_to_encode):
        f.write(chunk)

'''        

# eventsExtJSON = eventsExt.to_json(orient = 'split')
# with open("Data/eventsExt",'w') as fle:
#     for chunk in json.JSONEncoder().iterencode(eventsExtJSON):
#         fle.write(json.dumps(chunk))
# print('eventExt done')

# eventsOnBeamJSON = eventsOnBeam.to_json(orient='split')
# with open("Data/eventsOnBeam",'w') as fle: 
#     for chunk in json.JSONEncoder().iterencode(eventsOnBeamJSON):
#         fle.write(json.dumps(chunk))
# print('eventsOnBeam done')

# eventsOverlayJSON = eventsOverlay.to_json(orient='split')
# with open("Data/eventsOverlay",'w') as fle:
#     for chunk in json.JSONEncoder().iterencode(eventsOverlayJSON):
#         fle.write(json.dumps(chunk))
# print('eventsOverlay done')

# eventsDirtJSON = eventsDirt.to_json(orient='split')
# with open("Data/eventsDirt",'w') as fle:
#     for chunk in json.JSONEncoder().iterencode(eventsDirtJSON):
#         fle.write(json.dumps(chunk))
# print('eventsDirt done')

# eventsExtJSON = eventsExt.to_csv("Data/eventsExt.csv")
print('eventExt done')

# eventsOnBeamJSON = eventsOnBeam.to_csv("Data/eventsOnBeam.csv")
print('eventsOnBeam done')

eventsOverlayJSON = eventsOverlay.to_csv("Data/eventsOverlay.csv")
print('eventsOverlay done')

eventsDirtJSON = eventsDirt.to_csv("Data/eventsDirt.csv")
print('eventsDirt done')