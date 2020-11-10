import ROOT, math, sys, os
import uproot
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import csv
import re
import functools
import scipy as sci
from mpl_toolkits.mplot3d import Axes3D
from ROOT import TH1, TAxis, gROOT, TCanvas
from scipy import stats

#Custom
from PlotUtils import PlotUtils

####################################################################################################

#EventTypes = {"1 Pi0 0 PiP" : 0, "1 Pi0 N PiP" : 1, "CC 0 Pi0" : 2, "CC N Pi0" : 3, "CC Other" : 4, "NC 0 Pi0" : 5, "NC N Pi0" : 6, "NC Other" : 7, "OOFV" : 8, "Non numu" : 9}
IntNames = ("1 Pi0 0 PiP", "1 Pi0 N PiP", "CC 0 Pi0", "CC N Pi0", "CC Other", "NC 0 Pi0", "NC N Pi0", "NC Other", "OOFV", "Non numu")

Plotter = PlotUtils(IntNames, "mc_EventType")
 
def passesE(energy):
   minYPlaneEnergy   = 50.0
   return np.greater(energy, minYPlaneEnergy)

def passesESubLeading(energy):
   minYPlaneEnergy = 25.0
   return np.greater(energy, minYPlaneEnergy)

def passesRadialAngle(radialAng):
   minRadialAngle = 0.96
   return np.greater(radialAng, minRadialAngle)

def passesConversionD(dist3d, energy):
   minConversionDist = 6.0 
   maxConversionDist = 82.0
   minYPlaneEnergy   = 50.0
   maxYPlaneEnergy   = 300

   return( (np.greater(dist3d, minConversionDist) & np.less(dist3d, maxConversionDist) ) | ( np.greater(energy, minYPlaneEnergy) & np.less(energy, maxYPlaneEnergy) & np.less(dist3d, minConversionDist) ) ) 

def passesConversionDSubLeading(dist3d):
   minConversionDist = 1.0
   return np.greater(dist3d, minConversionDist)

def angleBetweenTwoShowers(shower1, shower2):
  return( shower1[0]*shower2[0] + shower1[1]*shower2[1] + shower1[2]*shower2[2] )

def LoadShowerWeights(InputFrame):
  weights = InputFrame['wgt']
  showerL = InputFrame['_fShowerEnergy2'].str.len()#This is a dummy variable to get the number of the showers
  InputFrame.insert(InputFrame.shape[1], "ShowerWeights", [np.full(N, wgt) for wgt, N in zip(weights.to_numpy(dtype=object), showerL.to_numpy(dtype=object) ) ] ) 


def LoadTrackWeights(InputFrame):
  weights = InputFrame['wgt']
  showerL = InputFrame['_fTrackPID'].str.len()#This is a dummy variable to get the number of tracks
  InputFrame.insert(InputFrame.shape[1], "TrackWeights", [np.full(N, wgt) for wgt, N in zip(weights.to_numpy(dtype=object), showerL.to_numpy(dtype=object) ) ] ) 

def TagLeadingShower(InputFrame):
  minRadialAngle    = 0.96
  signalMassLow    = 60.0
  signalMassHigh   = 180.0
  showerE = InputFrame['_fShowerEnergy2']
  showerDist = InputFrame['_fShowerDist3d']
  showerAngle = InputFrame['_fShowerAng3d']
  leadingShowers = [np.argmax(row) for row in showerE.to_numpy(dtype=object) ]
  passesEnergy   = [passesE(row) for row in showerE.to_numpy(dtype=object)]
  passesSubLeadingE = [passesESubLeading(row) for row in showerE.to_numpy(dtype=object)]
  passesAngle    = [passesRadialAngle(row) for row in showerAngle.to_numpy(dtype=object)]
  passesConerstionDist = [passesConversionD(rowA, rowB) for rowA, rowB in zip(showerDist.to_numpy(dtype=object), showerE.to_numpy(dtype=object))]
  passesSubLeadingConverstionDist = [passesConversionDSubLeading(rowA) for rowA in showerDist.to_numpy(dtype=object) ]
  InputFrame.insert(InputFrame.shape[1], "LeadingShowerIdx", leadingShowers)
  InputFrame.insert(InputFrame.shape[1], "passesEnergy", passesEnergy)
  InputFrame.insert(InputFrame.shape[1], "passesAngle", passesAngle)
  InputFrame.insert(InputFrame.shape[1], "passesConvD", passesConerstionDist)
  InputFrame.insert(InputFrame.shape[1], "passesSubLeadingEnergy", passesSubLeadingE)
  InputFrame.insert(InputFrame.shape[1], "passesSubLeadingConvD", passesSubLeadingConverstionDist)
  InputFrame.insert(InputFrame.shape[1], "NumPassConv", [sum(l) for l in  InputFrame['passesSubLeadingConvD'] ] )
  InputFrame.insert(InputFrame.shape[1], "NumPassE", [sum(l) for l in InputFrame['passesSubLeadingEnergy'] ] )
  InputFrame.insert(InputFrame.shape[1], "PassMass", [map(lambda l: l < signalMassHigh and l > signalMassLow, x) for x in  InputFrame['_fMultiPairPi0Mass'] ] )
  InputFrame.insert(InputFrame.shape[1], "NumPassMass", [sum(l) for l in  InputFrame['PassMass'] ] )

def TagSubLeadingShower(InputFrame):
  #calculate invariant mass for all pairs, then apply down seleciton later
  #As long as you have 2 pairs, you have an invariant mass one way or the other
  #And you can always apply shower selection / masking down the road
  #Get the mass for each leading / other shower pair. This may lead to multiple masses per event.
  #This is ok as in most cases we can just pick the mass that is closest to the pi0 mass
  #Ask / look up what happens if we have 2 track pairs and both are consistent w pi0 mass
  signalMassLow    = 60.0
  signalMassHigh   = 180.0
  subLeadingIdx    = 0 #Solve for this
  leadingShowers    = np.array( [shower[idx]  for idx, shower in zip(InputFrame["LeadingShowerIdx"].to_list(), InputFrame['_fShowerEnergy2'].to_list()) ] )
  subLeadingShowers = [np.delete(shower, idx) for idx, shower in zip(InputFrame["LeadingShowerIdx"].to_list(), InputFrame['_fShowerEnergy2'].to_list() ) ]
 
def TagNaNMIPs(InputFrame):
  '''
  The two MIP sideband has two problems:
  1) There are nan log(chi2) values, presumably from when the ratio is 0 that are tagged as two MIP events
  2) Sometimes, the primary muon fails the log(chi2) < -0.9. The upstream C code tags theve events correctly, but trying to re-do it in python will show 1 MIP events in the 2 MIP sample
  In case 1, we should remove the events. There is no good reason to believe a log(chi2) of nan (ie chi2 <= 0) is pion like
  In case 2, we should keep the events. As long as there is one other MIP like particle, the upstream C code counts the primary muon regardless of it's log(chi2) 
  '''
  maxMIPPID = -0.9
  InputFrame.insert(InputFrame.shape[1], "isMIP", [map(lambda p : p <= maxMIPPID, x) for x in InputFrame['_fTrackPID'] ] )
  InputFrame.insert(InputFrame.shape[1], "isNanMIP", [map(lambda p : math.isnan(p), x) for x in InputFrame['_fTrackPID'] ] )
  InputFrame.insert(InputFrame.shape[1], "primMuonIsMIP", [(lambda p : p <= maxMIPPID)(x) for x in InputFrame['_fCandidateMuonPID'] ]  )
  InputFrame.insert(InputFrame.shape[1], "NumMIP", [sum(l) for l in InputFrame['isMIP'] ] )
  InputFrame.insert(InputFrame.shape[1], "NumNanMIP", [sum(l) for l in InputFrame['isNanMIP'] ] )

def PlotLeadingShowerVarialbes(dataFrame):
  varList = ("_fShowerEnergy2", "_fShowerEnergy0")
  leadingShowers = dataFrame["LeadingShowerIdx"]
  #Create a mask of passing shower cuts
  passesEnergy = [shower[idx] for idx, shower in zip(leadingShowers, dataFrame["passesEnergy"].to_list() ) ]
  passesAngle  = [shower[idx] for idx, shower in zip(leadingShowers, dataFrame["passesAngle"].to_list() ) ]
  passesConvD  = [shower[idx] for idx, shower in zip(leadingShowers, dataFrame["passesConvD"].to_list() ) ]
  passesCuts   = np.array(passesEnergy and passesAngle and passesConvD)
  for var in varList:
    
    showerSeries = dataFrame[var].tolist()
    plotVar = np.array([shower[idx] for idx, shower in zip(leadingShowers, showerSeries) ] )
    plotVar = plotVar[passesCuts]

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


#InputFiles = ["/uboone/data/users/suprajab/MCC9/March2020/CCpi0Trees/bnb_eventtree.root", "/uboone/data/users/suprajab/MCC9/March2020/CCpi0Trees/ext_eventtree.root", "/uboone/data/users/suprajab/MCC9/March2020/CCpi0Trees/dirt_eventtree.root"]
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

#SAVE THIS
#print trackOverlay.loc[:100,['isContained', 'track_range_mom_mu', 'track_mcs_mom', 'track_mom_best']]

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

#Tag the leading showers
#TagLeadingShower(eventsWithShowers)
#TagSubLeadingShower(eventsWithShowers)
#Build up frames of the signal events and various sidebands
Plotter.setWeightName('wgt')
singalEventsOverlay = eventsWithShowers.query('_fNChargedPiCandidates == 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates == 1 and _fNOtherFancyPairs == 0')
singalEventsDirt    = dirtEventsWithShowers.query('_fNChargedPiCandidates == 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates == 1 and _fNOtherFancyPairs == 0')
singalEventsExt     = extEventsWithShowers.query('_fNChargedPiCandidates == 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates == 1 and _fNOtherFancyPairs == 0')
signalEventsData    = dataEventsWithShowers.query('_fNChargedPiCandidates == 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates == 1 and _fNOtherFancyPairs == 0')



InvariantMassRange = (0, 500)
axesLabels = ["Pi0 Invariant Mass (Two MIP Sideband)", "Invaraint Mass (MeV/c2)", "Number of Events"]
limits = {"xlimits" : (200, ), "Titles" : axesLabels }
Plotter.makeDataMCStack(singalEventsOverlay, singalEventsDirt, singalEventsExt, signalEventsData, '_fCandidatePi0Mass',(0, 500), 25, "SignalPi0Mass", {}, limits)


#2 Mip+
twoMIPEventsOverlay    = eventsWithShowers.query('_fNChargedPiCandidates > 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates == 1 and _fNOtherFancyPairs == 0')
twoMIPEventsDirt       = dirtEventsWithShowers.query('_fNChargedPiCandidates > 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates == 1 and _fNOtherFancyPairs == 0')
twoMIPEventsExt        = extEventsWithShowers.query('_fNChargedPiCandidates > 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates == 1 and _fNOtherFancyPairs == 0')
twoMIPEventsData       = dataEventsWithShowers.query('_fNChargedPiCandidates > 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates == 1 and _fNOtherFancyPairs == 0')

TagNaNMIPs(twoMIPEventsOverlay)
TagNaNMIPs(twoMIPEventsDirt)
TagNaNMIPs(twoMIPEventsExt)
TagNaNMIPs(twoMIPEventsData)


Plotter.makeDataMCStack(twoMIPEventsOverlay.query('NumNanMIP == 0'), twoMIPEventsDirt.query('NumNanMIP == 0'), twoMIPEventsExt.query('NumNanMIP == 0'), twoMIPEventsData.query('NumNanMIP == 0'), '_fTwoMIPPi0Mass', InvariantMassRange, 25, "TwoMIPPi0Mass", {}, limits)
limits["xlimits"] = ()
limits["Titles"]  = ["Number of MIPs (Two MIP Sideband)", "N MIPs", "Number of Events"]
NumPairsRange = (0,10)
Plotter.makeDataMCStack(twoMIPEventsOverlay.query('NumNanMIP == 0'), twoMIPEventsDirt.query('NumNanMIP == 0'), twoMIPEventsExt.query('NumNanMIP == 0'), twoMIPEventsData.query('NumNanMIP == 0'), 'NumMIP', NumPairsRange, 10, "TwoMIPNumMIPs", {}, limits)


#twoMIPEvents = twoMIPEvents['_fTrackPID'].fillna(9999.9)


#print twoMIPEvents.query('NumMIP == 1 and NumNanMIP == 0 and primMuonIsMIP == True' ).loc[:100000, ['_fTrackPID', 'NumNanMIP', '_fCandidateMuonPID', 'primMuonIsMIP']]
#raw_input()
#.stack().reset_index(drop=True).to_numpy()
#Confirm all 1 MIP events are coming from nans
#Then set NumMIP > 1 in query



Stack = Plotter.SmartStack(twoMIPEventsOverlay.query('NumNanMIP == 0'), 'NumMIP')
Plotter.makeMCStackedHistogram(Stack[0], Stack[1], (0,10), 10, "TwoMIPNumMIPs", legend=Stack[2])

#multi pair
multiPairEventsOverlay = eventsWithShowers.query('_fNChargedPiCandidates == 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates > 1 and _fNOtherFancyPairs == 0')
multiPairEventsDirt    = dirtEventsWithShowers.query('_fNChargedPiCandidates == 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates > 1 and _fNOtherFancyPairs == 0')
multiPairEventsExt     = extEventsWithShowers.query('_fNChargedPiCandidates == 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates > 1 and _fNOtherFancyPairs == 0')
multiPairEventsData    = dataEventsWithShowers.query('_fNChargedPiCandidates == 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNPairCandidates > 1 and _fNOtherFancyPairs == 0')

#Stack = Plotter.SmartStack(multiPairEvents, '_fMultiPairPi0Mass')
TagLeadingShower(multiPairEventsOverlay)
TagLeadingShower(multiPairEventsDirt)
TagLeadingShower(multiPairEventsExt)
TagLeadingShower(multiPairEventsData)
#print multiPairEvents.loc[:100000, ['passesSubLeadingEnergy', 'passesSubLeadingConvD', '_fMultiPairPi0Mass'] ]


#print multiPairEvents
#raw_input()

#Stack = Plotter.SmartStack(multiPairEvents, '_fMultiPairPi0Mass')
#Plotter.makeMCStackedHistogram(Stack[0], Stack[1], InvariantMassRange, 25, "MultiPairPi0Mass", **{"xlimits" : (0, 200), "legend" : Stack[2]} )

#Plotter.makeMCOnlyStack(multiPairEvents, '_fMultiPairPi0Mass', InvariantMassRange, 25, "MultiPairPi0Mass", {}, {"xlimits" : (0, 200)} )


axesLabels = ["Pass Subleading Conversion Distance Cut", "Number passing", "Number of Events"]
limits = {"xlimits" : NumPairsRange, "ylimits" : (15,), "Titles" : axesLabels }
Plotter.makeDataMCStack(multiPairEventsOverlay, multiPairEventsDirt, multiPairEventsExt, multiPairEventsData, 'NumPassConv', NumPairsRange, 10, "MultiPairPassConversion", {}, limits)
axesLabels[0] = "Pass Subleading Energy Cut"
limits["ylimits"] = (20,)
Plotter.makeDataMCStack(multiPairEventsOverlay, multiPairEventsDirt, multiPairEventsExt, multiPairEventsData, 'NumPassE', NumPairsRange, 10, "MultiPairPassEnergy", {}, limits)
axesLabels[0] = "Pass Invariant Mass Cut"
limits["ylimits"] = (25,)
Plotter.makeDataMCStack(multiPairEventsOverlay, multiPairEventsDirt, multiPairEventsExt, multiPairEventsData, 'NumPassMass', NumPairsRange, 10, "MultiPairPassMass", {}, limits)

#HiMass
highMassEvents = eventsWithShowers.query('_fNChargedPiCandidates == 0 and _fNShowers >= 2 and _fNLeadCandidates >= 1 and _fNSubleadCandidates >= 1 and _fNHiMassPairs == 1 and _fNOtherFancyPairs == 0')
Plotter.makeMCOnlyStack(highMassEvents, '_fHiMassPi0Mass',InvariantMassRange, 25, "HiMassPi0Mass", {}, {"xlimits" : (175, 400)})
#Stack = Plotter.SmartStack(highMassEvents, '_fHiMassPi0Mass')
#Plotter.makeMCStackedHistogram(Stack[0], Stack[1], InvariantMassRange, 25, "HiMassPi0Mass", legend=Stack[2], xlimits=(175, 400))

#print eventsWithShowers.loc[:1000, ['LeadingShowerIdx', 'passesEnergy', 'passesAngle', 'passesConvD']]


sys.exit()
