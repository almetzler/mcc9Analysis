import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

'''
This class defines some common plotting tasks and automates some common calculations
Stacked plots, data / MC comparisons etc.
'''

class PlotUtils:
   
   def __init__(self, eventList, mcKey):
      self.nEventTypes = len(eventList)
      self.eventDict   = dict(zip(eventList, range(self.nEventTypes)))    
      self.mcKey = mcKey
      self.good_colors = ["mediumblue", "darkgreen", "indianred", "maroon", "peru", "darkorange", "goldenrod", "yellowgreen", "darkcyan", "dodgerblue",  "indigo", "gainsboro", "dimgray" ]
      self.color_map = dict(zip(eventList, self.good_colors[:self.nEventTypes] ) )
      self.wgtName = ''
      # ("1 Pi0 0 PiP", "1 Pi0 N PiP", "CC 0 Pi0", "CC N Pi0", "CC Other", "NC 0 Pi0", "NC N Pi0", "NC Other", "OOFV", "Non numu")

   def getEventDict(self):
      return self.eventDict
   
   def setWeightName(self, name):
      self.wgtName = name

   def tagDuplicateEvents(inputFrame):
     inputFrame.insert(inputFrame.shape[1], "DuplicatedEvent", inputFrame.duplicated() )

   def loadMCEventInfo(inputFrame):
     tagDuplicateEvents(inputFrame)
     inputFrame.insert(inputFrame.shape[1], "mc_channel", [getChan(x, y) for x, y in zip(inputFrame['mc_nu_interaction_type'], inputFrame['mc_nu_ccnc'])] ) #Classify neutrino events based on CC / NC and event Type
     inputFrame.eval('mc_Ehad = mc_nu_energy - mc_nu_lepton_energy', inplace=True) #Insert the true energy transfer (nu)
     inputFrame.insert(inputFrame.shape[1], "mc_expQ2", [getQ2(x, y, z) for x, y, z in zip(inputFrame['mc_nu_energy'], inputFrame['mc_nu_lepton_energy'], inputFrame['mc_nu_lepton_theta'])] )
     inputFrame.insert(inputFrame.shape[1], "mc_expW", [getW(x, y, z) for x, y, z in zip(inputFrame['mc_Ehad'], inputFrame['mc_expQ2'], inputFrame['mc_channel'] ) ] )
     inputFrame.insert(inputFrame.shape[1], "mc_expXbj", [getXbj(x, y) for x, y in zip(inputFrame['mc_Ehad'], inputFrame['mc_expQ2'] ) ] )
     inputFrame.insert(inputFrame.shape[1], "mc_expY", [getInel(x, y) for x, y in zip(inputFrame['mc_Ehad'], inputFrame['mc_nu_energy'] ) ] )
     #inputFrame.insert(inputFrame.shape[1], "template_wgt", [getChanWeight(x, y) for x, y in zip(inputFrame['mc_nu_interaction_type'], inputFrame['mc_nu_ccnc'])] ) #Classify neutrino events based on CC / NC and event Type
     inputFrame.insert(inputFrame.shape[1], "pot", mcPOT)

   def loadMCTrackInfo(inputFrame):
     inputFrame.insert(inputFrame.shape[1], 'particle', [getParticle(x) for x in inputFrame['mc_pdg'] ])
  
   def loadTrackInfo(inputFrame, isMC=False):
     inputFrame.insert(inputFrame.shape[1], "DuplicatedEvent", inputFrame.duplicated())
     inputFrame.insert(inputFrame.shape[1], "phi", [getPhi(x, y) for x, y in zip(inputFrame['track_diry'], inputFrame['track_dirx'] ) ] )
     inputFrame.eval('track_chi2_ratio = track_chi2_proton / track_chi2_muon', inplace=True)
     inputFrame.insert(inputFrame.shape[1], "isContained", [isContained(x, y, z, a, b, c) for x, y, z, a, b, c in zip(inputFrame['vx'], inputFrame['vy'], inputFrame['vz'], inputFrame['track_endx'], inputFrame['track_endy'], inputFrame['track_endz'])])
     inputFrame.insert(inputFrame.shape[1], 'track_mom_best', [getBestP(x, y, z) for x, y, z in zip(inputFrame['isContained'], inputFrame['track_range_mom_mu'], inputFrame['track_mcs_mom'])])
     if(isMC):
      inputFrame.insert(inputFrame.shape[1], 'particle', [getParticle(x) for x in inputFrame['mc_pdg'] ])

   def makeMCHistogram(self, mc, channel, binRange, nBins, filename, Titles):
     dir_name = "PlotDir/MC_Only"
     colors = {0:'b', 1:'g', 2:'y', 3:'r', 4:'grey', 5:'magenta'}
   
     plt.hist(mc, bins=nBins, stacked=False, range=binRange, color = colors[channel])
     plt.legend([channel])
     try:
       plotTitle, xAxisTitle, yAxisTitle = Titles
     except(ValueError):
       print "Pleast provide three titles for plot name, x and y axis names" 
       plotTitle  = ""
       xAxisTitle = ""
       yAxisTitle = "" 
     plt.title(plotTitle)
     plt.xlabel(xAxisTitle)
     plt.ylabel(yAxisTitle)
     plt.savefig("%s/%s%d.png" % ( dir_name, filename, channel ) )
     plt.close()
   
   def makeMCStackedHistogram(self, mcList, mcWeights, binRange, nBins, filename, **kwargs):
     dir_name = "PlotDir/StackedPlots"
     legend = []
     xLimits = []
     yLimits = []
     if "legend" in kwargs:
      legend = kwargs["legend"]
     if "xlimits" in kwargs:
      xLimits = kwargs["xlimits"]
     if "ylimits" in kwargs:
      yLimits = kwargs["ylimits"]  

     if(len(xLimits) == 1):
       xLimits.insert(0, 0.0)
     if(len(yLimits) == 1):
       yLimits.insert(0, 0.0)

     plt.hist(mcList, bins=nBins, stacked=True, range=binRange, weights = mcWeights )
     plt.legend(legend)
     if "plotTitle" in kwargs:
       plt.title(kwargs[plotTitle])
     if "xAxisTitle" in kwargs:
       plt.xlabel(kwargs[xAxisTitle])
     if "yAxisTitle" in kwargs:
       plt.ylabel(kwargs[yAxisTitle])

     if(len(xLimits) == 2):
       plt.xlim(xLimits)
     if(len(yLimits) == 2):
       plt.ylim(yLimits)

     plt.savefig("%s/%s.png" % ( dir_name, filename) )
     plt.close()
        
   
   def make2DMCHistogram(mc, channel, binRange, nBins, filename, Titles):
     dir_name = "PlotDir"
     colors = {"QE":'b', "RES":'g', "DIS":'y', "2p2h":'r', "NC / Other":'grey', "Ext":'magenta'}
     zMin = 0.01
   
     try:
       xBins, yBins = binRange
   
     except(ValueError):
       print "Please provide a range of bins for each axis"
       return
   
     plt.hist2d(mc[0], mc[1], bins=nBins, range=binRange, cmin=zMin, normed=True)
     plt.legend(legend)
   
     try:
       plotTitle, xAxisTitle, yAxisTitle = Titles
     except(ValueError):
       print "Pleast provide three titles for plot name, x and y axis names" 
       plotTitle  = ""
       xAxisTitle = ""
       yAxisTitle = "" 
     plt.title(plotTitle)
     plt.xlabel(xAxisTitle)
     plt.ylabel(yAxisTitle)
     plt.savefig("%s/%s%s.png" % ( dir_name, filename, channel.replace(" / Other", "")) )
     plt.close()
   
   def makeDataMCHistogram(self, mcList, mcWeights, dataList, binRange, nBins, filename, legend, colors, **kwargs):
     dir_name = "PlotDir/DataMC"
     xLimits  = []
     yLimits  = []
     plotTitles = ["", "", ""]
     rRange = []
     if "Titles" in kwargs:
       plotTitles = kwargs["Titles"]

     if "xlimits" in kwargs:
       for x in kwargs["xlimits"]:
         xLimits.append(x)       
     
     if "ylimits" in kwargs:
       for y in kwargs["ylimits"]:
         yLimits.append(y) 
     
     if "rRange" in kwargs:
       rRange = kwargs["rRange"]

     if(len(xLimits) == 1):
       xLimits.insert(0, 0.0)
     if(len(yLimits) == 1):
       yLimits.insert(0, 0.0)
     
     plt.hist(mcList, bins=nBins, stacked=True, range=binRange, weights=mcWeights, color=colors )
     plt.legend(legend)
     
     

     plt.title(plotTitles[0])
     plt.xlabel(plotTitles[1])
     plt.ylabel(plotTitles[2])
     
     if(len(xLimits) == 2):
       plt.xlim(xLimits)
     if(len(yLimits) == 2):
       plt.ylim(yLimits)
     
     data_hist = self.dataify(dataList, nBins, binRange)
     plt.errorbar(data_hist[0], data_hist[1], yerr=data_hist[2], fmt='o', color='black')
     plt.savefig("%s/%s.png" % ( dir_name, filename) )
     plt.close()
   
     self.makeDataMCRatioHistogram(mcList, mcWeights, dataList, binRange, nBins, filename, plotTitles, xLimits, rRange)
   
   def makeDataMCRatioHistogram(self, mcList, mcWeights, dataList, binRange, nBins, filename, plotTitles, xRange = [], yRange = []):
     dir_name  = "PlotDir/DataMC"
     xLimits  = []
     yLimits  = []
     for x in xRange:
       xLimits.append(x)
     for y in yRange:
       yLimits.append(y)  
     if(len(xLimits) == 1):
       xLimits.insert(0, 0.0)
     if(len(yLimits) == 1):
       yLimits.insert(0, -1.0)
     mcSum = np.full(nBins, 0.0 )
     for mc, weight in zip(mcList, mcWeights):
        mc_hist   = np.histogram(mc, bins=nBins, range=binRange, weights = weight )
        np.add(mc_hist[0], mcSum, out=mcSum)
     data_hist = self.dataify(dataList, nBins, binRange)
     MCScalarSum   = np.sum(mcSum)
     DataScalarSum = np.sum(data_hist[1])
     sumRatio = DataScalarSum / MCScalarSum 
     ratio = np.divide(data_hist[1], mcSum)
     err   = np.multiply(ratio, np.divide(1.0, data_hist[2]))
     np.nan_to_num(ratio, copy=False)
     np.nan_to_num(err, copy=False)
   
     fig, axi = plt.subplots() #create subplots so I can put a textbox in
   
     plt.errorbar(data_hist[0], ratio, yerr=err, fmt='o', color='black') #This ignores MC stats.
     plt.title(plotTitles[0])
     plt.xlabel(plotTitles[1])
     plt.ylabel("Data / MC")
     
     if(len(xLimits) == 2):
       plt.xlim(xLimits)
     if(len(yLimits) == 2):
       plt.ylim(yLimits)
   
     text = r'$\int \frac{data}{MC} = %.3f$' % sumRatio
     props = dict(boxstyle='round', facecolor='lightsteelblue', alpha=0.5)
   
     axi.text(0.75, 1.1, text, transform=axi.transAxes, fontsize=14, verticalalignment='top', bbox=props)
     plt.savefig("%s/%sRatio.png" % ( dir_name, filename) )
     plt.close()
     
   def KEtoMomentum(T, restMass):
       TotalE = T + restMass
       return math.sqrt(TotalE**2 - restMass**2)
       '''
       0.2 + 0.1 = 0.3 
       sqrt(0.3^2 - 0.1^2) = sqrt(0.09 - 0.01)
       sqrt(0.08) = 0.28
       '''
   def dataify(self, array, bins, limits):
      counts, bin_edges = np.histogram(array, bins, range=limits)
      bin_centers       = (bin_edges[:-1] + bin_edges[1:]) / 2
      errs = np.sqrt(counts)
      return [bin_centers, counts, errs]
      
   def getChan(interaction, isNC):
       if (isNC):
         return "NC / Other"
       elif(interaction == 0):
         return "QE"
       elif(interaction == 1):
         return "RES"
       elif(interaction == 2):
         return "DIS"
       elif(interaction == 10):
         return "2p2h"
   
   def getParticle(pdg):
       if (pdg == 13):
         return "muon"
       elif(pdg == 2212):
         return "proton"
       elif(pdg == 211):
         return "pion"
       elif(pdg == 11):
         return "electron"
       elif(pdg == -13):
         return "muon+"  
       else:
         return "other"      
   
   def getChanWeight(interaction, isNC):
       templateFactors = [6.47365e-01, 1.20327e-07, 6.02801e-08, 2.71038e-08, 3.42514e-01, 1.00000e-02]
       if (isNC):
         return (1.0 + templateFactors[4])
       elif(interaction == 0):
         return (1.0 + templateFactors[0])
       elif(interaction == 1):
         return (1.0 + templateFactors[1])
       elif(interaction == 2):
         return (1.0 + templateFactors[2])
       elif(interaction == 10):
         return (1.0 + templateFactors[3])
   
   
   
   def dotProduct(df1, df2, dim1, dim2, dimList):
      returnSeries = pd.Series(0.0, index=df1.index.copy())
     # print df1['%s' % dim1]*df1['%s' % "track_dirx"]
     # print df2['%s' % dim2]*df2['%s' % "track_dirx"]
     # print df1['%s' % dim1]*df1['%s' % "track_dirx"] + df2['%s' % dim2]*df2['%s' % "track_dirx"]
      denomSeries = df1['%s' % dim1]*df2['%s' % dim2]
      #print denomSeries
   
      for dim in dimList:    
          returnSeries += (df1['%s' % dim1]*df1['%s' % dim]*df2['%s' % dim2]*df2['%s' % dim])/denomSeries
          #print returnSeries
      #print returnSeries
      return returnSeries.to_numpy()
   
   def getPhi(pY, pX):
       return (np.arctan2(pY, pX) / np.pi)
   
   def getTheta(pTotal, pZ):
       return (pZ / pTotal)
   
   def getQ2(nuE, nuMu, thetaMu):
       return 4*nuE*nuMu*math.pow(math.sin(thetaMu/2), 2)
   
   #Containment taken from mcc8 Nproton analysis, note 1065
   def isContained(xStart, yStart, zStart, xEnd, yEnd, zEnd):
       return (checkContained(xStart, yStart, zStart) and checkContained(xEnd, yEnd, zEnd))
   
   def isFiducial(x, y, z):
     if(x > 245.0 or x < 10.0):
         return False
     elif(y > 100.0 or y < -100.0):  
         return False
     elif(z > 970.0 or z <10.0):
         return False
     else:  
         return True
   
   def checkContained(x, y, z):
       if(x > 245.0 or x < 10.0):
         return False
       elif(y > 100.0 or y < -100.0):  
         return False
       elif(z > 1030.0 or z <10.0):
         return False
       else:  
         return True
   
   def getBestP(isContained, pOne, pTwo):
     return (pOne if isContained else pTwo)
   
   def getW2(Ehad, Q2, interaction):
       targetMass = 0.989
       #if it's a 2p 2 h, you have 2ps as your target! Wakka Wakka Wakka!
       if(interaction == "2p2h"):
         targetMass = 2*targetMass
       return 2*targetMass*Ehad + math.pow(targetMass, 2) - Q2
   
   def getW(Ehad, Q2, interaction):
       W2 = getW2(Ehad, Q2, interaction)
       if(W2 >= 0.0):
         return math.sqrt(W2)
       else:
         return -1.0
   
   def getXbj(Ehad, Q2):
       targetMass = 0.989
       if(Ehad > 0.0):
         return Q2 / (2*targetMass*Ehad)
       else:
         return -1.0
   
   def getInel(Ehad, Enu):
       if(Enu > 0.0):
         return (Ehad / Enu)
       else:
         return -1.0
 
   def SmartStack(self, dataFrame, var, **kwargs):
     plotList = {}
     weightList = {}
     sortPlots = []
     sortWeights = []
     legend = []
     colors = []
     preSort = {}
     signalIDX = 0
     signalName = ""
     if "signalName" in kwargs:
       signalName = kwargs[signalName]
     if "signalOnTop" in kwargs:
       signalOnTop = kwargs[signalOnTop]   

     q_attribute = self.mcKey
     for key in self.eventDict:
       eventNo = self.eventDict[key] 
       if key == signalName:
        signalIDX = eventNo
       call = '{} == {}'.format(q_attribute,eventNo)
       if not dataFrame.query(call)[var].empty :
         array  = dataFrame.query(call)[var].apply(pd.Series).stack().reset_index(drop=True).to_numpy()
       else:
         array = np.array([])

       if not self.wgtName:
        weights = np.full(array.size, 1.0)
       else:
        if not dataFrame.query(call)[var].empty :
         weights = dataFrame.query(call)[self.wgtName].apply(pd.Series).stack().reset_index(drop=True).to_numpy()
        else:
         weights = np.array([])
       preSort[key] = array.size
       plotList[key] = array #Maps 0 in dict to 0 in array
       weightList[key] = weights
     postSort = sorted(preSort.items(), key=lambda x: x[1], reverse=True)

     for chan, num in postSort:
        sortPlots.append(plotList[chan])
        sortWeights.append(weightList[chan])
        legend.append(chan)
        colors.append(self.color_map[chan])      
     #print sortPlots[3]
     if signalName:
      if signalOnTop:
         sortPlots.insert(0, sortPlots.pop(signalIDX))
         sortWeights.insert(0, sortWeights.pop(signalIDX))
         legend.insert(0, legend.pop(signalIDX))  

     return (sortPlots, sortWeights, legend, colors)
   
   def makeMCOnlyStack(self, dataframe, var, binRange, nBins, filename, stackArgs, plotArgs):
     Stack = self.SmartStack(dataframe, var, **stackArgs)
     plotArgs["legend"] = Stack[2]
     self.makeMCStackedHistogram(Stack[0], Stack[1], binRange, nBins, filename, **plotArgs)
     return Stack
   
   def makeDataMCStack(self, mcFrame, dirtFrame, extFrame, dataFrame, var, binRange, nBins, filename, stackArgs, plotArgs, **kwargs):
     Stack = self.SmartStack(mcFrame, var, **stackArgs)
     if not dirtFrame[var].empty :
         array  = dirtFrame[var].apply(pd.Series).stack().reset_index(drop=True).to_numpy()
         wgt    = dirtFrame[self.wgtName].apply(pd.Series).stack().reset_index(drop=True).to_numpy()
     else:
         array = np.array([])
         wgt = np.array([])
     Stack[0].append(array)
     Stack[1].append(wgt)
     Stack[2].append("Dirt")
     Stack[3].append("darkgoldenrod")

     if not extFrame[var].empty :
         array  = extFrame[var].apply(pd.Series).stack().reset_index(drop=True).to_numpy()
         wgt    = extFrame[self.wgtName].apply(pd.Series).stack().reset_index(drop=True).to_numpy()
     else:
         array = np.array([])
         wgt   = np.array([])
     
     if "extOnTop" in kwargs:
       Stack[0].insert(0, array)
       Stack[1].insert(0, wgt)
       Stack[2].insert(0, "Ext.")
       Stack[3].insert(0,"magenta")
     
     else:  
       Stack[0].append(array)
       Stack[1].append(wgt)
       Stack[2].append("Ext.")
       Stack[3].append("magenta")

     #plotArgs["legend"] = Stack[2]
     dataList = dataFrame[var].to_numpy()

     self.makeDataMCHistogram(Stack[0], Stack[1], dataList, binRange, nBins, filename, Stack[2], Stack[3], **plotArgs)
     return Stack

   def Stack(dataframe, dirtDF, extDF, variable, stackSize=0):
     retlist=[]
   
     q_attribute = self.mcKey
     
     dirt = dirtDF[variable].to_numpy()
     ext = extDF[variable].to_numpy()
     if(stackSize == 0):
        stackSize = self.nEventTypes
     value_list = range(stackSize)
   
     for value in value_list:
         call = '{} == "{}"'.format(q_attribute,value) + addons
         item = dataframe.query(call)[variable].to_numpy()
         retlist.append(item)
     retlist.append(dirt)
     retlist.append(ext)
     return retlist         
