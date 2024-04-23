"""
/***************************************************************************
Date                 : August/2023
email                : eric.w.neilson@nrcan-rncan.gc.ca
reference:



## run the following in teh console
py3_env
python3 -m pip install --upgrade pip
python3 -m pip install statsmodels
python3 -m pip install numpy
python3 -m pip install scipy
python3 -m pip install keyboard

py3_env
python3 "C:/Users/erneilso/OneDrive - NRCan RNCan/Collaborations/ROF/SpatialPowerAnalysis/reSample.py"
 ***************************************************************************/
/***************************************************************************
 some useful updates
 
 conda update -n base -c defaults conda
 conda install gdal
 
   - 
 ***************************************************************************/
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

import statsmodels as sm
import math
from osgeo import gdal, ogr, gdalconst
import numpy as np
import os
import sys
import keyboard

import json

Max_Distance= -1.0
options = []
if Max_Distance > 0: options.append('MAXDIST=' + str(Max_Distance))


# ## https://gis.stackexchange.com/questions/326968/ogr2ogr-error-1-proj-pj-obj-create-cannot-find-proj-db
# os.environ['PROJ_LIB'] = "C:/Users/erneilso/AppData/Local/anaconda3/Library/share/proj"
# os.environ['GDAL_DATA'] = "C:/Users/erneilso/AppData/Local/anaconda3/Lib/site-packages/osgeo"
# #GDAL_DATA="C:/Users/erneilso/AppData/Local/anaconda3/Lib/site-packages/osgeo/"

rdriver = gdal.GetDriverByName('GTiff')
vdriver = ogr.GetDriverByName('ESRI Shapefile')

# =============================================================================
# functions
# =============================================================================
    
def shapetoRasterSA(inPath,inputName,rType):

    global outPath; global nodata; global cell_size; global extent; global srs; global n_rows; global n_cols
    
    print("converting shapefile to raster")
    shpIn = vdriver.Open(inPath,0)
    layer = shpIn.GetLayer()   
    pathOut = outPath + "/" + inputName + "_" + str(cell_size) + "_m" + '.tif'
    rasterizedDS = rdriver.Create(pathOut, n_rows, n_cols, 1, gdal.GetDataTypeByName(rType))
    rasterizedDS.SetGeoTransform([extent[0], cell_size, 0, extent[3], 0, -1 * cell_size])    
    rasterizedDS.SetProjection(srs.ExportToWkt())
    rasterizedDS.GetRasterBand(1).SetNoDataValue(nodata)  
    gdal.RasterizeLayer(rasterizedDS, [1], layer)
    print("Shape to raster complete")
    
    band1 = rasterizedDS.GetRasterBand(1)
    list1 = band1.ReadAsArray()
    return list1

    layer = None
    shpIn = None

def shapeProxRaster(inPath,inputName):

    global outPath; global nodata; global cell_size; global extent; global srs; global n_rows; global n_cols
    
    print("converting shapefile to raster")
    shpIn = vdriver.Open(inPath,0)
    layer = shpIn.GetLayer()     
    pathOut = outPath + "/intermed_" + inputName + "_" + str(cell_size) + "_m" + '.tif'
    rasterizedDS = rdriver.Create(pathOut, n_rows, n_cols, 1, gdal.GDT_Byte)
    rasterizedDS.SetGeoTransform([extent[0], cell_size, 0, extent[3], 0, -1 * cell_size])    
    rasterizedDS.SetProjection(srs.ExportToWkt())
    rasterizedDS.GetRasterBand(1).SetNoDataValue(nodata)  
    gdal.RasterizeLayer(rasterizedDS, [1], layer)
    rasterizedDS = None
    layer = None
    shpIn = None
    
    print("converting raster to proximity") # proximity raster   https://gis.stackexchange.com/questions/220753/how-do-i-create-blank-geotiff-with-same-spatial-properties-as-existing-geotiff
    rasIn = gdal.Open(pathOut, gdalconst.GA_ReadOnly)
    prox_path = outPath + "/" + inputName + "_dist_" + str(cell_size) + "_m" + '.tif'
    proximityDs = rdriver.Create(prox_path, n_rows, n_cols, 1, gdal.GetDataTypeByName("Int32"))   
    proximityDs.SetGeoTransform([extent[0], cell_size, 0, extent[3], 0, -1 * cell_size])    
    proximityDs.SetProjection(srs.ExportToWkt())
    proximityDs.GetRasterBand(1).SetNoDataValue(nodata)
    proximityDs.GetRasterBand(1).Fill(nodata)
    gdal.ComputeProximity(rasIn.GetRasterBand(1), proximityDs.GetRasterBand(1), options, callback = None)
    print("finsihed proximity raster") 

    band1 = proximityDs.GetRasterBand(1)
    list1 = band1.ReadAsArray()
    return list1    
    
    pathOut = None
    rasIn = None

def rasterSA(inPath,inputName):

    global outPath; global nodata; global cell_size; global extent; global srs; global n_rows; global n_cols

    print("resampling raster to study area")
    rasIn = gdal.Open(inPath, gdalconst.GA_ReadOnly)
    inputProj = rasIn.GetProjection()
    #inputType = rasIn.GetMetadata()[0]
    inputType = "Int32"
    pathOut = outPath + "/" + inputName + "_" + str(cell_size) + "_m" + '.tif'
    rasterizedDS = rdriver.Create(pathOut, n_rows, n_cols, 1, gdal.GetDataTypeByName(inputType))
    rasterizedDS.SetGeoTransform([extent[0], cell_size, 0, extent[3], 0, -1 * cell_size])    
    rasterizedDS.SetProjection(srs.ExportToWkt())
    rasterizedDS.GetRasterBand(1).SetNoDataValue(nodata)
    gdal.ReprojectImage(rasIn,rasterizedDS,inputProj,srs.ExportToWkt(),gdalconst.GRA_NearestNeighbour)
    print("finsihed reprojecting raster") 
    
    band1 = rasterizedDS.GetRasterBand(1)
    list1 = band1.ReadAsArray()
    return list1
    
    layer = None
    rasIn = None

def normScaler (inArray):
    outArray = (inArray - np.min(inArray)) / (np.max(inArray) - np.min(inArray))
    return outArray
    
def proximityCalc (inArray):
    outArray = np.max(inArray) - inArray
    return outArray

def inputSpatial():

    global responseVar
    print('')
    print("You can enter a series of spatial data files, including a study area shapefile, to be standardized, transformed and stacked to conduct a power analysis for data design. The files must be entered using a json file with the following format.\n")    

    print("All analysis will be done using a raster stack. The input json file has a field for indicating the file path to a study area shapefile. The projection and extent of the rasters will match that of the input study area shapefile and the cell size is input by the user in the json file. The json also requires the path to an (existing) output folder to which all subseqeuent spatial analytical files will be written.")
    print('')
    print("The predictors list has entries for layers to be used in modeling " + responseVar + ", each of which requires a user defined variable name, the path to the input file, the type (vector or raster) and type of analysis transformation, includeing:\n (Scale) values scaled between 0 and 1 \n (Proximity) binary raster or shapefile converted to 'distance-to' \n Etc. \n")
    
    openJson = input(" (1) Open example json file")
    if openJson == "1":
        os.startfile(r'C:\Users\erneilso\OneDrive - NRCan RNCan\Collaborations\ROF\SpatialPowerAnalysis\inJson\inputFile1.json')
    else:
        print("Use only numbers available.")
        readQuestion(strIndex)
        
    print('')
    # print('')
    # inputFile = input("Enter path and file of input json: ")
    # if inputFile == "1":
        # inputFile = r'C:\Users\erneilso\OneDrive - NRCan RNCan\Collaborations\ROF\SpatialPowerAnalysis\inJson\inputFile1.json'
    inputFile = r'C:\Users\erneilso\OneDrive - NRCan RNCan\Collaborations\ROF\SpatialPowerAnalysis\inJson\inputFile1.json'
    
    # processing json
    f = open(inputFile)
    inData = json.load(f) ;f.close()
    global outPath
    outPath = inData["outPath"]

    ## raster info ##
    global nodata; global cell_size; global extent; global srs; global n_rows; global n_cols
    nodata = -9999
    cell_size = inData["cellSize"]
    SAshp = vdriver.Open(inData["SA"],0)  ## this made here: C:/Users/erneilso/OneDrive - NRCan RNCan/Collaborations/ROF/ROF_Map.qgz
    rlayer = SAshp.GetLayer()
    extent = rlayer.GetExtent()
    srs = rlayer.GetSpatialRef()
    n_rows = int(math.ceil(abs(extent[1] - extent[0]) / cell_size))
    n_cols = int(math.ceil(abs(extent[3] - extent[2]) / cell_size)) 
    
    print("The study area has " + str(extent)  + ". \nIt has " + str(n_rows) + " rows and "  + str(n_cols) + " columns." )

    ## raster info ##
    predictors = inData['predictors']
    
    global predictorValues
    predictorValues = {}
    
    for pred in predictors:
        
        if pred["Type"] == "Vector":
            if pred["Transform"] == "Proximity":
                tempArray = shapeProxRaster(pred["Path"],pred["Name"])
                predictorValues[pred["Name"]] = proximityCalc(tempArray) 
            
        elif pred["Type"] == "Raster":
            tempArray = rasterSA(pred["Path"],pred["Name"])    
            if pred["Transform"] == "Scale":
                predictorValues[pred["Name"]] = normScaler(tempArray)                 
                
        print('')

    print("Predictor names:", list(predictorValues.keys()))
    
def simulateReponse():
    
    global predictorValues
    global responseVar
    global outPath
    
    varArray = []
    for i in predictorValues:
        varArray.append(predictorValues[i])
    varArray = np.array(varArray) 
    
    print(" Calculating probability of use across study area")
    predictorValues["response"] = np.prod(varArray, axis=0)
    #print(predictorValues["response"])
    
    # for i in predictorValues:
        # print(i)
        # for j in range(10):            
            # print(predictorValues[i][j][range(5)])

    ## make response raster
    responsePath = outPath + "/" + responseVar + "_" + str(cell_size) + "_m" + '.tif'
    rasterizedDS = rdriver.Create(responsePath, n_rows, n_cols, 1, gdal.GetDataTypeByName("Float32"))
    rasterizedDS.SetGeoTransform([extent[0], cell_size, 0, extent[3], 0, -1 * cell_size])    
    rasterizedDS.SetProjection(srs.ExportToWkt());rasterizedDS.GetRasterBand(1).SetNoDataValue(nodata);rasterizedDS.GetRasterBand(1).Fill(nodata)
    rasterizedDS.GetRasterBand(1).WriteArray(predictorValues["response"])

def testPower():
    
    global nodata; global cell_size; global extent; global srs; global n_rows; global n_cols
    global predictorValues
    global prevPos
    global responseVar
    
    camN = input("Enter the number of cameras (number).")
    try:
        int(camN)
        print(camN + " cameras")
    except:
        print("Enter a number.")
        testPower()       
    print('')
    
    camConfig = input("Enter the configuration of cameras \n (1) Systematic\n (2) Random\n (3) Stratigied Random \n")    
    if camConfig.upper() == "B":
        readQuestion(prevPos)    
    if camConfig != "1":
        print("Not available.")
        testPower()
    print('')
    
    # camN = 20
    # camConfig = "Systematic"
    
    rArray = predictorValues["response"]
    rArray = normScaler(rArray) 
    rSize = rArray.size
    # print("shape is  " + str(rArray.shape))
    rSysN = int(round((rSize / int(camN))))
    # print("every " + str(rSysN))
    
    print("getting probabilities")
    detProb = []
    list2 = np.ndarray.flatten(rArray)
    detProb = list2[1::rSysN]
    print("The probabilities of use at cameras are \n" + str(detProb))    
      
    print("exporting map of cameras")
    simCamPath = outPath + "/" + camN + "_cam_" + str(cell_size) + "_m" + '.tif'
    rasterizedDS = rdriver.Create(simCamPath, n_rows, n_cols, 1, gdal.GetDataTypeByName("Int32"))
    rasterizedDS.SetGeoTransform([extent[0], cell_size, 0, extent[3], 0, -1 * cell_size])    
    rasterizedDS.SetProjection(srs.ExportToWkt());rasterizedDS.GetRasterBand(1).SetNoDataValue(nodata);rasterizedDS.GetRasterBand(1).Fill(0)
    band1 = rasterizedDS.GetRasterBand(1)
    list1 = band1.ReadAsArray()
    
    list2 = np.ndarray.flatten(list1)
    list2[1::rSysN] = 1
    list2 = np.reshape(list2,(list1.shape))
    rasterizedDS.GetRasterBand(1).WriteArray(list2)
  
def readQuestion(strIndex):
    
    print('')
    global currPos
    global prevPos
    global answer    
    
    currPos = iData[strIndex]    
    
    answer = input( currPos["Text"] )        
        
    if answer.upper() == "B":
        readQuestion(prevPos)
    
    if currPos["Function"] != "":        
        #print("Calling " + currPos["Function"] ) 
        method = globals() [ currPos["Function"] ]
        method()        
    
    if answer in currPos["posAnswers"]:
        prevPos = strIndex
        nextPos = currPos["posAnswers"].index(answer)
        readQuestion( currPos["ResponseIndices"][nextPos] )

    else:
        print("Use only numbers available.")
        readQuestion(strIndex)
        
def addResponseVar():

    global answer
    global responseVar
    
    if answer == "1":
        responseVar = "Use"
    elif answer == "2":
        responseVar = "Range"
    else:
        responseVar = "Occupancy"
        print("modeling occupancy")

def endTool():
    global answer    
    if answer == "1":
        sys.exit()  


###################
#  start

currPos = "0"
prevPos = "0"
answer = "0"    
    
inputFile = r'C:\Users\erneilso\OneDrive - NRCan RNCan\Collaborations\ROF\SpatialPowerAnalysis\inJson\CYOA.json'
f = open(inputFile)
iData = json.load(f)["DT_Branches"]; f.close()

readQuestion("0")  



# global responseVar
# responseVar = "Occupancy"
# inputSpatial()
# simulateReponse()



# outPath = "C:/Users/erneilso/OneDrive - NRCan RNCan/Collaborations/ROF/SpatialPowerAnalysis/Output" 
# cell_size=5000
# SAshp=vdriver.Open('C:/Users/erneilso/OneDrive - NRCan RNCan/Collaborations/ROF/Data/Shapes/ROF_SA.shp',0)  ## this made here: C:/Users/erneilso/OneDrive - NRCan RNCan/Collaborations/ROF/ROF_Map.qgz

# nodata = -9999
# rlayer = SAshp.GetLayer()
# extent = rlayer.GetExtent()
# srs = rlayer.GetSpatialRef()
# ##print(srs)
# n_rows = int(math.ceil(abs(extent[1] - extent[0]) / cell_size))
# n_cols = int(math.ceil(abs(extent[3] - extent[2]) / cell_size))
# ras1 = gdal.Open(r"C:/Users/erneilso/OneDrive - NRCan RNCan/Collaborations/ROF/SpatialPowerAnalysis/Output/Occupancy_5000_m.tif", gdalconst.GA_ReadOnly)
# band1 = ras1.GetRasterBand(1)
# list1 = band1.ReadAsArray()


# global predictorValues
# predictorValues = {}
# predictorValues["response"] = list1

# global prevPos
# prevPos = "0"

# testPower()












########################
## junk

# print(rasDict)
# ras1 = rasDict[0]
# band1 = ras1.GetRasterBand(1)
# list1 = band1.ReadAsArray()
# waterProx = np.max(list1) - list1
# waterProx = normScaler(waterProx)
# print(waterProx)

# print(dfsdf)

# ras1 = gdal.Open(path_prefix + '/Data/Rasters/waterbody1m_' + str(cell_size) + '.tif', gdalconst.GA_ReadOnly)
# band1 = ras1.GetRasterBand(1)
# water = band1.ReadAsArray()
# # print(water)

# ras1 = gdal.Open(path_prefix + '/Data/Rasters/cand3d30_DEM_ROF_' + str(cell_size) + '.tif', gdalconst.GA_ReadOnly)
# band1 = ras1.GetRasterBand(1)
# list1 = band1.ReadAsArray()
# ele = normScaler(list1)
# # print(ele)

# ############################
# ### occupancy surface
# ############################

# ## maunal inputs
# answer = input("Does the effect of proximity to water on occupancy decay with distance (Y or N? \n")
# if answer.upper() == "Y":
#     logwater = True
# else:
#     logwater = False
# answer = input("What is the strength of proximity to water on occupancy\n")
# waterCoef = float(answer)
# answer = input("What is the strength of elevation on occupancy?\n")
# eleCoef = float(answer) 

# if logwater:
#     waterProx = np.log(waterProx + 1)
    
    
# print("calculate occupancy")
# occArray = eleCoef*ele + waterCoef*waterProx  

# ## occupany surface
# outOccPath = path_prefix + '/Data/Rasters/Occupancy.tif'
# rasterizedDS = rdriver.Create(outOccPath, n_rows, n_cols, 1, gdal.GetDataTypeByName("Float32"))
# rasterizedDS.SetGeoTransform([extent[0], cell_size, 0, extent[3], 0, -1 * cell_size])    
# rasterizedDS.SetProjection(srs.ExportToWkt());rasterizedDS.GetRasterBand(1).SetNoDataValue(nodata);rasterizedDS.GetRasterBand(1).Fill(nodata)
# rasterizedDS.GetRasterBand(1).WriteArray(occArray)
# print('occupancy output')


# src = rasterio.open(path_prefix + '/Data/Rasters/Occupancy.tif')
# pyplot.imshow(src.read(1), cmap='pink')
# pyplot.show()
    
    



            
            
        




# answer = input("Manually enter inputs (Y or N)?")
# print('')
# if answer.upper() == "Y":
    
    # ## output directory
    # outPath = input("Enter path of output directory (no hanging '/'): ")
    # print(outPath)
    # print('')
    
    # cell_size = int(input("Enter grain size (m)?"))
    # print('')  
    # saPath = input("Enter path of Study Area shapefile. All proceeding spatial analysis will have the extent and CRS of SA.")
    # SAshp = vdriver.Open(saPath,0)  ## this made here: C:/Users/erneilso/OneDrive - NRCan RNCan/Collaborations/ROF/ROF_Map.qgz
    # nodata = -9999
    # rlayer = SAshp.GetLayer()
    # extent = rlayer.GetExtent()
    # srs = rlayer.GetSpatialRef()
    # n_rows = int(math.ceil(abs(extent[1] - extent[0]) / cell_size))
    # n_cols = int(math.ceil(abs(extent[3] - extent[2]) / cell_size))
    # print('')
    
    # rasDict = {}
    # vecDict = {}
    # print("Enter other input layers.")
    # inp = True
    # while inp:      
        
        # inputName = input("Enter a variable name for this input.")
        # inputType = input("Enter the file type of the input; shapefile (V) or GTiff (R).")
        
        # if inputType.upper() == "V":
            # var = vdriver.Open(saPath,0)
            
        # elif inputType.upper() == "R":
            # inputPath = input("Enter the full path to the raster.")            
            # rasDict[inputName] = rasterSA(inputPath,inputName)        

        # addMore = input("Add more? yes (Y) or no (N).") 
        # if addMore.upper() == "N":
            # inp = False

        # print('')
        # print('')



# else:
    # # outPath = "C:/Users/erneilso/Documents/LocalProjects/OutFolder"
    # outPath = "C:/Users/erneilso/OneDrive - NRCan RNCan/Collaborations/ROF/SpatialPowerAnalysis/Output" 
    # cell_size=1000
    # SAshp=vdriver.Open('C:/Users/erneilso/OneDrive - NRCan RNCan/Collaborations/ROF/Data/Shapes/ROF_SA.shp',0)  ## this made here: C:/Users/erneilso/OneDrive - NRCan RNCan/Collaborations/ROF/ROF_Map.qgz



    
  
  
  
# outPath = "C:/Users/erneilso/OneDrive - NRCan RNCan/Collaborations/ROF/SpatialPowerAnalysis/Output" 
# cell_size=1000
# SAshp=vdriver.Open('C:/Users/erneilso/OneDrive - NRCan RNCan/Collaborations/ROF/Data/Shapes/ROF_SA.shp',0)  ## this made here: C:/Users/erneilso/OneDrive - NRCan RNCan/Collaborations/ROF/ROF_Map.qgz

# nodata = -9999
# rlayer = SAshp.GetLayer()
# extent = rlayer.GetExtent()
# srs = rlayer.GetSpatialRef()
# ##print(srs)
# n_rows = int(math.ceil(abs(extent[1] - extent[0]) / cell_size))
# n_cols = int(math.ceil(abs(extent[3] - extent[2]) / cell_size))

  
# inPath = "D:/LocalProjects/Data/Water/Waterbody_1M.shp"
# inputName = "WaterProx"  
  
  
# print("converting shapefile to raster")
# shpIn = vdriver.Open(inPath,0)
# layer = shpIn.GetLayer()     
# pathOut = outPath + "/" + inputName + "_" + str(cell_size) + "_m" + '.tif'
# rasterizedDS = rdriver.Create(pathOut, n_rows, n_cols, 1, gdal.GDT_Byte)
# rasterizedDS.SetGeoTransform([extent[0], cell_size, 0, extent[3], 0, -1 * cell_size])    
# rasterizedDS.SetProjection(srs.ExportToWkt())
# rasterizedDS.GetRasterBand(1).SetNoDataValue(nodata)  
# gdal.RasterizeLayer(rasterizedDS, [1], layer)
# shpIn = None
# rasterizedDS = None

# print("converting raster to proximity")
# input = gdal.Open( pathOut, gdalconst.GA_ReadOnly)
# prox_path = outPath + "/" + inputName + "_prox_" + str(cell_size) + "_m" + '.tif'
# proximityDs = rdriver.Create(prox_path, n_rows, n_cols, 1, gdal.GetDataTypeByName("Int32"))   
# proximityDs.SetGeoTransform([extent[0], cell_size, 0, extent[3], 0, -1 * cell_size])    
# proximityDs.SetProjection(srs.ExportToWkt())
# proximityDs.GetRasterBand(1).SetNoDataValue(nodata)
# proximityDs.GetRasterBand(1).Fill(nodata)
# gdal.ComputeProximity(input.GetRasterBand(1), proximityDs.GetRasterBand(1), options, callback = None)
# print("finsihed proximity raster") 

# rasIn = gdal.Open(pathOut,0)
# prox_path = outPath + "/" + inputName + "_prox_" + str(cell_size) + "_m" + '.tif'
# proximityDs = rdriver.Create(prox_path, n_rows, n_cols, 1, gdal.GDT_Float32) 
# proximityDs.SetGeoTransform([extent[0], cell_size, 0, extent[3], 0, -1 * cell_size])    
# proximityDs.SetProjection(srs.ExportToWkt())
# proximityDs.GetRasterBand(1).SetNoDataValue(nodata)
# proximityDs.GetRasterBand(1).Fill(nodata)
# gdal.ComputeProximity(rasIn.GetRasterBand(1), proximityDs.GetRasterBand(1))
# print("finsihed proximity raster") 
    
    

# rasPath = "D:/LocalProjects/Data/Elevation/can3d30_tiff.tif"
# rasIn = gdal.Open(rasPath, gdalconst.GA_ReadOnly)
# inputProj = rasIn.GetProjection()
# pathOut = outPath + "/DEM10_" + str(cell_size) + "_m" + '.tif'


# rasterizedDS = rdriver.Create(pathOut, n_rows, n_cols, 1, gdal.GetDataTypeByName("Int32"))
# rasterizedDS.SetGeoTransform([extent[0], cell_size, 0, extent[3], 0, -1 * cell_size])    
# rasterizedDS.SetProjection(srs.ExportToWkt())
# rasterizedDS.GetRasterBand(1).SetNoDataValue(nodata)
# gdal.ReprojectImage(rasIn,rasterizedDS,inputProj,srs.ExportToWkt(),gdalconst.GRA_NearestNeighbour)
# print("finsihed reprojecting raster") 

# rasIn = None
# pathOut = None
# outPath = None
# rasterizedDS = None
# clipOut = None

