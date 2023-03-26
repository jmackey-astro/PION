
#args = ("-np", "1")
#OpenComputeEngine("localhost", args)
#SetWindowArea(0,0,800,800)
#DeleteAllPlots()

fdir="./"
fbase="ecc05_d2l3n128"
list_files=[]
list_plots=[]
for i in range(0,3):
  #print(i,str(i).zfill(2))
  files="localhost:"+fdir+fbase+"_level"+str(i).zfill(2)+"_0000.*.silo database"
  OpenDatabase(files, 0)
  AddPlot("Pseudocolor", "Density", 1, 1)
  list_files.append(files)
  list_plots.append(i)

CreateDatabaseCorrelation("Correlation01", list_files, 0)
DrawPlots()
TimeSliderNextState()
SetActivePlots((0, 1, 2))
PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.scaling = PseudocolorAtts.Log  # Linear, Log, Skew
PseudocolorAtts.skewFactor = 1
PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
PseudocolorAtts.minFlag = 1
PseudocolorAtts.min = 1e-27
PseudocolorAtts.maxFlag = 1
PseudocolorAtts.max = 3e-22
PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
PseudocolorAtts.colorTableName = "viridis"
PseudocolorAtts.invertColorTable = 0
PseudocolorAtts.legendFlag = 1
PseudocolorAtts.lightingFlag = 1
SetPlotOptions(PseudocolorAtts)
TimeSliderNextState()
OpenGUI()


