
#args = ("-np", "1")
#OpenComputeEngine("localhost", args)
#SetWindowArea(0,0,800,800)
#DeleteAllPlots()

OpenDatabase("localhost:V444Cyg_d2l3n320_level00_0000.*.silo database", 0)
AddPlot("Pseudocolor", "Density", 1, 1)
OpenDatabase("localhost:V444Cyg_d2l3n320_level01_0000.*.silo database", 0)
AddPlot("Pseudocolor", "Density", 1, 1)
OpenDatabase("localhost:V444Cyg_d2l3n320_level02_0000.*.silo database", 0)
AddPlot("Pseudocolor", "Density", 1, 1)
CreateDatabaseCorrelation("Correlation01",("localhost:V444Cyg_d2l3n320_level00_0000.*.silo database", "localhost:V444Cyg_d2l3n320_level01_0000.*.silo database", "localhost:V444Cyg_d2l3n320_level02_0000.*.silo database"), 0)
DrawPlots()
TimeSliderNextState()
SetActivePlots((0, 1, 2))
PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.scaling = PseudocolorAtts.Log  # Linear, Log, Skew
PseudocolorAtts.skewFactor = 1
PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
PseudocolorAtts.minFlag = 1
PseudocolorAtts.min = 2e-16
PseudocolorAtts.maxFlag = 1
PseudocolorAtts.max = 2e-11
PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
PseudocolorAtts.colorTableName = "magma"
PseudocolorAtts.invertColorTable = 0
PseudocolorAtts.legendFlag = 1
PseudocolorAtts.lightingFlag = 1
SetPlotOptions(PseudocolorAtts)
TimeSliderNextState()
#OpenGUI()


