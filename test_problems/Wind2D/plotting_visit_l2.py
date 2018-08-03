# python script to open silo files for each level of the refined
# grid, and create a database correlation for them so that all files
# step through in time together.
# It starts the GUI and quits, leaving the GUI in control of
# plotting.  It just sets up the environment so that the simulation
# is easy to visualise.

#OPDIR="/mnt/data/jm/nested_pion/Wind2D"
OPDIR="/vol/aibn129/data1/jmackey/scratch/Wind2D"

DB_level00="localhost:"+OPDIR+"/OA2_n0128_l2_level00.*.silo database"
DB_level01="localhost:"+OPDIR+"/OA2_n0128_l2_level01.*.silo database"

SetDatabaseCorrelationOptions(0,0)
OpenComputeEngine("localhost", ("-np", "2"))


OpenDatabase(DB_level00, 0)
OpenDatabase(DB_level01, 0)

ActivateDatabase(DB_level00)
AddPlot("Pseudocolor", "Density", 1, 0)
ActivateDatabase(DB_level01)
AddPlot("Pseudocolor", "Density", 1, 0)

CreateDatabaseCorrelation("Correlation01",(DB_level00, DB_level01), 0)

SetActivePlots((0, 1))

PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.scaling = PseudocolorAtts.Log  # Linear, Log, Skew
PseudocolorAtts.skewFactor = 1
PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
PseudocolorAtts.minFlag = 1
PseudocolorAtts.min = 1e-26
PseudocolorAtts.maxFlag = 1
PseudocolorAtts.max = 1e-22
PseudocolorAtts.centering = PseudocolorAtts.Zonal
PseudocolorAtts.colorTableName = "magma"
PseudocolorAtts.invertColorTable = 0
PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque
PseudocolorAtts.opacityVariable = ""
PseudocolorAtts.opacity = 1
PseudocolorAtts.opacityVarMin = 0
PseudocolorAtts.opacityVarMax = 1
PseudocolorAtts.opacityVarMinFlag = 0
PseudocolorAtts.opacityVarMaxFlag = 0
PseudocolorAtts.pointSize = 0.05
PseudocolorAtts.pointType = PseudocolorAtts.Point
PseudocolorAtts.pointSizeVarEnabled = 0
PseudocolorAtts.pointSizeVar = "default"
PseudocolorAtts.pointSizePixels = 2
PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID
PseudocolorAtts.lineType = PseudocolorAtts.Line
PseudocolorAtts.lineWidth = 0
PseudocolorAtts.tubeResolution = 10
PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox
PseudocolorAtts.tubeRadiusAbsolute = 0.125
PseudocolorAtts.tubeRadiusBBox = 0.005
PseudocolorAtts.tubeRadiusVarEnabled = 0
PseudocolorAtts.tubeRadiusVar = ""
PseudocolorAtts.tubeRadiusVarRatio = 10
PseudocolorAtts.tailStyle = PseudocolorAtts.None
PseudocolorAtts.headStyle = PseudocolorAtts.None
PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox
PseudocolorAtts.endPointRadiusAbsolute = 0.125
PseudocolorAtts.endPointRadiusBBox = 0.05
PseudocolorAtts.endPointResolution = 10
PseudocolorAtts.endPointRatio = 5
PseudocolorAtts.endPointRadiusVarEnabled = 0
PseudocolorAtts.endPointRadiusVar = ""
PseudocolorAtts.endPointRadiusVarRatio = 10
PseudocolorAtts.renderSurfaces = 1
PseudocolorAtts.renderWireframe = 0
PseudocolorAtts.renderPoints = 0
PseudocolorAtts.smoothingLevel = 0
PseudocolorAtts.legendFlag = 1
PseudocolorAtts.lightingFlag = 1
PseudocolorAtts.wireframeColor = (0, 0, 0, 0)
PseudocolorAtts.pointColor = (0, 0, 0, 0)
SetPlotOptions(PseudocolorAtts)


AddOperator("Transform", 0)
TransformAtts = TransformAttributes()
TransformAtts.doRotate = 0
TransformAtts.doScale = 1
TransformAtts.scaleOrigin = (0, 0, 0)
TransformAtts.scaleX = 3.24e-19
TransformAtts.scaleY = 3.24e-19
TransformAtts.scaleZ = 3.24e-19
TransformAtts.doTranslate = 0
TransformAtts.transformType = TransformAtts.Similarity  # Similarity, Coordinate, Linear
TransformAtts.inputCoordSys = TransformAtts.Cartesian  # Cartesian, Cylindrical, Spherical
TransformAtts.outputCoordSys = TransformAtts.Spherical  # Cartesian, Cylindrical, Spherical
TransformAtts.continuousPhi = 0
TransformAtts.invertLinearTransform = 0
TransformAtts.vectorTransformMethod = TransformAtts.AsDirection  # None, AsPoint, AsDisplacement, AsDirection
TransformAtts.transformVectors = 1
SetOperatorOptions(TransformAtts, 0)


ActivateDatabase(DB_level00)
AddPlot("Subset", "domains", 1, 0)
ActivateDatabase(DB_level01)
AddPlot("Subset", "domains", 1, 0)

SetActivePlots((2, 3))
AddOperator("Transform", 0)
SetOperatorOptions(TransformAtts, 0)

SetActivePlots((2, 3))
SubsetAtts = SubsetAttributes()
SubsetAtts.colorType = SubsetAtts.ColorBySingleColor
SubsetAtts.colorTableName = "Default"
SubsetAtts.invertColorTable = 0
SubsetAtts.legendFlag = 1
SubsetAtts.lineStyle = SubsetAtts.SOLID
SubsetAtts.lineWidth = 0
SubsetAtts.singleColor = (255, 255, 255, 255)
SubsetAtts.SetMultiColor(0, (255, 0, 0, 255))
SubsetAtts.subsetNames = ("0")
SubsetAtts.opacity = 1.0
SubsetAtts.wireframe = 1
SubsetAtts.drawInternal = 0
SubsetAtts.smoothingLevel = 0
SubsetAtts.pointSize = 0.05
SubsetAtts.pointType = SubsetAtts.Point
SubsetAtts.pointSizeVarEnabled = 0
SubsetAtts.pointSizeVar = "default"
SubsetAtts.pointSizePixels = 2
SetPlotOptions(SubsetAtts)


SetActivePlots((0,2))
AddOperator("Reflect", 0)
ReflectAtts = ReflectAttributes()
ReflectAtts.octant = ReflectAtts.PXPYPZ
ReflectAtts.useXBoundary = 1
ReflectAtts.specifiedX = 0
ReflectAtts.useYBoundary = 1
ReflectAtts.specifiedY = 0
ReflectAtts.useZBoundary = 1
ReflectAtts.specifiedZ = 0
ReflectAtts.reflections = (1, 1, 1, 1, 0, 0, 0, 0)
SetOperatorOptions(ReflectAtts, 0)

SetActivePlots((1,3))
AddOperator("Reflect", 0)
ReflectAtts = ReflectAttributes()
ReflectAtts.octant = ReflectAtts.PXPYPZ
ReflectAtts.useXBoundary = 1
ReflectAtts.specifiedX = 0
ReflectAtts.useYBoundary = 1
ReflectAtts.specifiedY = 0
ReflectAtts.useZBoundary = 1
ReflectAtts.specifiedZ = 0
ReflectAtts.reflections = (1, 0, 1, 1, 0, 0, 0, 0)
SetOperatorOptions(ReflectAtts, 0)


DrawPlots()
OpenGUI()
exit()

