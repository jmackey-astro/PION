#!/bin/bash
#
# 2013.06.14 JM: Make plots for Pascal Tremblin's test problem.

FDIR=$1
FBASE=$2
SAVEDIR=$3
SAVEFILE=$4
SIM_TITLE=$5
XMIN=$6
XMAX=$7
YMIN=$8
YMAX=$9
N_MIN=${10}
N_MAX=${11}
VEL=${12}    # Velocity of grid with respect to star.
ANGLE=${13}  # ANGLE in degrees of velocity wrt -ve x-axis.

mkdir $SAVEDIR
VISITFILE=visit_${FBASE}_NVY.py
#VISITCMD=/vol/aibn129/aibn129_1/jmackey/extra_libraries/visit_bin/bin/visit
VISITCMD=/Applications/VisIt.app/Contents/MacOS/VisIt

cat <<EOF > ${VISITFILE}
# Define expressions
# First peculiar velocity in km/s
DefineVectorExpression("dv2d", "{1.0e-5*{VelocityX+${VEL}*cos(0.0174533*${ANGLE}),VelocityY+${VEL}*sin(0.0174533*${ANGLE})}}")
DefineScalarExpression("dv2d_mag", "magnitude(dv2d)")
DefineScalarExpression("NH", "Density/1.67e-24")
DefineScalarExpression("LogNH", "log10(NH)")

OpenComputeEngine("localhost", ("-nn", "1", "-np", "2"))
SetWindowArea(0,0,800,800)

DeleteAllPlots()

SetDatabaseCorrelationOptions(2,0)

# OpenDatabase("hostname:path/to/files",timestep)
#OpenDatabase("localhost:${FDIR}/${FBASE}.*.silo database",0)
OpenDatabase("${FDIR}/${FBASE}.*.silo database",0)

AddPlot("Pseudocolor", "LogNH", 1, 0)
PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.legendFlag = 1
PseudocolorAtts.lightingFlag = 0
PseudocolorAtts.minFlag = 1
PseudocolorAtts.maxFlag = 1
PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
PseudocolorAtts.min = ${N_MIN}
PseudocolorAtts.max = ${N_MAX}
PseudocolorAtts.pointSize = 0.05
PseudocolorAtts.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Point, Sphere
PseudocolorAtts.skewFactor = 1
PseudocolorAtts.opacity = 1
PseudocolorAtts.colorTableName = "caleblack"
PseudocolorAtts.invertColorTable = 0
PseudocolorAtts.smoothingLevel = 0
PseudocolorAtts.pointSizeVarEnabled = 0
PseudocolorAtts.pointSizeVar = "default"
PseudocolorAtts.pointSizePixels = 2
PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
PseudocolorAtts.lineWidth = 0
PseudocolorAtts.opacityType = PseudocolorAtts.Explicit  # Explicit, ColorTable
SetPlotOptions(PseudocolorAtts)

AddPlot("Contour", "Tr000_H1p___", 1, 1)
ContourAtts = ContourAttributes()
ContourAtts.defaultPalette.smoothing = ContourAtts.defaultPalette.None  # None, Linear, CubicSpline
ContourAtts.defaultPalette.equalSpacingFlag = 1
ContourAtts.defaultPalette.discreteFlag = 1
ContourAtts.defaultPalette.externalFlag = 0
ContourAtts.changedColors = (0,1,2,3,4, 5)
ContourAtts.colorType = ContourAtts.ColorByMultipleColors  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
ContourAtts.colorTableName = "Default"
ContourAtts.invertColorTable = 0
ContourAtts.legendFlag = 1
ContourAtts.lineStyle = ContourAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
ContourAtts.lineWidth = 5
ContourAtts.singleColor = (0, 255, 0, 255)
ContourAtts.SetMultiColor(0, (255, 255, 255, 255))
ContourAtts.SetMultiColor(1, (255, 255, 255, 235))
ContourAtts.SetMultiColor(2, (255, 255, 255, 215))
ContourAtts.SetMultiColor(3, (255, 255, 255, 195))
ContourAtts.SetMultiColor(4, (255, 255, 255, 175))
ContourAtts.SetMultiColor(5, (255, 255, 255, 155))
ContourAtts.SetMultiColor(6, (255, 255, 255, 135))
ContourAtts.SetMultiColor(7, (255, 255, 255, 115))
ContourAtts.SetMultiColor(8, (255, 255, 255,  95))
ContourAtts.SetMultiColor(9, (255, 68, 68, 255))
ContourAtts.contourNLevels = 5
ContourAtts.contourValue = ()
ContourAtts.contourPercent = ()
ContourAtts.contourMethod = ContourAtts.Level  # Level, Value, Percent
ContourAtts.minFlag = 1
ContourAtts.maxFlag = 1
ContourAtts.min = 0.01
ContourAtts.max = 0.99
ContourAtts.scaling = ContourAtts.Linear  # Linear, Log
ContourAtts.wireframe = 0
SetPlotOptions(ContourAtts)


##########################################################################################
##### TRANSFORM BOTH PLOTS TO PARSECS, ROTATED SO THAT MOTION IS IN X-DIRECTION  #########
##########################################################################################
SetActivePlots((0,1))
AddOperator("Transform",1)

TransformAtts = TransformAttributes()
TransformAtts.doRotate = 1
TransformAtts.rotateOrigin = (0, 0, 0)
TransformAtts.rotateAxis = (0, 0, 1)
TransformAtts.rotateAmount = -${ANGLE}
TransformAtts.rotateType = TransformAtts.Deg  # Deg, Rad
TransformAtts.doScale = 1
TransformAtts.scaleOrigin = (0, 0, 0)
TransformAtts.scaleX = 3.24e-19
TransformAtts.scaleY = 3.24e-19
TransformAtts.scaleZ = 3.24e-19
TransformAtts.doTranslate = 0
TransformAtts.translateX = 0
TransformAtts.translateY = 0
TransformAtts.translateZ = 0
TransformAtts.transformType = TransformAtts.Similarity  # Similarity, Coordinate, Linear
TransformAtts.inputCoordSys = TransformAtts.Cartesian  # Cartesian, Cylindrical, Spherical
TransformAtts.outputCoordSys = TransformAtts.Spherical  # Cartesian, Cylindrical, Spherical
TransformAtts.m00 = 1
TransformAtts.m01 = 0
TransformAtts.m02 = 0
TransformAtts.m03 = 0
TransformAtts.m10 = 0
TransformAtts.m11 = 1
TransformAtts.m12 = 0
TransformAtts.m13 = 0
TransformAtts.m20 = 0
TransformAtts.m21 = 0
TransformAtts.m22 = 1
TransformAtts.m23 = 0
TransformAtts.m30 = 0
TransformAtts.m31 = 0
TransformAtts.m32 = 0
TransformAtts.m33 = 1
TransformAtts.invertLinearTransform = 0
TransformAtts.vectorTransformMethod = TransformAtts.AsDirection  # None, AsPoint, AsDisplacement, AsDirection
TransformAtts.transformVectors = 1
SetOperatorOptions(TransformAtts,0,1)

##########################################################################################
####### SET SOME ANNOTATION ATTRIBUTES
##########################################################################################

A = AnnotationAttributes()
A.axes2D.visible = 1
A.axes2D.autoSetTicks = 0
A.axes2D.autoSetScaling = 0
A.axes2D.lineWidth = 4
A.axes2D.tickLocation = A.axes2D.Both  # Inside, Outside, Both
A.axes2D.tickAxes = A.axes2D.BottomLeft  # Off, Bottom, Left, BottomLeft, All
A.axes2D.xAxis.title.visible = 0
A.axes2D.xAxis.title.font.font = A.axes2D.xAxis.title.font.Times  # Arial, Courier, Times
A.axes2D.xAxis.title.font.scale = 1
A.axes2D.xAxis.title.font.useForegroundColor = 1
A.axes2D.xAxis.title.font.color = (0, 0, 0, 255)
A.axes2D.xAxis.title.font.bold = 1
A.axes2D.xAxis.title.font.italic = 1
A.axes2D.xAxis.title.userTitle = 0
A.axes2D.xAxis.title.userUnits = 0
A.axes2D.xAxis.title.title = ""
A.axes2D.xAxis.title.units = ""
A.axes2D.xAxis.label.visible = 1
A.axes2D.xAxis.label.font.font = A.axes2D.xAxis.label.font.Times  # Arial, Courier, Times
A.axes2D.xAxis.label.font.scale = 2
A.axes2D.xAxis.label.font.useForegroundColor = 1
A.axes2D.xAxis.label.font.color = (0, 0, 0, 255)
A.axes2D.xAxis.label.font.bold = 1
A.axes2D.xAxis.label.font.italic = 0
A.axes2D.xAxis.label.scaling = 0
A.axes2D.xAxis.tickMarks.visible = 1
A.axes2D.xAxis.tickMarks.majorMinimum = -100
A.axes2D.xAxis.tickMarks.majorMaximum = 100
A.axes2D.xAxis.tickMarks.minorSpacing = 0.1
A.axes2D.xAxis.tickMarks.majorSpacing = 1.0
A.axes2D.xAxis.grid = 0
A.axes2D.yAxis.title.visible = 0
A.axes2D.yAxis.title.font.font = A.axes2D.yAxis.title.font.Times  # Arial, Courier, Times
A.axes2D.yAxis.title.font.scale = 1
A.axes2D.yAxis.title.font.useForegroundColor = 1
A.axes2D.yAxis.title.font.color = (0, 0, 0, 255)
A.axes2D.yAxis.title.font.bold = 1
A.axes2D.yAxis.title.font.italic = 1
A.axes2D.yAxis.title.userTitle = 0
A.axes2D.yAxis.title.userUnits = 0
A.axes2D.yAxis.title.title = ""
A.axes2D.yAxis.title.units = ""
A.axes2D.yAxis.label.visible = 1
A.axes2D.yAxis.label.font.font = A.axes2D.yAxis.label.font.Times  # Arial, Courier, Times
A.axes2D.yAxis.label.font.scale = 2
A.axes2D.yAxis.label.font.useForegroundColor = 1
A.axes2D.yAxis.label.font.color = (0, 0, 0, 255)
A.axes2D.yAxis.label.font.bold = 1
A.axes2D.yAxis.label.font.italic = 0
A.axes2D.yAxis.label.scaling = 0
A.axes2D.yAxis.tickMarks.visible = 1
A.axes2D.yAxis.tickMarks.majorMinimum = -100
A.axes2D.yAxis.tickMarks.majorMaximum = 100
A.axes2D.yAxis.tickMarks.minorSpacing = 0.1
A.axes2D.yAxis.tickMarks.majorSpacing = 1.0
A.axes2D.yAxis.grid = 0
A.userInfoFlag = 0
A.userInfoFont.font = A.userInfoFont.Times  # Arial, Courier, Times
A.userInfoFont.scale = 1
A.userInfoFont.useForegroundColor = 1
A.userInfoFont.color = (0, 0, 0, 255)
A.userInfoFont.bold = 0
A.userInfoFont.italic = 0
A.databaseInfoFlag = 0
A.databaseInfoFont.font = A.databaseInfoFont.Times  # Arial, Courier, Times
A.databaseInfoFont.scale = 1
A.databaseInfoFont.useForegroundColor = 1
A.databaseInfoFont.color = (0, 0, 0, 255)
A.databaseInfoFont.bold = 0
A.databaseInfoFont.italic = 0
A.databaseInfoExpansionMode = A.File  # File, Directory, Full, Smart, SmartDirectory
A.databaseInfoTimeScale = 3.1645570e-11
A.databaseInfoTimeOffset = 0
A.legendInfoFlag = 1
A.backgroundColor = (255, 255, 255, 255)
A.foregroundColor = (0, 0, 0, 255)
A.backgroundMode = A.Solid  # Solid, Gradient, Image, ImageSphere
A.axesArray.visible = 1
A.axesArray.ticksVisible = 1
A.axesArray.autoSetTicks = 1
A.axesArray.autoSetScaling = 1
A.axesArray.lineWidth = 4
A.axesArray.axes.title.visible = 1
A.axesArray.axes.title.font.font = A.axesArray.axes.title.font.Times  # Arial, Courier, Times
A.axesArray.axes.title.font.scale = 1
A.axesArray.axes.title.font.useForegroundColor = 1
A.axesArray.axes.title.font.color = (0, 0, 0, 255)
A.axesArray.axes.title.font.bold = 0
A.axesArray.axes.title.font.italic = 0
A.axesArray.axes.title.userTitle = 0
A.axesArray.axes.title.userUnits = 0
A.axesArray.axes.title.title = ""
A.axesArray.axes.title.units = ""
A.axesArray.axes.label.visible = 1
A.axesArray.axes.label.font.font = A.axesArray.axes.label.font.Times  # Arial, Courier, Times
A.axesArray.axes.label.font.scale = 1
A.axesArray.axes.label.font.useForegroundColor = 1
A.axesArray.axes.label.font.color = (0, 0, 0, 255)
A.axesArray.axes.label.font.bold = 0
A.axesArray.axes.label.font.italic = 0
A.axesArray.axes.label.scaling = 0
A.axesArray.axes.tickMarks.visible = 1
A.axesArray.axes.tickMarks.majorMinimum = 0
A.axesArray.axes.tickMarks.majorMaximum = 1
A.axesArray.axes.tickMarks.minorSpacing = 0.02
A.axesArray.axes.tickMarks.majorSpacing = 0.2
A.axesArray.axes.grid = 0
SetAnnotationAttributes(A)
##########################################################################################
##########################################################################################



##########################################################################################
##########  GET THE VIEW AND IMAGE RESOLUTION THAT WE WANT ##############
##########################################################################################
DrawPlots()
ResizeWindow(1,1876,1536)

sw = SaveWindowAttributes()
sw.fileName = "test"
sw.width = 1876
sw.height= 1536
sw.outputToCurrentDirectory = 0
sw.outputDirectory = "${SAVEDIR}"
sw.family = 0
sw.format = sw.PNG # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK
sw.resConstraint = sw.NoConstraint

View2DAtts = View2DAttributes()
View2DAtts.windowCoords = (${XMIN}, ${XMAX}, ${YMIN}, ${YMAX})
View2DAtts.viewportCoords = (0.075, 0.9, 0.075, 0.98)
View2DAtts.fullFrameActivationMode = View2DAtts.Off  # On, Off, Auto
View2DAtts.fullFrameAutoThreshold = 100
View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
SetView2D(View2DAtts)
##########################################################################################

##########################################################################################
# Add a time slider in the lower right corner
slider = CreateAnnotationObject("TimeSlider")
slider.height = 0.1
slider.width = 0.4
slider.text = "Time=\$time kyr"
slider.position = (0.55, 0.1)
# Make the slider transparent, so that only the time shows.
slider.startColor = (0, 255, 255, 0)
slider.endColor = (255, 255, 255, 0)
##########################################################################################


##########################################################################################
# Incorporate AIfA logo image (aifa_logo.jpeg) as an annotation
#image = CreateAnnotationObject("Image")
#image.image = "/vol/aibn129/aibn129_1/jmackey/active/Talks/Talks_figs/Logos/aifa_logo.jpeg"
#image.position = (0.002, 0.002)
#image.width = 15
#image.height = 15
#image.maintainAspectRatio = 1
#image.transparencyColor = (255,255,255, 255)
#image.useTransparencyColor = 0
##########################################################################################


##########################################################################################
# Add Text label describing simulation
text=CreateAnnotationObject("Text2D")
text.visible = 1
text.active = 1
text.position = (0.6, 0.94)
text.height = 0.03
text.textColor = (0, 0, 0, 255)
text.useForegroundForTextColor = 1
text.text = "${SIM_TITLE}"
text.fontFamily = text.Times  # Arial, Courier, Times
text.fontBold = 1
text.fontItalic = 0
text.fontShadow = 0
##########################################################################################
##########################################################################################
# Add Text label claiming copyright!
#text2=CreateAnnotationObject("Text2D")
#text2.visible = 1
#text2.active = 1
#text2.position = (0.77, 0.005)
#text2.width = 0.23
#text2.textColor = (0, 0, 0, 255)
#text2.useForegroundForTextColor = 1
#text2.text = "Jonathan Mackey, 2012"
#text2.fontFamily = text.Times  # Arial, Courier, Times
#text2.fontBold = 0
#text2.fontItalic = 0
#text2.fontShadow = 0
##########################################################################################

#####################################################################
##########  ADD LABELS TO THE AXES, AT THE POSITION WE WANT #########
#####################################################################
t1=CreateAnnotationObject("Text2D")
t1.visible = 1
t1.active = 1
t1.position = (0.475, 0.015)
t1.height = 0.03
t1.textColor = (0, 0, 0, 255)
t1.useForegroundForTextColor = 1
t1.text = "x1 (pc)"
t1.fontFamily = t1.Times  # Arial, Courier, Times
t1.fontBold = 1
t1.fontItalic = 0
t1.fontShadow = 0

t2=CreateAnnotationObject("Text2D")
t2.visible = 1
t2.active = 1
t2.position = (0.02, 0.65)
t2.height = 0.03
t2.textColor = (0, 0, 0, 255)
t2.useForegroundForTextColor = 1
t2.text = "x2"
t2.fontFamily = t2.Times  # Arial, Courier, Times
t2.fontBold = 1
t2.fontItalic = 0
t2.fontShadow = 0

t3=CreateAnnotationObject("Text2D")
t3.visible = 1
t3.active = 1
t3.position = (0.01, 0.6)
t3.height = 0.03
t3.textColor = (0, 0, 0, 255)
t3.useForegroundForTextColor = 1
t3.text = "(pc)"
t3.fontFamily = t2.Times  # Arial, Courier, Times
t3.fontBold = 1
t3.fontItalic = 0
t3.fontShadow = 0

##########################################################################################
## MAKE A NICE LEGEND ALONG THE RIGHT SIDE OF THE FIGURE
##########################################################################################
plotName = GetPlotList().GetPlots(0).plotName 
legend = GetAnnotationObject(plotName)
legend.xScale = 1
legend.yScale = 2.0
legend.drawBoundingBox = 0
legend.boundingBoxColor = (200,200,200,230)
legend.orientation = legend.VerticalRight
legend.managePosition = 0
legend.position = (0.83,0.98)
#InvertBackgroundColor()
#legend.useForegroundForTextColor = 0
#legend.textColor = (255, 0, 0, 255)
legend.numberFormat = "%1.1f"
legend.numTicks = 5
legend.fontFamily = legend.Times
legend.fontBold = 1
legend.fontItalic = 0
legend.fontHeight = 0.04
# turning off the labels.
legend.drawLabels = 1
legend.drawMinMax = 1
# turning off the title.
legend.drawTitle = 0

plotName = GetPlotList().GetPlots(1).plotName 
legend2 = GetAnnotationObject(plotName)
legend2.numTicks = 5
legend2.numberFormat = "%1.2f"
legend2.xScale = 1
legend2.yScale = 2
legend2.orientation = legend.VerticalRight
legend2.managePosition = 0
legend2.position = (0.83,0.4)
legend2.fontFamily = legend.Times
legend2.fontBold = 1
legend2.fontItalic = 0
legend2.fontHeight = 0.04
legend2.drawLabels = 1
legend2.drawMinMax = 1
legend2.drawTitle = 0
legend2.drawBoundingBox = 1

SetActivePlots(1)
SetPlotOptions(ContourAtts)
SetActivePlots((0,1))
##########################################################################################


##########################################################################################
##########################################################################################
##########################################################################################


for state in range(0,TimeSliderGetNStates(),1):
	SetTimeSliderState(state)
        sw.width = 1876
        sw.height= 1536
        sw.format = sw.PNG
        sw.fileName = "%s_%03d" %("${SAVEFILE}", state)
	SetSaveWindowAttributes(sw)
	SaveWindow()

quit()
EOF


${VISITCMD} -cli -nowin -s ${VISITFILE}
#${VISITCMD} -cli -s ${VISITFILE}


