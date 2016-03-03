#!/bin/bash
#
# Generate PNG figures of the final frame of DMR test problems.
#
# Modifications:
# - 2015.03.24 JM: edited to write png files.

# input a filename to read as $1
filedir=$1
filename=$2
visit_cmd=$3

cat <<EOF > plot_DMR_results.py
OpenComputeEngine("localhost", ("serial"))
SetWindowArea(0,0,1200,800)

DeleteAllPlots()

OpenDatabase("localhost:${filedir}/${filename}.*.silo database", 0)

SetDatabaseCorrelationOptions(0,0)

AddPlot("Contour", "Density", 1, 1)

ca = ContourAttributes()
ca.contourNLevels = 30
ca.singleColor = (0,0,0,255)
ca.contourMethod = ca.Level
ca.minFlag=1
ca.maxFlag=1
ca.min=1.4
ca.max=22.5
ca.scaling = ca.Linear
SetPlotOptions(ca)

#AddOperator("Transform", 0)
#T = TransformAttributes()
#T.doScale = 0
#T.scaleOrigin = (0, 0, 0)
#T.doTranslate = 1
#T.translateX = 0
#T.translateY = 0
#T.translateZ = 0
#SetOperatorOptions(T, 0)

#
# Now Annotations:
#
A = GetAnnotationAttributes()
A.axes2D.visible = 1
A.axes2D.autoSetTicks = 1
A.axes2D.autoSetScaling = 1
A.axes2D.lineWidth = 0
A.axes2D.tickLocation = A.axes2D.Outside  # Inside, Outside, Both
A.axes2D.tickAxes = A.axes2D.BottomLeft  # Off, Bottom, Left, BottomLeft, All
A.axes2D.xAxis.title.visible = 0
A.axes2D.xAxis.label.visible = 1
A.axes2D.xAxis.label.font.font = A.axes2D.xAxis.label.font.Courier  # Arial, Courier, Times
A.axes2D.xAxis.label.font.scale = 2
A.axes2D.xAxis.grid = 0
A.axes2D.yAxis.title.visible = 0
A.axes2D.yAxis.label.visible = 1
A.axes2D.yAxis.label.font.font = A.axes2D.yAxis.label.font.Courier  # Arial, Courier, Times
A.axes2D.yAxis.label.font.scale = 2
A.axes2D.yAxis.grid = 0
A.userInfoFlag = 0
A.databaseInfoFlag = 0
#A.gradientColor1 = (0, 0, 255, 255)
#A.gradientColor2 = (0, 255, 255, 255)
A.backgroundMode = A.Solid
A.legendInfoFlag=0

SetAnnotationAttributes(A)


DrawPlots()
#state = TimeSliderGetNStates()
#SetTimeSliderState(state-1)
SetTimeSliderState(4)

sw = SaveWindowAttributes()
sw.fileName = "${filename}"
sw.width = 1182
sw.height= 600
sw.outputToCurrentDirectory = 0
sw.outputDirectory = "${filedir}"
sw.family = 0
sw.format = sw.PNG # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK
sw.resConstraint = sw.NoConstraint
SetSaveWindowAttributes(sw)
SaveWindow()

quit()
EOF

$visit_cmd -cli -nowin -s plot_DMR_results.py

exit
