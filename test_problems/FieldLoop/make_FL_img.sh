#!/bin/bash
#
# Make figures for the Field Loop test problem, using VisIt 1.12.0 and 
# the ImageMagick command-line tool 'convert'
#
#  - JM 2009-12-15: Added cropping of figures; plus annotations in python files.
#  - JM 2009-12-14: Writing file; can make Magnetic pressure and Current Density 
#    plots now.
#  - JM 2010.10.19: Changed so they are contour plots, not pseudocolour
# - 2016.03.06 JM: updated to fix a couple of bugs.

data_dir=$1
infile=$2
outfile=$3
visit_cmd=$4

#
# First plot magnetic pressure
#
cat <<EOF > plot_FL_results.py
OpenComputeEngine("localhost", ("serial"))
SetWindowArea(0,0,1200,800)

DeleteAllPlots()

OpenDatabase("localhost:${data_dir}/${infile}.*.silo database", 0)

SetDatabaseCorrelationOptions(0,0)

DefineScalarExpression("S_MagP","0.5*(sqr(MagneticFieldX)+sqr(MagneticFieldY)+sqr(MagneticFieldZ))")
#AddPlot("Pseudocolor", "S_MagP", 1, 1)

#P = PseudocolorAttributes()
#P.minFlag = 1
#P.maxFlag = 1
#P.min = 0
#P.max = 5.7e-7
#P.scaling = P.Linear
#P.colorTableName = "gray"
#SetPlotOptions(P)

AddPlot("Contour", "S_MagP", 1, 1)
ca = ContourAttributes()
ca.contourNLevels = 15
ca.colorType = ca.ColorBySingleColor # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
ca.singleColor = (0,0,0,255)
ca.contourMethod = ca.Level
ca.minFlag=1
ca.maxFlag=1
ca.min=0.0
ca.max=5.5e-7
ca.scaling = ca.Linear
ca.lineWidth=1
SetPlotOptions(ca)

va=GetView2D()
va.windowCoords = (-0.5,0.5, -0.5, 0.5)
SetView2D(va)

A = GetAnnotationAttributes()
A.axes2D.visible = 1
A.axes2D.autoSetTicks = 0
A.axes2D.autoSetScaling = 0
A.axes2D.lineWidth = 2
A.axes2D.tickLocation = A.axes2D.Outside  # Inside, Outside, Both
A.axes2D.tickAxes = A.axes2D.Off  #A.axes2D.BottomLeft  # Off, Bottom, Left, BottomLeft, All
A.axes2D.xAxis.title.visible = 0
A.axes2D.xAxis.label.visible = 1
A.axes2D.xAxis.label.font.font = A.axes2D.xAxis.label.font.Times  # Arial, Courier, Times
A.axes2D.xAxis.label.font.scale = 1.5
A.axes2D.xAxis.label.font.useForegroundColor = 1
A.axes2D.xAxis.label.font.color = (0, 0, 0, 255)
A.axes2D.xAxis.label.font.bold = 0
A.axes2D.xAxis.label.font.italic = 0
A.axes2D.xAxis.label.scaling = 0
A.axes2D.xAxis.tickMarks.visible = 1
A.axes2D.xAxis.tickMarks.majorMinimum = -0.5
A.axes2D.xAxis.tickMarks.majorMaximum = 0.5
A.axes2D.xAxis.tickMarks.minorSpacing = 0.1
A.axes2D.xAxis.tickMarks.majorSpacing = 0.5
A.axes2D.xAxis.grid = 0
A.axes2D.yAxis.title.visible = 0
A.axes2D.yAxis.label.visible = 1
A.axes2D.yAxis.label.font.font = A.axes2D.yAxis.label.font.Times  # Arial, Courier, Times
A.axes2D.yAxis.label.font.scale = 1.5
A.axes2D.yAxis.label.font.useForegroundColor = 1
A.axes2D.yAxis.label.font.color = (0, 0, 0, 255)
A.axes2D.yAxis.label.font.bold = 0
A.axes2D.yAxis.label.font.italic = 0
A.axes2D.yAxis.label.scaling = 0
A.axes2D.yAxis.tickMarks.visible = 1
A.axes2D.yAxis.tickMarks.majorMinimum = -0.5
A.axes2D.yAxis.tickMarks.majorMaximum = 0.5
A.axes2D.yAxis.tickMarks.minorSpacing = 0.1
A.axes2D.yAxis.tickMarks.majorSpacing = 0.5
A.axes2D.yAxis.grid = 0
A.userInfoFlag = 0
A.databaseInfoFlag = 0
#A.gradientColor1 = (0, 0, 255, 255)
#A.gradientColor2 = (0, 255, 255, 255)
A.backgroundMode = A.Solid
SetAnnotationAttributes(A)

tkyr=0
txtlabel=CreateAnnotationObject("Text2D")
txtlabel.visible = 1
txtlabel.active = 1
txtlabel.position = 0.215, 0.9
#txtlabel.width = 0.14
txtlabel.height = 0.04
txtlabel.textColor = (0, 0, 0, 255)
txtlabel.useForegroundForTextColor = 1
txtlabel.text = "Time =  %1.3f" %tkyr
txtlabel.fontFamily = txtlabel.Times  # Arial, Courier, Times
txtlabel.fontBold = 0
txtlabel.fontItalic = 0
txtlabel.fontShadow = 0

DrawPlots()

###################
## Make a Legend ##
###################
plotName = GetPlotList().GetPlots(0).plotName 
legend = GetAnnotationObject(plotName)
## See if we can scale the legend.
#legend.xScale = 1.0
#legend.yScale = 0.5
## the bounding box (draw=0 no, =1 yes)
#legend.drawBoundingBox = 0
#legend.boundingBoxColor = (200,200,200,230)
## number format
legend.numberFormat = "%1.2e"
## the font.
legend.fontFamily = legend.Times
#legend.fontBold = 1
#legend.fontItalic = 0
legend.fontHeight = 0.04
## turning off the title.
legend.drawTitle = 0
## Make it horizontal
#legend.orientation = legend.VerticalLeft
## moving the legend
legend.managePosition = 0
legend.position = (0.5,1.2)
#legend.position = (0.75,1.575)
## text color
##InvertBackgroundColor()
##legend.useForegroundForTextColor = 0
##legend.textColor = (255, 0, 0, 255)
## turning off the labels.
#legend.drawLabels = 1
#legend.drawMinMax = 1
#legend.numTicks = 8

sw = SaveWindowAttributes()
sw.fileName = "test"
sw.width = 1200
sw.height= 700
sw.outputToCurrentDirectory = 0
sw.outputDirectory = "${data_dir}"
sw.family = 0
sw.format=sw.PNG
sw.resConstraint = sw.NoConstraint



for state in range(TimeSliderGetNStates()):
	SetTimeSliderState(state)
        Query("Time")
        ttt=GetQueryOutputValue()
        tkyr=ttt/1.0
        stst = "Time =  %1.3f" %tkyr
	txtlabel.text = stst
	sw.fileName = "%s_%02d" %("${outfile}_MagP", state)
	SetSaveWindowAttributes(sw)
	SaveWindow()
quit()
EOF

$visit_cmd -cli -nowin -s plot_FL_results.py
sleep 15
#exit

#
# Now plot Curl(B) = current density.
#
cat <<EOF > plot_FL_results.py
OpenComputeEngine("localhost", ("serial"))
SetWindowArea(0,0,1200,800)

DeleteAllPlots()

OpenDatabase("localhost:${data_dir}/${infile}.*.silo database", 0)

SetDatabaseCorrelationOptions(0,0)

DefineVectorExpression("S_B2D","{MagneticFieldX,MagneticFieldY}")
DefineScalarExpression("S_CurlB2D","curl(S_B2D)")
DefineScalarExpression("S_MagP", "0.5*(sqr(MagneticFieldX)+sqr(MagneticFieldY)+sqr(MagneticFieldZ))")

#AddPlot("Pseudocolor", "S_CurlB2D", 1, 1)
#P = PseudocolorAttributes()
#P.minFlag = 0
#P.min = 0.0
#P.maxFlag = 0
#P.max = 0.1
#P.scaling = P.Linear
#P.colorTableName = "gray"
#SetPlotOptions(P)

AddPlot("Contour", "S_CurlB2D", 1, 1)
ca = ContourAttributes()
ca.contourNLevels = 15
ca.colorType = ca.ColorBySingleColor # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
ca.singleColor = (0,0,0,255)
ca.contourMethod = ca.Level
ca.minFlag=0
ca.maxFlag=0
ca.min=0.0
ca.max=1.0
ca.scaling = ca.Linear
ca.lineWidth=1
SetPlotOptions(ca)

va=GetView2D()
va.windowCoords = (-0.5,0.5, -0.5, 0.5)
SetView2D(va)

A = GetAnnotationAttributes()
A.axes2D.visible = 1
A.axes2D.autoSetTicks = 0
A.axes2D.autoSetScaling = 0
A.axes2D.lineWidth = 2
A.axes2D.tickLocation = A.axes2D.Outside  # Inside, Outside, Both
A.axes2D.tickAxes = A.axes2D.Off  #A.axes2D.BottomLeft  # Off, Bottom, Left, BottomLeft, All
A.axes2D.xAxis.title.visible = 0
A.axes2D.xAxis.label.visible = 1
A.axes2D.xAxis.label.font.font = A.axes2D.xAxis.label.font.Times  # Arial, Courier, Times
A.axes2D.xAxis.label.font.scale = 1.5
A.axes2D.xAxis.label.font.useForegroundColor = 1
A.axes2D.xAxis.label.font.color = (0, 0, 0, 255)
A.axes2D.xAxis.label.font.bold = 0
A.axes2D.xAxis.label.font.italic = 0
A.axes2D.xAxis.label.scaling = 0
A.axes2D.xAxis.tickMarks.visible = 1
A.axes2D.xAxis.tickMarks.majorMinimum = -0.5
A.axes2D.xAxis.tickMarks.majorMaximum = 0.5
A.axes2D.xAxis.tickMarks.minorSpacing = 0.1
A.axes2D.xAxis.tickMarks.majorSpacing = 0.5
A.axes2D.xAxis.grid = 0
A.axes2D.yAxis.title.visible = 0
A.axes2D.yAxis.label.visible = 1
A.axes2D.yAxis.label.font.font = A.axes2D.yAxis.label.font.Times  # Arial, Courier, Times
A.axes2D.yAxis.label.font.scale = 1.5
A.axes2D.yAxis.label.font.useForegroundColor = 1
A.axes2D.yAxis.label.font.color = (0, 0, 0, 255)
A.axes2D.yAxis.label.font.bold = 0
A.axes2D.yAxis.label.font.italic = 0
A.axes2D.yAxis.label.scaling = 0
A.axes2D.yAxis.tickMarks.visible = 1
A.axes2D.yAxis.tickMarks.majorMinimum = -0.5
A.axes2D.yAxis.tickMarks.majorMaximum = 0.5
A.axes2D.yAxis.tickMarks.minorSpacing = 0.1
A.axes2D.yAxis.tickMarks.majorSpacing = 0.5
A.axes2D.yAxis.grid = 0
A.userInfoFlag = 0
A.databaseInfoFlag = 0
#A.gradientColor1 = (0, 0, 255, 255)
#A.gradientColor2 = (0, 255, 255, 255)
A.backgroundMode = A.Solid
SetAnnotationAttributes(A)

tkyr=0
txtlabel=CreateAnnotationObject("Text2D")
txtlabel.visible = 1
txtlabel.active = 1
txtlabel.position = 0.215, 0.9
#txtlabel.width = 0.14
txtlabel.height = 0.04
txtlabel.textColor = (0, 0, 0, 255)
txtlabel.useForegroundForTextColor = 1
txtlabel.text = "Time =  %1.3f" %tkyr
txtlabel.fontFamily = txtlabel.Times  # Arial, Courier, Times
txtlabel.fontBold = 0
txtlabel.fontItalic = 0
txtlabel.fontShadow = 0

DrawPlots()

###################
## Make a Legend ##
###################
plotName = GetPlotList().GetPlots(0).plotName 
legend = GetAnnotationObject(plotName)
## See if we can scale the legend.
#legend.xScale = 1.0
#legend.yScale = 0.5
## the bounding box (draw=0 no, =1 yes)
#legend.drawBoundingBox = 0
#legend.boundingBoxColor = (200,200,200,230)
## number format
legend.numberFormat = "%1.2e"
## the font.
legend.fontFamily = legend.Times
#legend.fontBold = 1
#legend.fontItalic = 0
legend.fontHeight = 0.04
## turning off the title.
legend.drawTitle = 0
## Make it horizontal
#legend.orientation = legend.VerticalLeft
## moving the legend
legend.managePosition = 0
legend.position = (0.5,1.2)
#legend.position = (0.75,1.575)
## text color
##InvertBackgroundColor()
##legend.useForegroundForTextColor = 0
##legend.textColor = (255, 0, 0, 255)
## turning off the labels.
#legend.drawLabels = 1
#legend.drawMinMax = 1
#legend.numTicks = 8

sw = SaveWindowAttributes()
sw.fileName = "test"
sw.width = 1200
sw.height= 700
sw.outputToCurrentDirectory = 0
sw.outputDirectory = "${data_dir}"
sw.family = 0
sw.format=sw.PNG
sw.resConstraint = sw.NoConstraint



for state in range(TimeSliderGetNStates()):
	SetTimeSliderState(state)
        Query("Time")
        ttt=GetQueryOutputValue()
        tkyr=ttt/1.0
        stst = "Time =  %1.3f" %tkyr
	txtlabel.text = stst
	sw.fileName = "%s_%02d" %("${outfile}_CurlB2D", state)
	SetSaveWindowAttributes(sw)
	SaveWindow()
quit()
EOF

$visit_cmd -cli -nowin -s plot_FL_results.py
sleep 15

#
# Now resize the images
#
convert -crop 635x600+180+25 ${data_dir}/${outfile}_CurlB2D_00.png ${data_dir}/${outfile}_CurlB2D_00.png
convert -crop 635x600+180+25 ${data_dir}/${outfile}_MagP_00.png    ${data_dir}/${outfile}_MagP_00.png
convert -crop 635x600+180+25 ${data_dir}/${outfile}_CurlB2D_01.png ${data_dir}/${outfile}_CurlB2D_01.png
convert -crop 635x600+180+25 ${data_dir}/${outfile}_MagP_01.png    ${data_dir}/${outfile}_MagP_01.png
convert -crop 635x600+180+25 ${data_dir}/${outfile}_CurlB2D_02.png ${data_dir}/${outfile}_CurlB2D_02.png
convert -crop 635x600+180+25 ${data_dir}/${outfile}_MagP_02.png    ${data_dir}/${outfile}_MagP_02.png
#
# Old pseudocolour plots:
#convert -crop 452x453+464+144 ${data_dir}/${outfile}_CurlB2D_00.png ${data_dir}/${outfile}_CurlB2D_00.png
#convert -crop 452x453+464+144 ${data_dir}/${outfile}_MagP_00.png    ${data_dir}/${outfile}_MagP_00.png
#convert -crop 452x453+464+144 ${data_dir}/${outfile}_CurlB2D_01.png ${data_dir}/${outfile}_CurlB2D_01.png
#convert -crop 452x453+464+144 ${data_dir}/${outfile}_MagP_01.png    ${data_dir}/${outfile}_MagP_01.png
#convert -crop 452x453+464+144 ${data_dir}/${outfile}_CurlB2D_02.png ${data_dir}/${outfile}_CurlB2D_02.png
#convert -crop 452x453+464+144 ${data_dir}/${outfile}_MagP_02.png    ${data_dir}/${outfile}_MagP_02.png

exit
