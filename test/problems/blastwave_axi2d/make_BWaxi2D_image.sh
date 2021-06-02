#!/bin/bash
#
# call with ./make_BWaxi2D_image.sh $test_dir/blastwave_axi2d $code_dir $data_dir $visit_cmd $infile_base $outfile_base

test_dir=$1
code_dir=$2
data_dir=$3
visit_cmd=$4
infile=$5
outfile=$6

#test_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test/problems/blastwave_axi2d
#code_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
#data_dir=/vol/aibn129/aibn129_1/jmackey/data_etc/code_tests/blastwave_axi2d
#visit_cmd=/vol/aibn129/aibn129_1/jmackey/extra_libraries/visit_bin/bin/visit
#infile=BWaxi2D_HalfPlane_n128_Roe_FKJav01
#outfile=BWaxi2D_HalfPlane_n128_Roe_FKJav01

#
# Plot Gas Density
#
cat <<EOF > plot_BWaxi2d_results.py
OpenComputeEngine("localhost", ("serial"))
SetWindowArea(0,0,1200,800)

DeleteAllPlots()

OpenDatabase("localhost:${data_dir}/${infile}.*.silo database", 0)

SetDatabaseCorrelationOptions(0,0)

DefineScalarExpression("ScaledDensity1e22","Density*1.0e+22")

#AddPlot("Pseudocolor", "S_MagP", 1, 1)
#P = PseudocolorAttributes()
#P.minFlag = 1
#P.maxFlag = 1
#P.min = 0
#P.max = 5.7e-7
#P.scaling = P.Linear
#P.colorTableName = "gray"
#SetPlotOptions(P)

AddPlot("Contour", "ScaledDensity1e22", 1, 1)
ca = ContourAttributes()
ca.contourNLevels = 17
ca.colorType = ca.ColorByMultipleColors # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
ca.singleColor = (0,0,0,255)
ca.contourMethod = ca.Level
ca.minFlag=1
ca.maxFlag=1
ca.min=1.0e-3
ca.max=1.0e+1
ca.scaling = ca.Log
ca.lineWidth=2
SetPlotOptions(ca)

AddOperator("Transform", 0)
TA = TransformAttributes()
TA.doRotate = 0
TA.doScale = 1
TA.scaleOrigin = (0, 0, 0)
TA.scaleX = 3.24e-19
TA.scaleY = 3.24e-19
TA.scaleZ = 3.24e-19
TA.doTranslate = 0
TA.translateX = 0
TA.translateY = 0
TA.translateZ = 0
SetOperatorOptions(TA, 1)

A = GetAnnotationAttributes()
A.axes2D.visible = 1
A.axes2D.autoSetTicks = 0
A.axes2D.autoSetScaling = 0
A.axes2D.lineWidth = 2
A.axes2D.tickLocation = A.axes2D.Outside  # Inside, Outside, Both
A.axes2D.tickAxes = A.axes2D.Off  #A.axes2D.BottomLeft  # Off, Bottom, Left, BottomLeft, All
A.axes2D.xAxis.title.visible = 1
A.axes2D.xAxis.label.visible = 1
A.axes2D.xAxis.label.font.font = A.axes2D.xAxis.label.font.Times  # Arial, Courier, Times
A.axes2D.xAxis.label.font.scale = 1.5
A.axes2D.xAxis.label.font.useForegroundColor = 1
A.axes2D.xAxis.label.font.color = (0, 0, 0, 255)
A.axes2D.xAxis.label.font.bold = 0
A.axes2D.xAxis.label.font.italic = 0
A.axes2D.xAxis.label.scaling = 0
A.axes2D.xAxis.tickMarks.visible = 1
A.axes2D.xAxis.tickMarks.majorMinimum = -10
A.axes2D.xAxis.tickMarks.majorMaximum = 10
A.axes2D.xAxis.tickMarks.minorSpacing = 0.5
A.axes2D.xAxis.tickMarks.majorSpacing = 2
A.axes2D.xAxis.grid = 0
A.axes2D.yAxis.title.visible = 1
A.axes2D.yAxis.label.visible = 1
A.axes2D.yAxis.label.font.font = A.axes2D.yAxis.label.font.Times  # Arial, Courier, Times
A.axes2D.yAxis.label.font.scale = 1.5
A.axes2D.yAxis.label.font.useForegroundColor = 1
A.axes2D.yAxis.label.font.color = (0, 0, 0, 255)
A.axes2D.yAxis.label.font.bold = 0
A.axes2D.yAxis.label.font.italic = 0
A.axes2D.yAxis.label.scaling = 0
A.axes2D.yAxis.tickMarks.visible = 1
A.axes2D.yAxis.tickMarks.majorMinimum = 0.0
A.axes2D.yAxis.tickMarks.majorMaximum = 10.0
A.axes2D.yAxis.tickMarks.minorSpacing = 0.5
A.axes2D.yAxis.tickMarks.majorSpacing = 2.0
A.axes2D.yAxis.grid = 0

A.axes2D.xAxis.title.font.font = A.axes2D.yAxis.title.font.Times # Arial, Courier, Times
A.axes2D.xAxis.title.font.scale = 1.5
A.axes2D.xAxis.title.font.useForegroundColor = 1
A.axes2D.xAxis.title.font.color = (0, 0, 0, 255)
A.axes2D.xAxis.title.font.bold = 0
A.axes2D.xAxis.title.font.italic = 0
A.axes2D.xAxis.title.userTitle = 1
A.axes2D.xAxis.title.userUnits = 1
A.axes2D.xAxis.title.title = "Z"
A.axes2D.xAxis.title.units = "pc"

A.axes2D.yAxis.title.font.font = A.axes2D.yAxis.title.font.Times # Arial, Courier, Times
A.axes2D.yAxis.title.font.scale = 1.5
A.axes2D.yAxis.title.font.useForegroundColor = 1
A.axes2D.yAxis.title.font.color = (0, 0, 0, 255)
A.axes2D.yAxis.title.font.bold = 0
A.axes2D.yAxis.title.font.italic = 0
A.axes2D.yAxis.title.userTitle = 1
A.axes2D.yAxis.title.userUnits = 1
A.axes2D.yAxis.title.title = "R"
A.axes2D.yAxis.title.units = "pc"

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
txtlabel.position = 0.55, 0.55
txtlabel.height = 0.04
txtlabel.textColor = (0, 0, 0, 255)
txtlabel.useForegroundForTextColor = 1
txtlabel.text = "Time =  %1.3f" %tkyr
txtlabel.fontFamily = txtlabel.Times  # Arial, Courier, Times
txtlabel.fontBold = 0
txtlabel.fontItalic = 0
txtlabel.fontShadow = 0

txtlbl2=CreateAnnotationObject("Text2D")
txtlbl2.visible = 1
txtlbl2.active = 1
txtlbl2.position = 0.6, 0.65
txtlbl2.height = 0.05
txtlbl2.textColor = (0, 0, 0, 255)
txtlbl2.useForegroundForTextColor = 1
txtlbl2.text = "$outfile"
txtlbl2.fontFamily = txtlbl2.Times  # Arial, Courier, Times
txtlbl2.fontBold = 0
txtlbl2.fontItalic = 0
txtlbl2.fontShadow = 0


DrawPlots()
va=GetView2D()
#va.windowCoords = (0.15,0.85, 0.1, 1.0)
va.viewportCoords = (0.15,0.85, 0.1, 1.0)
SetView2D(va)


###################
## Make a Legend ##
###################
plotName = GetPlotList().GetPlots(0).plotName 
legend = GetAnnotationObject(plotName)
legend.managePosition = 0
legend.orientation = legend.VerticalRight
legend.xScale = 1.0
legend.yScale = 1.5
legend.position = (0.85,0.49)
legend.drawTitle = 0
legend.drawLabels = 1
legend.fontFamily=legend.Times
legend.fontHeight=0.02
legend.position = (0.85,0.7)

### See if we can scale the legend.
##legend.xScale = 1.0
##legend.yScale = 0.5
### the bounding box (draw=0 no, =1 yes)
##legend.drawBoundingBox = 0
##legend.boundingBoxColor = (200,200,200,230)
### number format
#legend.numberFormat = "%1.2e"
### the font.
#legend.fontFamily = legend.Times
##legend.fontBold = 1
##legend.fontItalic = 0
#legend.fontHeight = 0.04
### turning off the title.
#legend.drawTitle = 1
### Make it horizontal
##legend.orientation = legend.VerticalLeft
### moving the legend
#legend.managePosition = 0
#legend.position = (0.75,1.05)
##legend.position = (0.75,1.575)
### text color
###InvertBackgroundColor()
###legend.useForegroundForTextColor = 0
###legend.textColor = (255, 0, 0, 255)
### turning off the labels.
#legend.drawLabels = 1
#legend.drawMinMax = 1
#legend.numTicks = 4

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
	tkyr=ttt/3.16e10
	stst = "Time =  %1.3f kyr" %tkyr
	txtlabel.text = stst
	sw.fileName = "%s_%02d" %("${outfile}_Dens", state)
	SetSaveWindowAttributes(sw)
	SaveWindow()
quit()

#SetTimeSliderState(15)
#ttt=GetQueryOutputValue()
#tkyr=ttt/3.16e10
#stst = "Time =  %1.3f kyr" %tkyr
#txtlabel.text = stst
EOF

$visit_cmd -cli -nowin -s plot_BWaxi2d_results.py

#exit

