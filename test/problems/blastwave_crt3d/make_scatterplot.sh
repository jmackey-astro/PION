#!/bin/bash
#
# 2016.08.25 JM: Make a scatterplot of density vs. radius.

FDIR=${1}
FBASE=${2}
SAVEDIR=${3}
SAVEFILE=${4}
SIM_TITLE=${5}


STEP=${6}   # only plot every STEPth snapshot
#STEP=1

mkdir -p $SAVEDIR
VISITFILE=visit_${FBASE}_NVB.py

case $HOSTNAME in
  rvs[0-9])
    echo "Running on SuperMUC visualization cluster"
    source /etc/profile.d/modules.sh
    module load visit
    VISITCMD="rvglrun visit"
    NCORE=32
  ;;
esac

if [ "${HOSTNAME}" = 'w26' ]; then
  echo "Running on w26@koeln"
  VISITCMD="/home/jmac/software/bin/bin/visit"
  NCORE=1

elif [ "${HOSTNAME}" = 'aibn129' ]; then
  echo "Running on aibn129@aifa"
  VISITCMD="/vol/aibn129/data1/jmackey/extra_libraries/visit_bin/bin/visit"
  NCORE=4

elif [ "${HOSTNAME}" = 'nova' ]; then
  echo "Running on nova@dunsink"
  VISITCMD="vglrun /usr/local/bin/visit"
  #VISITCMD="/usr/local/bin/visit"
  NCORE=1
fi

rm ${SAVEDIR}/${SAVEFILE}*.png


cat <<EOF > ${VISITFILE}
DefineVectorExpression("pos3d",    "{coord(MultiMesh)[0],coord(MultiMesh)[1],coord(MultiMesh)[2]}")
DefineScalarExpression("Rad3D_pc",     "magnitude(pos3d)*3.24e-19")
DefineScalarExpression("SiloLogDensity",     "log10(Density)")
DefineScalarExpression("D24",     "Density*cell_constant(MultiMesh,1.0e24)")


args = ("-np", "$NCORE")
OpenComputeEngine("localhost", args)
SetWindowArea(0,0,800,800)

DeleteAllPlots()

# OpenDatabase("localhost:files",timestep)
OpenDatabase("localhost:${FDIR}/${FBASE}.*.silo database",0)


#################################
### Scatter plot of density ##
#################################
AddPlot("Scatter", "Density", 1, 0)
ScatterAtts = ScatterAttributes()
ScatterAtts.var1 = "Rad3D_pc"
ScatterAtts.var1Role = ScatterAtts.Coordinate0  # Coordinate0, Coordinate1, Coordinate2, Color, None
ScatterAtts.var1MinFlag = 0
ScatterAtts.var1MaxFlag = 0
ScatterAtts.var1Min = 0
ScatterAtts.var1Max = 1
ScatterAtts.var1Scaling = ScatterAtts.Linear  # Linear, Log, Skew
ScatterAtts.var1SkewFactor = 1
ScatterAtts.var2 = "Density"
ScatterAtts.var2Role = ScatterAtts.Coordinate1  # Coordinate0, Coordinate1, Coordinate2, Color, None
ScatterAtts.var2MinFlag = 0
ScatterAtts.var2MaxFlag = 0
ScatterAtts.var2Min = 0
ScatterAtts.var2Max = 1
ScatterAtts.var2Scaling = ScatterAtts.Linear  # Linear, Log, Skew
ScatterAtts.var2SkewFactor = 1
ScatterAtts.var3Role = ScatterAtts.None  # Coordinate0, Coordinate1, Coordinate2, Color, None
ScatterAtts.var3 = "default"
ScatterAtts.var3MinFlag = 0
ScatterAtts.var3MaxFlag = 0
ScatterAtts.var3Min = 0
ScatterAtts.var3Max = 1
ScatterAtts.var3Scaling = ScatterAtts.Linear  # Linear, Log, Skew
ScatterAtts.var3SkewFactor = 1
ScatterAtts.var4Role = ScatterAtts.None  # Coordinate0, Coordinate1, Coordinate2, Color, None
ScatterAtts.var4 = "default"
ScatterAtts.var4MinFlag = 0
ScatterAtts.var4MaxFlag = 0
ScatterAtts.var4Min = 0
ScatterAtts.var4Max = 1
ScatterAtts.var4Scaling = ScatterAtts.Linear  # Linear, Log, Skew
ScatterAtts.var4SkewFactor = 1
ScatterAtts.pointSize = 0.05
ScatterAtts.pointSizePixels = 4
ScatterAtts.pointType = ScatterAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
ScatterAtts.scaleCube = 0
ScatterAtts.colorType = ScatterAtts.ColorByForegroundColor  # ColorByForegroundColor, ColorBySingleColor, ColorByColorTable
ScatterAtts.singleColor = (255, 0, 0, 255)
ScatterAtts.colorTableName = "Default"
ScatterAtts.invertColorTable = 0
ScatterAtts.legendFlag = 0
SetPlotOptions(ScatterAtts)


DrawPlots()



#####################################################################
#####################################################################

A = AnnotationAttributes()
A.axes2D.visible = 1
A.axes2D.autoSetTicks = 1
A.axes2D.autoSetScaling = 1
A.axes2D.lineWidth = 1
A.axes2D.tickLocation = A.axes2D.Both  # Inside, Outside, Both
A.axes2D.tickAxes = A.axes2D.BottomLeft  # Off, Bottom, Left, BottomLeft, All
A.axes2D.xAxis.title.visible = 1
A.axes2D.xAxis.title.font.font = A.axes2D.xAxis.title.font.Arial  # Arial, Courier, Times
A.axes2D.xAxis.title.font.scale = 2
A.axes2D.xAxis.title.font.useForegroundColor = 1
A.axes2D.xAxis.title.font.color = (0, 0, 0, 255)
A.axes2D.xAxis.title.font.bold = 0
A.axes2D.xAxis.title.font.italic = 0
A.axes2D.xAxis.title.userTitle = 1
A.axes2D.xAxis.title.userUnits = 1
A.axes2D.xAxis.title.title = "r"
A.axes2D.xAxis.title.units = "pc"
A.axes2D.xAxis.label.visible = 1
A.axes2D.xAxis.label.font.font = A.axes2D.xAxis.label.font.Arial  # Arial, Courier, Times
A.axes2D.xAxis.label.font.scale = 2
A.axes2D.xAxis.label.font.useForegroundColor = 1
A.axes2D.xAxis.label.font.color = (0, 0, 0, 255)
A.axes2D.xAxis.label.font.bold = 1
A.axes2D.xAxis.label.font.italic = 0
A.axes2D.xAxis.label.scaling = 0
A.axes2D.xAxis.tickMarks.visible = 1
A.axes2D.xAxis.tickMarks.majorMinimum = -200
A.axes2D.xAxis.tickMarks.majorMaximum = 200
A.axes2D.xAxis.tickMarks.minorSpacing = 1.0
A.axes2D.xAxis.tickMarks.majorSpacing = 10.0
A.axes2D.xAxis.grid = 0
A.axes2D.yAxis.title.visible = 1
A.axes2D.yAxis.title.font.font = A.axes2D.yAxis.title.font.Arial  # Arial, Courier, Times
A.axes2D.yAxis.title.font.scale = 2
A.axes2D.yAxis.title.font.useForegroundColor = 1
A.axes2D.yAxis.title.font.color = (0, 0, 0, 255)
A.axes2D.yAxis.title.font.bold = 0
A.axes2D.yAxis.title.font.italic = 0
A.axes2D.yAxis.title.userTitle = 1
A.axes2D.yAxis.title.userUnits = 1
A.axes2D.yAxis.title.title = "Density"
A.axes2D.yAxis.title.units = "g/cm^3"
A.axes2D.yAxis.label.visible = 1
A.axes2D.yAxis.label.font.font = A.axes2D.yAxis.label.font.Arial  # Arial, Courier, Times
A.axes2D.yAxis.label.font.scale = 2
A.axes2D.yAxis.label.font.useForegroundColor = 1
A.axes2D.yAxis.label.font.color = (0, 0, 0, 255)
A.axes2D.yAxis.label.font.bold = 1
A.axes2D.yAxis.label.font.italic = 0
A.axes2D.yAxis.label.scaling = 0
A.axes2D.yAxis.tickMarks.visible = 1
A.axes2D.yAxis.tickMarks.majorMinimum = -200
A.axes2D.yAxis.tickMarks.majorMaximum = 200
A.axes2D.yAxis.tickMarks.minorSpacing = 1.0
A.axes2D.yAxis.tickMarks.majorSpacing = 10.0
A.axes2D.yAxis.grid = 0
A.userInfoFlag = 0
A.userInfoFont.font = A.userInfoFont.Arial  # Arial, Courier, Times
A.userInfoFont.scale = 1
A.userInfoFont.useForegroundColor = 1
A.userInfoFont.color = (0, 0, 0, 255)
A.userInfoFont.bold = 0
A.userInfoFont.italic = 0
A.databaseInfoFlag = 0
A.databaseInfoFont.font = A.databaseInfoFont.Arial  # Arial, Courier, Times
A.databaseInfoFont.scale = 1
A.databaseInfoFont.useForegroundColor = 1
A.databaseInfoFont.color = (0, 0, 0, 255)
A.databaseInfoFont.bold = 0
A.databaseInfoFont.italic = 0
A.databaseInfoExpansionMode = A.File  # File, Directory, Full, Smart, SmartDirectory
A.databaseInfoTimeScale = 3.1685678e-14
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
A.axesArray.axes.title.font.font = A.axesArray.axes.title.font.Arial  # Arial, Courier, Times
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
A.axesArray.axes.label.font.font = A.axesArray.axes.label.font.Arial  # Arial, Courier, Times
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

#tkyr=0
text=CreateAnnotationObject("Text2D")
text.visible = 1
text.active = 1
text.position = 0.08, 0.965
text.height = 0.025
text.textColor = (0, 0, 0, 255)
text.useForegroundForTextColor = 1
text.text = "${SIM_TITLE} : Density as function of Radius"
text.fontFamily = text.Arial  # Arial, Courier, Times
text.fontBold = 0
text.fontItalic = 0
text.fontShadow = 0


#####################################################################
##########  GET THE VIEW AND IMAGE RESOLUTION THAT WE WANT ##########
#####################################################################
DrawPlots()
ResizeWindow(1,900,800)

sw = SaveWindowAttributes()
sw.fileName = "test"
sw.width = 1200
sw.height= 1000
#sw.width = 2816
#sw.height= 3072
sw.outputToCurrentDirectory = 0
sw.outputDirectory = "${SAVEDIR}"
sw.family = 0
sw.format = sw.PNG # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK
sw.resConstraint = sw.NoConstraint

#View2DAtts = View2DAttributes()
#View2DAtts.windowCoords = ( -45, 45, -45, 45)
#View2DAtts.viewportCoords = (0.1, 0.9, 0.1, 0.9)
#View2DAtts.fullFrameActivationMode = View2DAtts.Off  # On, Off, Auto
#View2DAtts.fullFrameAutoThreshold = 100
#View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
#View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
#SetView2D(View2DAtts)

#####################################################################



#####################################################################
# Add a time slider in the upper left corner
slider = CreateAnnotationObject("TimeSlider")
slider.height = 0.08
slider.width = 0.4
slider.text = "\$time Myr"
slider.position = (0.6, 0.2)
# Make the slider transparent, so that only the time shows.
slider.startColor = (0, 255, 255, 0)
slider.endColor = (255, 255, 255, 0)
slider.textColor = (0, 0, 0, 255)
slider.useForegroundForTextColor = 1
slider.timeFormatString = "%.3g"
#####################################################################



#####################################################################
##########  LOOP OVER ALL TIMESTEPS AND SAVE IMAGES         #########
#####################################################################

for state in range(0,TimeSliderGetNStates(),${STEP}):
	SetTimeSliderState(state)
	sw.fileName = "%s_%03d" %("${SAVEFILE}", state)
	SetSaveWindowAttributes(sw)
	SaveWindow()
        #CloseComputeEngine()
quit()
EOF
#
${VISITCMD} -cli -nowin -s ${VISITFILE}
#${VISITCMD} -cli -s ${VISITFILE}

exit


rm ${SAVEDIR}/${SAVEFILE}_lores.mp4
ffmpeg  -r 4.0 -f image2  -i ${SAVEDIR}/${SAVEFILE}_%03d.png \
 -q:v 0 -s 768x400 -pix_fmt yuv420p -vcodec h264 ${SAVEDIR}/${SAVEFILE}_lores.mp4

rm ${SAVEDIR}/${SAVEFILE}_ani.mp4
ffmpeg  -r 4.0 -f image2  -i ${SAVEDIR}/${SAVEFILE}_%03d.png \
 -q:v 0 -pix_fmt yuv420p -vcodec h264 ${SAVEDIR}/${SAVEFILE}_ani.mp4


exit





