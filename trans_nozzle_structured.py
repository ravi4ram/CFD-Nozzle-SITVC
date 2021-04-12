#
# script to design and generate mesh for conical nozzle with SITVC port 
# in .su2 format for CFD analysis.   
#
# Inputs
#     conical nozzle : area_ratio, throat_radius
# Outputs
#     conical nozzle : conical_nozzle_cgrid.su2
# 
# Author
#     Ravi
#

import gmsh
import math
import warnings
import numpy as np
from bisect import bisect_left

#  area_ratio, throat_radius
def mesh_conical_nozzle(aratio, Rt):
	# half cone angle
	half_cone_angle = 15
	alpha           = math.radians(half_cone_angle)
	# exit radius
	Re = Rt * math.sqrt(aratio)	
	# entrant functions
	entrant_angle  	= -135 # -135 -> yi = 601.64
	ea_radian 		= math.radians(entrant_angle)	
	R1 = 2.0 * Rt
	ea_start 		= ea_radian
	ea_end 			= -math.pi/2	
	x1 = ( R1 * math.cos(ea_start) )
	y1 = ( R1 * math.sin(ea_start) + 3.0 * Rt )
	x2 = ( R1 * math.cos(ea_end) )
	y2 = ( R1 * math.sin(ea_end) + 3.0 * Rt )	
	# exit section
	R1 = 1.5 * Rt	
	L1 = R1 * math.sin(alpha)
	Rn = Rt + R1 * ( 1 - math.cos(alpha) )
	Ln = (Re - Rn) / math.tan(alpha)
	L = L1 + Ln	
	# sitvc location (where A/At = 2.5)	
	Rs = math.sqrt(2.5) * Rt
	Ls = L1 + ( (Rs - Rn) / (Re - Rn) ) * Ln
	# create a port
	Rs1 = math.sqrt(2.485) * Rt
	Ls1 = L1 + ( (Rs1 - Rn) / (Re - Rn) ) * Ln
	Rs2 = math.sqrt(2.515) * Rt
	Ls2 = L1 + ( (Rs2 - Rn) / (Re - Rn) ) * Ln

	# gmsh packages
	model = gmsh.model
	geo = model.geo	
	mesh = model.mesh
	field = mesh.field	
	# init gmsh
	gmsh.initialize()
	model.add("conical_nozzle")
	# mesh size
	h = 100
	# z co-ordinate
	z = 0
	# inlet origin	
	geo.addPoint( x1, 0, z, h * 0.01, 301 )
	# entry arc points	
	geo.addPoint( x1, y1, z, h * 0.01, 1 )
	geo.addPoint( x2, y2, z, h * 0.01, 2 )
	# down
	geo.addPoint( x1, -y1, z, h * 0.01, 5 )
	geo.addPoint( x2, -y2, z, h * 0.01, 6 )	
	# inlet arc origin	
	geo.addPoint( 0, 3.0 * Rt, z, h * 0.01, 302 )
	geo.addPoint( 0, -3.0 * Rt, z, h * 0.01, 312 )
	# inlet arc
	geo.addCircleArc( 1, 302, 2, 61 )
	geo.addCircleArc( 5, 312, 6, 63 )
	# throat origin point	
	geo.addPoint( 0, 2.5 * Rt, z, h * 0.01, 303 )
	geo.addPoint( 0, -2.5 * Rt, z, h * 0.01, 313 )
	# throat curve (connect inlet last point with bell first point)
	geo.addPoint( L1, Rn, z, h * 0.01, 3 )
	geo.addPoint( L1, -Rn, z, h * 0.01, 7 )
	geo.addCircleArc( 2, 303, 3, 62 )
	geo.addCircleArc( 6, 313, 7, 64 )
	# symmetry point underneath throat origin
	geo.addPoint( 0, 0, z, h * 0.01, 304 )
	geo.addPoint( L1, 0, z, h * 0.01, 305 )
	# sitvc port
	geo.addPoint( Ls1, Rs1, z, h * 0.01, 307 )
	geo.addPoint( Ls2, Rs2, z, h * 0.01, 308 )
	# end of cone
	geo.addPoint( L, Re, z, h * 0.01, 4 )
	geo.addPoint( L, -Re, z, h * 0.01, 8 )
	# cone lines
	geo.addLine( 3, 307, 50 )
	geo.addLine( 307, 308, 51 )
	geo.addLine( 308, 4, 52 )	
	geo.addLine( 7, 8, 57 )
	# symmetry point underneath throat exit	
	geo.addPoint( L, 0, z, h * 0.01, 306 )
	# top loop lines
	geo.addLine( 301, 1, 41 )
	geo.addLine( 2, 304, 42 )
	geo.addLine( 304, 301, 43 )
	geo.addLine( 3, 305, 44 )
	geo.addLine( 305, 304, 45 )
	geo.addLine( 4, 306, 46 )	
	geo.addLine( 306, 305, 47 )
	# bottom loop lines
	geo.addLine( 301, 5, 53 )
	geo.addLine( 6, 304, 54 )
	geo.addLine( 7, 305, 55 )
	geo.addLine( 8, 306, 56 )
	# top Curveloop and Surface
	geo.addCurveLoop( [-41, -61, -42, -43], 601 )
	geo.addPlaneSurface( [601], 801 )
	geo.addCurveLoop( [42, -62, -44, -45], 602 )
	geo.addPlaneSurface( [602], 802 )	
	geo.addCurveLoop( [44, -50, -51, -52, -46, -47], 603 )
	geo.addPlaneSurface( [603], 803 )	
	# bottom Curveloop and Surface
	geo.addCurveLoop( [53, 43, 54, 63], 604 )
	geo.addPlaneSurface( [604], 804 )
	geo.addCurveLoop( [-54, 45, 55, 64], 605 )
	geo.addPlaneSurface( [605], 805 )	
	geo.addCurveLoop( [-55, 47, 56, 57], 606 )
	geo.addPlaneSurface( [606], 806 )	

	# synchronize
	geo.synchronize()

	# transfinite
	numCellsX = 25 
	gradingX = 1.05; gradingY = 1.0;

	# entrant top
	numCellsY = 5
	mesh.setTransfiniteCurve(41, numCellsX, meshType="Progression", coef=-gradingX)
	mesh.setTransfiniteCurve(42, numCellsX, meshType="Progression", coef=1/gradingX)
	
	mesh.setTransfiniteCurve(61, numCellsY, meshType="Progression", coef=gradingY)
	mesh.setTransfiniteCurve(43, numCellsY, meshType="Progression", coef=gradingY)		
	mesh.setTransfiniteSurface( 801, cornerTags=[301, 304, 2, 1] ) 
	# entrant bottom
	mesh.setTransfiniteCurve(53, numCellsX, meshType="Progression", coef=-gradingX)
	mesh.setTransfiniteCurve(54, numCellsX, meshType="Progression", coef=1/gradingX)
	
	mesh.setTransfiniteCurve(63, numCellsY, meshType="Progression", coef=gradingY)
	mesh.setTransfiniteCurve(43, numCellsY, meshType="Progression", coef=gradingY)		
	mesh.setTransfiniteSurface( 804, cornerTags=[301, 304, 6, 5] ) 
	
	# throat top
	numCellsY = 3
	mesh.setTransfiniteCurve(42, numCellsX, meshType="Progression", coef=gradingX)
	mesh.setTransfiniteCurve(44, numCellsX, meshType="Progression", coef=1/gradingX)
	
	mesh.setTransfiniteCurve(62, numCellsY, meshType="Progression", coef=gradingY)
	mesh.setTransfiniteCurve(45, numCellsY, meshType="Progression", coef=gradingY)		
	mesh.setTransfiniteSurface( 802, cornerTags=[304, 305, 3, 2] ) 	
	# throat bottom
	mesh.setTransfiniteCurve(54, numCellsX, meshType="Progression", coef=gradingX)
	mesh.setTransfiniteCurve(55, numCellsX, meshType="Progression", coef=1/gradingX)
	
	mesh.setTransfiniteCurve(64, numCellsY, meshType="Progression", coef=gradingY)
	mesh.setTransfiniteCurve(45, numCellsY, meshType="Progression", coef=gradingY)		
	mesh.setTransfiniteSurface( 805, cornerTags=[304, 305, 7, 6] )
	
	# cone top	
	numCellsY = 20 
	mesh.setTransfiniteCurve(44, numCellsX, meshType="Progression", coef=gradingX)
	mesh.setTransfiniteCurve(46, numCellsX, meshType="Progression", coef=1/-gradingX)
	
	numCellsY_1 = int(numCellsY/2); numCellsY_2 = 2; numCellsY_3 = int(numCellsY/2);
	numCellsY_4 = (numCellsY_1 + numCellsY_2 + numCellsY_3 - 2)		
	mesh.setTransfiniteCurve(50, numCellsY_1, meshType="Progression", coef=gradingY)
	mesh.setTransfiniteCurve(51, numCellsY_2, meshType="Progression", coef=gradingY)
	mesh.setTransfiniteCurve(52, numCellsY_3, meshType="Progression", coef=gradingY)
	mesh.setTransfiniteCurve(47, numCellsY_4, meshType="Progression", coef=gradingY)		
	mesh.setTransfiniteSurface( 803, cornerTags=[305, 306, 4, 3] ) 	
	# cone bottom
	mesh.setTransfiniteCurve(55, numCellsX, meshType="Progression", coef=gradingX)
	mesh.setTransfiniteCurve(56, numCellsX, meshType="Progression", coef=1/-gradingX)
	
	mesh.setTransfiniteCurve(57, numCellsY_4, meshType="Progression", coef=gradingY)
	mesh.setTransfiniteCurve(47, numCellsY_4, meshType="Progression", coef=gradingY)			
	mesh.setTransfiniteSurface( 806, cornerTags=[305, 306, 8, 7] ) 	

	# boundaries
	dim = 1
	gmsh.model.addPhysicalGroup(dim, [61, 62, 63, 64, 50, 51, 52, 57], 1001)
	gmsh.model.setPhysicalName(dim, 1001, "WALL")
	gmsh.model.addPhysicalGroup(dim, [41, 53], 1002)
	gmsh.model.setPhysicalName(dim, 1002, "INFLOW1")
	gmsh.model.addPhysicalGroup(dim, [51], 1003)
	gmsh.model.setPhysicalName(dim, 1003, "INFLOW2")
	gmsh.model.addPhysicalGroup(dim, [46, 56], 1004)
	gmsh.model.setPhysicalName(dim, 1004, "OUTFLOW")
	dim = 2
	gmsh.model.addPhysicalGroup(dim, [801, 802, 803, 804, 805, 806], 2001)
	gmsh.model.setPhysicalName(dim, 2001, "Plane surface")
	
	# generate mesh
	gmsh.option.setNumber("Mesh.RecombineAll", 1)
	gmsh.option.setNumber("General.Terminal", 1)
	gmsh.option.setNumber("Mesh.Smoothing", 100)
	gmsh.option.setNumber("Mesh.Algorithm", 5) # delquad
	mesh.generate(2)
	mesh.refine() # Refine the mesh of the current model by uniformly splitting the elements.
	gmsh.write('conical_nozzle_cgrid.su2')
	gmsh.fltk.run()
	gmsh.finalize()	
	# return		
	return
	
		
# __main method__
if __name__=="__main__":

	# PS1 - S139 ("Development of nozzle for PSLV Booster.pdf)
	throat_radius = 836/2.0 
	exit_radius = 2377/2.0
	aratio = math.pow(exit_radius,2) / math.pow(throat_radius,2) # Ae / At
	
	# design and mesh conical_nozzle_contour
	mesh_conical_nozzle(aratio, throat_radius)



