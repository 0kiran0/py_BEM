#***********************************************************************************
# Iterative code for calculating axial induction factor (a)
#  and angular induction factor (a') using Blade Element Momentum Theory.
# ---- Kiran Siddappaji
# ---- University of Cincinnati
#-----------------------------------------------------------------------------------
# Inputs: Fluid properties,Flow Velocity,TSR,Rhub,Rtip, nbrows, nblades,nairfoils
#       : Airfoil properties (Re,type,AoA,Cl(AoA),Cd(AoA)), spanwise chord (option)
#
# Outputs: a, a', spanwise properties, Total Power and Thrust, Efficiency
#        : 3dbgbinput file, radial plots of the spanwise properties
#***********************************************************************************
import sys,os, errno
import math, timeit
import numpy as np
import pylab as py
from matplotlib.pyplot import *
from functions import * 
from functions2 import * 
from scipy.optimize import fsolve
from scipy.optimize import brentq
import mpl_toolkits.mplot3d.axes3d as p3
import subroutines
from spanwiseplots import *


#timing the code
start = timeit.default_timer()
#
workdir = os.getcwd()
print
print "Working Directory: ", workdir
print
## Directory for plots
plotlogic = sys.argv[2]
if plotlogic == 'T':
	plot = True
elif plotlogic == 'F':
	plot = False 
####
##Propeller or Turbine
if device == "Propeller" or device == "propeller":
	propeller = True
elif device == "Turbine" or device == "turbine":
	propeller = False
print
print "The device is a ", case,"",device,"." 
if nbrow == 2:
	print "The device is a COUNTER ROTATING ", device, "."
print
#---------------------------------------------------------
# Setting co/counter rotation system as false to begin with
counter_rotation = False
co_rotation = False
# Setting multirow as flase to begin with
multirow = False
#
# Setting currentrow as 1 to start with
currentrow = 1
# Setting the axial gap between rotors to 2 for single row configurations
if nbrow == 1:
	gap = 2.0#float(sys.argv[3])
elif nbrow == 2:
	gap = float(sys.argv[3])
#---------------------------------------------------------
# Radial stations
r = np.linspace(Rhub, Rtip, num=nsect)
r1 = r
#print r
#---------------------------------------------------------
# OMEGA
print "Vztip: ",Vz1[nsect-1]
if inputRPM: # RPM is input and need to calculate TSR
	print "input RPM  : ", RPMin
	Lambda = calcTSR(RPMin,Vz,Rtip)
	OMEGA  = Lambda*Vz1[nsect-1]/Rtip 
	RPM = RPMin
else: # TSR is input and need to calculate RPM
	OMEGA = Lambda*Vz1[nsect-1]/Rtip
	RPM   = OMEGA*(60./(2.*math.pi))
#
print "OMEGA: ", OMEGA, "rad/s"
print
if not inputRPM:
	print "RPM  : ", RPM
print
print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print
#----------------------------------------------------------
#----------------------------------------------------------
# Co/Counter rotation system activation
if nbrow == 2:
	if (RPM > 0.) and (RPM2read > 0.):
		co_rotation = True
	elif (RPM > 0.) and (RPM2read < 0.):
		counter_rotation = True
	elif (RPM2read == 0.):
		counter_rotation = True
		
	print " co_rotation: ", co_rotation
	print " counter_rotation: ", counter_rotation 

#----------------------------------------------------------
#----------------------------------------------------------
# Vzfront = np.zeros(nsect)
for i in range(nsect):
	span[i] = (r[i]-Rhub)/(Rtip-Rhub)          
			 
	# Lambda at each radial station
	lambda_r[i] = OMEGA*r[i]/Vz1[i]
	# lambda_r[i] = Lambda*(r[i]/Rtip)
	# print lambda_r[i]

	# Local Wheel speed [m/s]
	U[i] = lambda_r[i]*Vz1[i] 
	# Vzfront[i] = 15.0
# end of for loop
#------------------------------------------------------------------
#------------------------------------------------------------------
# Chord definition
# Input chord at each section
if not chordspline:
	if definedchord:
		print chord1defn
	chrdmplier = [1.0]*nsect
# Spline chord through control points
if chordspline:
	#---------------------------------------------------------
	# Control points for spanwise distribution of chord multiplier 
	#---------------------------------------------------------
	print
	print " Control points span vs chord mutliplier"
	print "   spanCP     chrdMPCP"
	for i in range(ncp):
		print "  ",'%1.6f'%spancp[i]," ",'%1.6f'%chrd_multipcp[i]
	print
	#---------------------------------------------------------
	#---------------------------------------------------------
	chrdmplier = subroutines.cubicspline(span,nsect,chrd_multipcp,spancp,ncp)
	print " Span vs Chord multiplier spline : "
	print "    span       chrdMP"
	for i in range(nsect):
		print "  ",'%1.6f'%span[i]," ",'%1.6f'%chrdmplier[i]
#---------------------------------------------------------
#---------------------------------------------------------
if spanwise:
	print "--------------------------------------------------"
	#---------------------------------------------------------
	# Control points for spanwise distribution of AoA
	#---------------------------------------------------------
	print
	print " Control points span vs AoA"
	print "   spanCP     alpha1CP"
	for i in range(alpha1_ncp):
		print "  ",'%1.6f'%alpha1span_cp[i]," ",'%1.6f'%alpha1_cp[i]
	print
	#---------------------------------------------------------
	#---------------------------------------------------------
	alpha1span = subroutines.cubicspline(span,nsect,alpha1_cp,alpha1span_cp,alpha1_ncp)
	print " Span vs AoA spline: "
	print "    span       alpha1"
	for i in range(nsect):
		print "  ",'%1.6f'%span[i]," ",'%1.2f'%alpha1span[i]

	if foilname == 'e857m':
		airfoilname = 'e857'  
	else:
		airfoilname = foilname
	#Clspan, Cdspan, L_over_D = xfoilrun(alpha1span,nsect,REYN,airfoilname) 
#	Clspan,Cdspan =	Cl_Cd_rangeREYN(foilname,alpha,propeller,Re_lookup)
#	print Clspan
#	print Cdspan
	# print
	# print Clspan, Cdspan, L_over_D       
	# print
# stop
#---------------------------------------------------------
#---------------------------------------------------------
# Defined tolerance
deftol     = 1.0e-6
#
print
print "Iterative procedure to calculate a and a' for all 2D airfoils:"
print
print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print
# Bisection method
# upper and lower limits for phi
if propeller:
	phi_lower = 1.*math.pi/180.
	phi_upper = 89.*math.pi/180.
else: #WT/HKT
	phi_lower = 1.*math.pi/180.
	phi_upper = 89.*math.pi/180.   
tol = 1e-6
k = 0
while k < 2: 
	for i in range(nsect):
		print "------------------------"
		print "Airfoil   :", i+1
		print "------------------------"
		print
		## Brent's method of convergence to find a and a'
		#print "Using Brentq's method for calculating a and a'..."
		#
		L_over_D[i] = Clspan[i]/float(Cdspan[i])
		#
		if varyingREYN and spanwise and k == 1:
			print " variable REYN and AoA spanwise" 
			alpha = alpha1span[i]
			Cl      = Clspan[i]
			Cd     = Cdspan[i]   
		elif varyingREYN and k == 1:
			print " variable REYN spanwise and constant AoA"
			Cl     = Clspan[i]
			Cd    = Cdspan[i]   
			alpha1span[i] = alpha
			print 'constant AoA :', alpha1span[i]
		elif spanwise:
			print " variable AoA spanwise "
			alpha = alpha1span[i]
			Cl      = Clspan[i]
			Cd     = Cdspan[i]
		##############
		if definedchord:
			chordr[i] = chord1defn[i]
			internalchord = False 
		######    
		####
		if internalchord and not definedchord:
			chordr[i]  = calcChord(r[i],Phirad[i],nblade,Cl,lambda_r[i])# optimal chord (obtained for max Power)
			chordr[i]  = math.fabs(chordr[i])
		####
		chordr[i]   = chordr[i]*chrdmplier[i]
		# print "chordr[i] : ", chordr[i]        
		#print "Cl,Cd:",Cl,Cd
		# print "Cl,Cd before the routine:",Cl,Cd	
		if not propeller:
		# Brentq method
			# Phirad[i],axlIFnew[i],angIFnew[i],floss[i],Ca[i],Cr[i],sigmar[i],chordr[i] = calc_a_aprime(currentrow,co_rotation,counter_rotation,propeller,file_length,
																							   # chordr[i],chrdmplier[i],lambda_r[i],nbrow,
																							   # r[i],Rhub,Rtip,nblade,Cl,Cd,gap)
		# ###
		# Fixed Point Iteration method###
			Phirad[i],axlIFnew[i],angIFnew[i],floss[i],Ca[i],Cr[i],sigmar[i],chordr[i]	= 	fixedPointIteration(currentrow,co_rotation,counter_rotation,propeller,
																								   file_length,chordr[i],chrdmplier[i],lambda_r[i],nbrow,r[i],
																								   Rhub,Rtip,nblade,Cl,Cd,gap)
		# #
		elif propeller: 
			rcurrent      = r[i]
			chord_current = chordr[i]
			Vz_current    = Vz1[i]
			if nbrow == 1 or gap == 2.0:
				arguments_in  = propeller,nblade,Lambda,Rtip,Vz_current,OMEGA,chord_current,rcurrent,Cl,Cd
				# Phirad[i] = bisection(phi_lower,phi_upper,arguments_in,tol) 
				Phirad[i] = fsolve(aPhi,0.1,[propeller,nblade,Lambda,Rtip,Vz_current,OMEGA,chord_current,rcurrent,Cl,Cd])	
			elif nbrow == 2 and gap < 2.0: # multirow
				print " Using axial gap value of ", gap
				arguments_in  = nblade,Lambda,Rtip,Vz_current,OMEGA,chord_current,rcurrent,Cl,Cd,gap        
				Phirad[i] = root_bisection(phi_lower,phi_upper,arguments_in,tol)         
			
			if Phirad[i] == 0.:
				Phirad[i] = 0.000001  
			if (i == nsect-1): # last section extrapolation
				Phirad[i] = Phirad[i-1] - abs(Phirad[i-1] - Phirad[i-2])           
			
			print "Phi[i] :",(Phirad[i]*180./math.pi)
				   
			Ca[i],Cr[i]  = Ca_Cr(Cl,Cd,Phirad[i],propeller)       
			floss[i]     = total_loss(nblade,lambda_r[i],Phirad[i],r[i],Rhub,Rtip)	
			sigmar[i]    = nblade*chordr[i]/float(2*math.pi*r[i])
			Term1        = 1./(lambda_r[i]*math.tan(Phirad[i]))
			Nr           = sigmar[i]*Vz1[i]*Cr[i]
			Dr           = 4.*floss[i]*OMEGA*(r[i])*(math.sin(Phirad[i]))**2
			if Dr == 0. :
				Dr = 0.000001
			Term2        = Nr/Dr
			if propeller:
				bterm           = (1.0/float(Term1 + Term2))    
				
				angIFnew[i]  = 1.0 - bterm/(lambda_r[i]*math.tan(Phirad[i])) 
				axlIFnew[i]   = bterm-1.0
			elif not propeller: #turbine WT/HKT
				bterm           = (1.0/float(Term1 - Term2))    
				
				angIFnew[i]  = (bterm*Term1) - 1
				axlIFnew[i]   = -bterm+1.0 
			# #####
			# if abs(axlIFnew[i]) > 0.40: # Using correction to calculate a from Glauert and Spera correction
				# axlIFnew[i] = axlIF_Glauert(axlIFnew[i],floss[i],Phirad[i],sigmar[i],Ca[i])        
		# print 'Term1, Term2 in b:', Term1, Term2                
			####            
			# if (i == nsect-1): #extrapolating at the tip section
				# if propeller:
					# axlIFnew[i]  = axlIFnew[i-1] - (axlIFnew[i-2] - axlIFnew[i-1])
					# angIFnew[i]  = angIFnew[i-1] - (angIFnew[i-2] - angIFnew[i-1])   
				# else:
					# axlIFnew[i]  = axlIFnew[i-1] + (abs(axlIFnew[i-2]) - abs(axlIFnew[i-1]))
					# angIFnew[i]  = angIFnew[i-1] + (abs(angIFnew[i-2]) - abs(angIFnew[i-1]))		
				###
		# print 'bterm,(lambda_r[i] x math.tan(Phirad[i])):',bterm,(lambda_r[i]*math.tan(Phirad[i]))        
		print "phi,axlIF,angIF: ",(Phirad[i]*180./math.pi),axlIFnew[i],angIFnew[i]
		#Calculating final velocities at various axial stations
		#=================================================
		# if not propeller:
		  # axlIFnew[i]  = - axlIFnew[i]
		  # angIFnew[i] = - angIFnew[i]		  
		V2[i],V3[i],Vexit[i] = velocitycalc(propeller,Vz1[i],axlIFnew[i])
		#=================================================	
		if (i == nsect-1): # last section extrapolation 
			V2[i]     = V2[i-1] + (V2[i-1] - V2[i-2]) 
			V3[i]     = V3[i-1] + (V3[i-1] - V3[i-2]) 
			Vexit[i]  = Vexit[i-1] + (Vexit[i-1] - Vexit[i-2])             
	#-------------------------------------------------------------------------------
	# end of for loop
	#-------------------------------------------------------------------------------
	#-------------------------------------------------------------------------------            
	for i in range(nsect):            
		#-------------------------------------------------------------------------------
		### Final Geometry parameters
		#-------------------------------------------------------------------------------  
		angular_term = 0.0
		Phi[i],inbeta[i],pitch[i],stggr[i],Wr[i] = geometry_param(propeller,Phirad[i],alpha,Vz1[i],axlIFnew[i],
																  angIFnew[i],lambda_r[i],currentrow,co_rotation,counter_rotation,gap,angular_term)
		# ## For stggr, removing anomalies in the curve
	if (i == nsect-1) and propeller: # last section extrapolation
		stggr[i]  = stggr[i-1] + (stggr[i-1] - stggr[i-2])
		# angIFnew[i] = -1.0
		# axlIFnew[i] = -1.0        
		#-------------------------------------------------------------------------------
	# end of for loop
	#--------------------------------------------------------------------
	# Reynolds_Number Calculation
	Re_lookup = np.zeros(nsect)
	Re_lookup,REYNfront = Reyn_Calc(Wr,chordr,nju,propeller,alpha,foilname,plot,span,currentrow,k)
	k = k + 1
	print "k:",k
	#### Discarding very high REYN cases for optimization. Temporary fix.
	for j in range(len(Re_lookup)):
		if Re_lookup[j] > 6.0:
			print " Too large a Reynolds number [e6]:",Re_lookup[j]
			varyingREYN = False
			#k = 2  
			break
	if spanwise or varyingREYN and k == 1:
		print "k is 1:",k	
		if varyingREYN and not spanwise: # constant AoA
			alpha1span = [float(alpha)]*nsect
		# print "alpha going in is ", alpha1span
		Clspan,Cdspan =	Cl_Cd_rangeREYN(foilname,alpha1span,propeller,Re_lookup)
		print Clspan
		print Cdspan
######################################################################	
### End of inductionfactor and REYN loop
print "Clspan:",Clspan
print
print "Cdspan:",Cdspan
print
print "Ca:",Ca
print
print "Cr:",Cr
# average radius
r_avg = sum(r)/float(len(r))
print
print "Average radius:",r_avg," m."
#---------------------------------------------------
# Radial Distribution of properties
#---------------------------------------------------	
print "dThrustBETdr[i],dAForceBMTdr[i], dTorqueBETdr[i], dTForceBMTdr[i], dPowerBETdr[i],dPowerBMTdr[i] :"
for i in range(nsect):
	# if not propeller:
	  # axlIFnew[i]  = - axlIFnew[i]
	  # angIFnew[i] = - angIFnew[i]
	  
	Vtheta1slipstream = np.zeros(nsect)      

	Vtheta1slipstream[i] = 2.0*angIFnew[i]*OMEGA*r[i]
	#-------------------------------------------------------------------------------
	# Blade Element theory definitions
	#-------------------------------------------------------------------------------
	dmdotdr[i], dThrustBETdr[i], dTorqueBETdr[i], dPowerBETdr[i], Vthetar[i] = radial_distr_BET(multirow,currentrow,co_rotation,counter_rotation,propeller,rho,Vz1[i],
																								sigmar[i],axlIFnew[i],r[i],Phirad[i],Ca[i],Cr[i],OMEGA,nbrow)
	# print "dThrustBETdr[i], dTorqueBETdr[i], dPowerBETdr[i]: ",dThrustBETdr[i], dTorqueBETdr[i], dPowerBETdr[i]
	#---------------------------------
	if (i == nsect-1):
		Vthetar[i] = Vthetar[i-1] + (Vthetar[i-1] - Vthetar[i-2])  
	# Power in the swirl downstream of the rotor
	dPower_swirldr[i] = 0.5*rho*Vz1[i]*(Vthetar[i]**2)*2.*math.pi*r[i] 
	#-------------------------------------------------------------------------------
	# Blade Momentum theory definitions
	#-------------------------------------------------------------------------------  
	value = 0.
	dAForceBMTdr[i], dTForceBMTdr[i], dPowerBMTdr[i] = radial_distr_BMT(multirow,currentrow,co_rotation,counter_rotation,propeller,rho,Vz1[i],
																		axlIFnew[i],angIFnew[i],value,r[i],floss[i],OMEGA,nbrow,gap)
	#############
	if (i == nsect-1): #Forcing the last section for the same thrust
		dAForceBMTdr[i] = dThrustBETdr[i] = 0.0
		dTForceBMTdr[i] = dTorqueBETdr[i] = 0.0
		dPowerBMTdr[i]  = dPowerBETdr[i]  = 0.0      
	#############
	print '%.3f'%dThrustBETdr[i],'%.3f'%dAForceBMTdr[i], '%.3f'%dTorqueBETdr[i], '%.3f'%dTForceBMTdr[i], '%.3f'%dPowerBETdr[i],'%.3f'%dPowerBMTdr[i] 
	#-------------------------------------------------------------------------------	

# end of for loop...iteration to calculate radial distribution of BMT and BET properties.
#-------------------------------------------------------
#-------------------------------------------------------
negative_velocity = False
# Print velocities
if nbrow == 1:
	print
	print "-----------------------------------------------------------"
	print "span   Vin[m/s]    V2[m/s]    V3[m/s]    Vexit[m/s]"
	print "-----------------------------------------------------------"
	for i in range(nsect):
		print '%.10f'%span[i],'%.10f'%Vz1[i],'%.10f'%V2[i],'%.10f'%V3[i],'%.10f'%Vexit[i]
	print
	for i in range (nsect):
		V_inlet = Vz1[i]
		V_exit  = Vexit[i] 
		if float(V_exit/V_inlet) <= 0.01 :
			negative_velocity = True
	###############################################
	if negative_velocity:
		print "Negative Exit Velocity"
#	else:
#		continue
	print	
#-------------------------------------------------------
#-------------------------------------------------------
# Total power and thrust calculation in kW(Trapezoidal integration)
# CP      = trapez_integral(dCP,lambda_r,nsect) 
ThrustBET   = trapez_integral(dThrustBETdr,r)
AForceBMT   = trapez_integral(dAForceBMTdr,r) 
TorqueBET   = trapez_integral(dTorqueBETdr,r) 
TForceBMT   = trapez_integral(dTForceBMTdr,r) 
Power_swirl = trapez_integral(dPower_swirldr,r)
PowerBET    = trapez_integral(dPowerBETdr,r)#OMEGA*TorqueBET
PowerBMT    = trapez_integral(dPowerBMTdr,r)#OMEGA*TForceBMT
print 'ThrustBET,AForceBMT : ',ThrustBET,AForceBMT
#-----------------------------------------------
Vztipv = Vz1[nsect-1]
#
Power_fluid = oneDvalues(Rhub,Rtip,axlIFnew,r_avg,r,Vztipv,rho,TorqueBET,TForceBMT,propeller)
#################################################################
C_P,CT_formula,eta_Prop,eta,eta_HKT = calcEfficiencies(rho,Rhub,Rtip,RPM,Lambda,
													   Vztipv,AForceBMT,
													   ThrustBET,TorqueBET,
													   PowerBET,Power_fluid,propeller)
#---------------------------------------------------------
#---------------------------------------------------------
# Non dimensional chord
#---------------------------------------------------------
for i in range(nsect):
	chrdr_nd[i] = chordr[i]/r[i]
	# print r[i],chordr[i]
values = np.zeros(nsect)
#---------------------------------------------------------
radialPrintProperties(multirow,nsect,span,floss,axlIFnew,
					  angIFnew,values,chrdr_nd,U,Wr,
					  Vthetar,inbeta,stggr,Phi,lambda_r)
#---------------------------------------------------------
#---------------------------------------------------------
# Creating 3DBGB input file:
if case == "CONTRA_AMP":
	case = "CRAMProtor1"
	nbrow1 = nbrow - 1
elif not propeller and nbrow == 2:
	case = "CRHKTr1" 
	nbrow1 = nbrow - 1
else:
	nbrow1 = nbrow
#
bsf = Rtip
xte1 = np.zeros(nsect)
xte1 = create_3dbgbinput(foilname,nbrow1,nblade,nsect,Rhub,Rtip,chrdr_nd,stggr,Wr,case,xte1,propeller,alpha,currentrow,bsf,gap)
#---------------------------------------------------------
#Sweeping lambda from 2 to 11 for eta_HKT
# lambda_sweep_for_etaHKT()
#---------------------------------------------------------
printFinalValues(Rhub,Rtip,Vz1[nsect-1],Lambda,nblade,ThrustBET,AForceBMT,TorqueBET,TForceBMT,
					 PowerBET,PowerBMT,Power_swirl,Power_fluid,C_P,CT_formula,
					 eta_Prop,eta,eta_HKT,poundforce,hp,propeller,RPM,currentrow)

avg_chrd = sum(chordr)/float(len(chordr))
if nsect == 21:
	frontpitch  = 90.0 - abs(stggr[15])
else:
	frontpitch  = 90.0 - abs(stggr[int(0.75*nsect)])	
print
print "span		a		a'		chord[m]	stagger		dTorque/dr[N]		dThrust/dr[N/m]		dPower/dr[W/m]		Vin[m/s]	V2[m/s]		V3[m/s]		Vexit[m/s]	Vthetaslipstream[m/s]"		
for i in range(nsect):
	print '%.10f'%span[i],'%.10f'%axlIFnew[i],'%.10f'%angIFnew[i],'%.10f'%chordr[i],'%.10f'%stggr[i],'%.10f'%dTorqueBETdr[i],'%.10f'%dThrustBETdr[i],'%.10f'%dPowerBETdr[i],'%.10f'%Vz1[i],'%.10f'%V2[i],'%.10f'%V3[i],'%.10f'%Vexit[i],'%.10f'%Vtheta1slipstream[i]		 
print

print "Total_Thrust         (N):",'%.5f'%(ThrustBET)
print "Total_Thrust       (lbf):",'%.5f'%((ThrustBET)*poundforce)
print "Total_Power         (kW):",'%.5f'%((PowerBET)/1000.)  
print "Total_Torque        (Nm):",'%.5f'%TorqueBET
print "OMEGA            (rad/s):",'%.5f'%OMEGA
print "Thrust_over_Power  (s/m):",'%.5f'%((ThrustBET)/(PowerBET))
print "Power_over_Thrust  (m/s):",'%.5f'%((PowerBET)/(ThrustBET))    
#####
print

geom_pitch = 2.0*math.pi*Rtip*math.tan(frontpitch*math.pi/180)

print "------------------------------------------------------------"    
print " Blade Properties "
print "------------------------------------------------------------"
print " Properties   | Units ||    Solo Prop    | "
print " Hub Diameter |   m   || ",'%.5f'%(2.*Rhub)," | "
print " Tip Diameter |   m   || ",'%.5f'%(2.*Rtip)," | "    
print " Vz           |   m/s || ",'%.5f'%Vz1[1]," | "
print " Blade Count  |   -   || ",nblade," | "  
print " Blade Pitch  |  deg  || ",'%.5f'%frontpitch," | "   
print " Geom  Pitch  |   m   || ",'%.5f'%geom_pitch," | " 
print " Aspect Ratio |   -   || ",'%.5f'%(Rtip/avg_chrd)," | "      
print " TSR          |   -   || ",'%.5f'%Lambda," | "  
print " RPM          |   rpm || ",'%.5f'%RPM," | "
print " Tip Mach     |   -   || ",'%.5f'%(OMEGA*Rtip/340.29)," | "
print " Thrust       |   N   || ",'%.5f'%ThrustBET," | "
print " Thrust       |   lbf || ",'%.5f'%(ThrustBET*poundforce)," | "
print " Torque       |   Nm  || ",'%.5f'%TorqueBET," | "
print " Power        |   kW  || ",'%.5f'%(PowerBET/1000.)," | " 
if propeller:
	print " Efficiency   |   -   || ",'%.5f'%eta_Prop," | "  
else: # turbine
	print " Efficiency   |   -   || ",'%.5f'%eta_HKT," | " 
print " Thrust/Power |   N/W || ",'%.5f'%(ThrustBET/PowerBET)," | "   
print " Thrust/Area  |  N/m2 || ",'%.5f'%(ThrustBET/(math.pi*Rtip**2))," | "  
print

# Writing ANALYSIS input file for the design point	
if nbrow == 1:	
	if not spanwise:
		for i in range(nsect):
			Clspan[i] = Cl
			Cdspan[i] = Cd
	arg1 = Rhub,Rtip,ThrustBET,TorqueBET,PowerBET,RPM,Vz1[1],eta_Prop
	arg2 = r,axlIFnew,angIFnew,chordr,REYNfront,stggr,Phi,Clspan,Cdspan,alpha,dAForceBMTdr,dThrustBETdr,dTForceBMTdr,dTorqueBETdr,V2,V3,Vexit
	write_analysis_input(arg1,arg2,case)

##########################################################################################                    
##********************************************************************************
# AFT ROTOR
##********************************************************************************
# bounds for bisection method
phi2_lower = -85.*math.pi/180.
phi2_upper = -1.*math.pi/180.

if nbrow == 2:
	# Co/Counter rotation system activation
	if (RPM > 0.) and (RPM2read > 0.):
		co_rotation = True
	elif (RPM > 0.) and (RPM2read < 0.):
		counter_rotation = True
	print " co_rotation: ", co_rotation
	print " counter_rotation: ", counter_rotation    
	##################   
	currentrow = 2    
	multirow = True   
	##################    
	# Define the gap distance by a multiplier to in V4 = Vz(1-gap*a) or Vz(1+gap*a)
	gap = float(sys.argv[3]) # default is 1.0
	##################
	print
	print "--------------------------------------------------------------------------"
	print "  AFT ROTOR"
	print "--------------------------------------------------------------------------"
	# Radial stations
	#
	r = np.linspace(Rhub2, Rtip2, num=nsect2)
	r2 = r
	# Inlet Velocity
	# Vz2 = Vz(1+2a1)
	print "     r        Vz2"
	for i in range(nsect):
		if propeller:
			Vz2[i] = Vz1[i]*(1 + gap*axlIFnew[i])	
		else: # turbine
			Vz2[i] = Vz1[i]*(1 - gap*axlIFnew[i])
		#######
		print '%.10f'%r[i],'%.10f'%Vz2[i]        
	#
	#-------------------------------------------------
	#Tip Speed Ratio
	if RPM2read == 0:
		OMEGA2 = -OMEGA
		RPM2 = OMEGA2*(60./(2.*math.pi))
		print "OMEGA2 : ", OMEGA2
	else:
		RPM2 = RPM2read
		OMEGA2 = RPM2*(math.pi)/30.
		print "OMEGA2 : ", OMEGA2
	print
	print "RPM2  : ", RPM2
	print
	Lambda2 = (OMEGA2*Rtip2)/Vz2[nsect-1]
	# if Vz2[nsect-1] < 0.:
		# Lambda2 = -Lambda2
	print "TSR2  : ", Lambda2
	print
	#--------------------------------------------------	
	for i in range(nsect2):
		span2CR[i] = (r[i]-Rhub2)/(Rtip2-Rhub2)          				 
		# Lambda at each radial station
		lambda_r2[i] = OMEGA2*r2[i]/Vz2[i]
		# lambda_r2[i] = Lambda2*(r[i]/Rtip2)
		# print lambda_r2[i]
		# Local Wheel speed [m/s]
		U2[i] = lambda_r2[i]*Vz2[i] 
	# end of for loop  
	# print Vz2
	if chordspline2:
		#---------------------------------------------------------
		# Control points for spanwise distribution of chord multiplier 
		#---------------------------------------------------------	
		print
		print " Control points span vs chord mutliplier for ROW 2"
		print "   spanCP2     chrdMPCP2"
		for i in range(ncp2):
			print "  ",'%1.6f'%spancp2[i]," ",'%1.6f'%chrd_multipcp2[i]
		print
		#---------------------------------------------------------
		#---------------------------------------------------------
		chrdmplier2 = subroutines.cubicspline(span2CR,nsect2,chrd_multipcp2,spancp2,ncp2)
		print " Chord multiplier spline vs span for ROW 2: "
		print "    span2CR       chrdMP2"
		for i in range(nsect2):
			print "  ",'%1.6f'%span2CR[i]," ",'%1.6f'%chrdmplier2[i]
		#---------------------------------------------------------
		#---------------------------------------------------------  
	# Spanwise AoA distributition
	if spanwise2:
		print "--------------------------------------------------"
		#---------------------------------------------------------
		# Control points for spanwise distribution of AoA
		#---------------------------------------------------------
		print
		print " Control points span vs AoA"
		print "   spanCP     alpha2CP"
		for i in range(alpha2_ncp):
			print "  ",'%1.6f'%alpha2span_cp[i]," ",'%1.2f'%alpha2_cp[i]
		print
		#---------------------------------------------------------
		#---------------------------------------------------------
		alpha2span = subroutines.cubicspline(span2CR,nsect2,alpha2_cp,alpha2span_cp,alpha2_ncp)
		print " Span vs AoA spline: "
		print "    span       alpha2"
		for i in range(nsect2):
			print "  ",'%1.6f'%span2CR[i]," ",'%1.2f'%alpha2span[i]

		if foilname == 'e857m':
			airfoilname = 'e857'  
		else:
			airfoilname = foilname
		#Clspan, Cdspan, L_over_D = xfoilrun(alpha1span,nsect,REYN,airfoilname) 
		#Cl2span,Cd2span =	Cl_Cd_rangeREYN(foilname,alpha2,propeller,Re_lookup)
		#print Cl2span
		#print Cd2span
	# print
	# print Clspan, Cdspan, L_over_D       
	# print
# stop
#---------------------------------------------------------
#---------------------------------------------------------
	# Defined tolerance
	deftol     = 1.0e-6
	#
	print
	print "Iterative procedure to calculate a and a' for all 2D airfoils:"
	print
	print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	print
	k = 0
	while k < 2: 	
		for i in range(nsect2):
			print "------------------------"
			print "Airfoil   :", i+1
			print "------------------------"
			print
			L_over_D2[i] = Cl2span[i]/float(Cd2span[i])
			if varyingREYN2 and spanwise2 and k == 1:
				print " variable REYN and AoA spanwise"
				alpha2  = alpha2span[i]
				Cl2     = Cl2span[i]
				Cd2     = Cd2span[i]      
			elif varyingREYN2 and k == 1:
				print " variable REYN spanwise and constant AoA"
				Cl2     = Cl2span[i]
				Cd2     = Cd2span[i]                
			elif spanwise2:
				print " variable AoA spanwise "
				alpha2 = alpha2span[i]
				Cl2     = Cl2span[i]
				Cd2     = Cd2span[i]	
			##############	
			if definedchord:
				chordr2[i] = chord2defn[i]	
				internalchord = False 			
			## Bisection method to find a and a'
			chordr2[i]   = chordr2[i]*chrdmplier2[i]
			print "chordr[i] : ", chordr2[i]                   
			#
			print "Cl,Cd before the routine:",Cl2,Cd2			
			# Phirad2[i],axlIFnew2[i],angIFnew2[i],angIFnew12[i],floss2[i],Ca2[i],Cr2[i],sigmar2[i],chordr2[i] = calc_a_aprime_aftrotor(currentrow,co_rotation,counter_rotation,propeller,file_length,
																																	  # chordr2[i],chrdmplier2[i],lambda_r2[i],
																																	  # nbrow,r[i],Rhub2,Rtip2,
																																	  # nblade2,Cl,Cd,gap)
			chordr2[i]  = math.fabs(chordr2[i])
			rcurrent2      = r2[i]
			chord_current2 = chordr2[i]
			Vz_current2    = Vz2[i]
			####
			if gap != 1.0:
				arguments2_in = propeller,nblade2,Lambda2,Rtip2,Vz_current2,(OMEGA2),chord_current2,rcurrent2,Cl2,Cd2     
				Phirad2[i]    = bisection(phi2_lower,phi2_upper,arguments2_in,tol) # phi-function is same as for single rotor            
				####
				if Phirad2[i] == 0.:
					Phirad2[i] = 0.000001  
				if (i == nsect-1): # last section extrapolation
					Phirad2[i] = Phirad2[i-1] + (Phirad2[i-1] - Phirad2[i-2])                
				####        
				Ca2[i],Cr2[i] = Ca_Cr(Cl2,Cd2,Phirad2[i],propeller)       
				floss2[i]     = total_loss(nblade2,lambda_r2[i],(Phirad2[i]),r2[i],Rhub2,Rtip2)	
				sigmar2[i]    = nblade2*chordr2[i]/float(2*math.pi*r2[i])
				Term1aft      = 1./(lambda_r2[i]*math.tan(Phirad2[i]))
				Nr2           = sigmar2[i]*Vz2[i]*Cr2[i]
				Dr2           = 4.*floss2[i]*(OMEGA2)*(r2[i])*(math.sin(Phirad2[i]))**2
				if Dr2 == 0. :
					Dr2 = 0.000001
				Term2aft      = Nr2/Dr2
				bterm2        = (1.0/float(Term1aft + Term2aft))    
				#
				if propeller:
					axlIFnew2[i]  = 1./float(((4.*floss2[i]*(math.sin(Phirad2[i]))**2)/(sigmar2[i]*Ca2[i])) - 1.)
					bterm2        = 1. + axlIFnew2[i]
				else:
					print " define axlIFnew2 for non propeller cases..."
				# axlIFnew2[i]  = bterm2 - 1.
				print 'axlIFnew2, b2-1 :',axlIFnew2[i], bterm2-1
				##
				# if (i == nsect2-1):
					# angIFnew2[i]  = -1.0
					# axlIFnew2[i]  = -1.0
				#
				# print 'bterm2,(lambda_r2[i] x math.tan(Phirad2[i])):',bterm2,(lambda_r2[i]*math.tan(Phirad2[i]))
				term3[i]      = (bterm2)/float(lambda_r2[i]*math.tan(Phirad2[i])) 
				term4         = 1.0 - term3[i] #bterm2-1.0
				if term4 == 0.:
					term4 = 1e-6
				term5         = ((bterm2-1.0)*(bterm2))/(lambda_r2[i]**2)
				# print 'term3,term4,term5 :', term3[i], term4, term5
				print 'term4,Term2aft*bterm2,(Term2aft*(1+axlIFnew2[i])) :',term4,(Term2aft*bterm2),(Term2aft*(1+axlIFnew2[i]))
				angIFnew12[i] = -(0.5/term4)*(term5 + term4**2 + term4)
				angIFnew2[i]  =  angIFnew12[i] + 1 - term4 #((bterm2)*Term2)#   #
				# angIFnew2[i]   = 1.0 - bterm2/(lambda_r2[i]*math.tan(Phirad2[i]))  
				
			elif gap == 1.0:
				Vz_current2   = Vz1[i] # Using Vz1
				angIFnew_current = angIFnew[i]
				axlIF_current = axlIFnew[i]
				arguments2_in = propeller,nblade2,Lambda2,Rtip2,Vz_current2,angIFnew_current,axlIF_current,OMEGA,OMEGA2,chord_current2,rcurrent2,Cl2,Cd2,gap            
				Phirad2[i]    = bisection_closely_placed_rotors(phi2_lower,phi2_upper,arguments2_in,tol) # phi-function is same as for single rotor            
				####
				if Phirad2[i] == 0.:
					Phirad2[i] = 0.000001  
				if (i == nsect-1): # last section extrapolation
					Phirad2[i] = Phirad2[i-1] + (Phirad2[i-1] - Phirad2[i-2])                
				####        
				Ca2[i],Cr2[i] = Ca_Cr(Cl2,Cd2,Phirad2[i],propeller)       
				floss2[i]     = total_loss(nblade2,lambda_r2[i],(Phirad2[i]),r2[i],Rhub2,Rtip2)	
				sigmar2[i]    = nblade2*chordr2[i]/float(2*math.pi*r2[i])
				#Term1aft      = Vz1[i]/(r2[i]*math.tan(Phirad2[i]))
				#Nr2           = sigmar2[i]*Vz1[i]*Cr2[i]
				#Dr2           = 4.*floss2[i]*(r2[i])*(math.sin(Phirad2[i]))**2
				#if Dr2 == 0. :
				#	Dr2 = 0.000001
				#Term2aft      = Nr2/Dr2
				#bterm2        = (OMEGA2 - (3.0*angIFnew[i]*OMEGA))/float(Term1aft + Term2aft)  
				#
				if propeller:  # Equating a1 to a2.
					axlIFnew2[i]  = axlIFnew[i]
					bterm2        = 1. + axlIFnew2[i]
				else:
					print " define axlIFnew2 for non propeller cases..."
				# axlIFnew2[i]  = bterm2 - 1.
				print 'axlIFnew2, b2-1 :',axlIFnew2[i], bterm2-1
				#
				term3[i]     = (Vz1[i]*(1+axlIFnew2[i]))/float(r2[i]*math.tan(Phirad2[i])) 
				angIFnew2[i] = 1 - (2*angIFnew[i]*OMEGA/float(OMEGA2)) - (term3[i]/float(OMEGA2)) #((angIFnew[i]*OMEGA) + OMEGA2 - term3[i])/float(OMEGA2)
			###### end if loop for gap=1 or not case
			print "lambdar2,phirad2,axlIF2,angIF2,angIF12: ",lambda_r2[i],(Phirad2[i]*180./math.pi),axlIFnew2[i],angIFnew2[i],angIFnew12[i]            
			#Calculating final velocities at various axial stations after correcting axlIF1
			#=================================================
			V2CR[i],V3CR[i],V4CR[i],V5CR[i],V6CR[i],VexitCR[i] = velocitycalc_CR(propeller,Vz1[i],axlIFnew[i],axlIFnew2[i],gap)
			#=================================================			
			if (i == nsect2-1):
				V2CR[i]    = V2CR[i-1] + (V2CR[i-1] - V2CR[i-2]) 
				V3CR[i]    = V3CR[i-1] + (V3CR[i-1] - V3CR[i-2]) 
				V4CR[i]    = V4CR[i-1] + (V4CR[i-1] - V4CR[i-2])  
				V5CR[i]    = V5CR[i-1] + (V5CR[i-1] - V5CR[i-2]) 
				V6CR[i]    = V6CR[i-1] + (V6CR[i-1] - V6CR[i-2]) 
				VexitCR[i] = VexitCR[i-1] + (VexitCR[i-1] - VexitCR[i-2])                 
			### End of loop
		#########################################################################################
		#-------------------------------------------------------------------------------
		# end of for loop
		#-------------------------------------------------------------------------------                      
		for i in range(nsect2):      			
			#-------------------------------------------------------------------------------
			### Final Geometry parameters
			#-------------------------------------------------------------------------------
			array2    = np.zeros(nsect2)
			# Phirad2[i] = -math.pi + Phirad2[i]
			array2[i] = angIFnew2[i] - angIFnew12[i] # this is done in order to use the same function for calculating geometry parameters
			if gap == 1.0:
				angular_term = r2[i]*(OMEGA2 - (OMEGA2*angIFnew2[i]) - (2.0*OMEGA*angIFnew[i]))
			else:
				angular_term = 0.0
			Phi2[i],inbeta2[i],pitch2[i],stggr2[i],Wr2[i] = geometry_param(propeller,Phirad2[i],alpha2,Vz2[i],axlIFnew2[i],
																		   array2[i],lambda_r2[i],currentrow,co_rotation,counter_rotation,gap,angular_term)
			# if propeller:
				# stggr2[i] = math.pi-stggr2[i] # Temp fix. Need a pmnt. fix.
		# ## For stggr, removing anomalies in the curve
		# if (i == nsect2-1): # last section extrapolation
			# stggr2[i]  = stggr2[i-1] + (stggr2[i-1] - stggr2[i-2])              
		#-------------------------------------------------------------------------------  
		#-------------------------------------------------------------------------------        
		# end of for loop			
		# Reynolds_Number Calculation
		Re_lookup = np.zeros(nsect2)       
		Re_lookup,REYNrear = Reyn_Calc(Wr2,chordr2,nju,propeller,alpha2,foilname,plot,span2CR,currentrow,k)
		k = k + 1
		print "k:",k    
		#### Discarding very high REYN cases for optimization. Temporary fix.
		for j in range(len(Re_lookup)):
			if Re_lookup[j] > 8.0:
				print " Too large a Reynolds number [e6]:",Re_lookup[j]
				print " Setting REYN to 6e6 for all cases above 8e6..."
				Re_lookup[j] = 6.0
				#varyingREYN2 = False
				#k = 2  
				#break    
		if spanwise2 or varyingREYN2:
			print "k is 1:",k	
			Cl2span,Cd2span =	Cl_Cd_rangeREYN(foilname,alpha2span,propeller,Re_lookup)
			print Cl2span
			print Cd2span
	######################################################################	
	### End of inductionfactor and REYN loop			
	print "Ca2:",Ca2
	print
	print "Cr2:",Cr2
	# average radius
	r_avg = sum(r)/float(len(r))
	print
	print "Average radius:",r_avg," m."
	#-------------------------------------------------------
	# Radial Distribution of properties
	#---------------------------------------------------
	print 'dThrustBET2dr[i],dAForceBMT2dr[i], dTorqueBET2dr[i], dTForceBMT2dr[i],dPowerBET2dr[i],dPowerBMT2dr[i]'     
	for i in range(nsect2):
		Vtheta2slipstream = np.zeros(nsect2)  
		Vtheta2slipstream[i] = -2.0*angIFnew2[i]*OMEGA2*r2[i] - 2.0*angIFnew[i]*OMEGA*r[i]
		#-------------------------------------------------------------------------------
		# Blade Element theory definitions
		#-------------------------------------------------------------------------------
		if gap == 1.0: # Using Vz1
			dmdotdr[i], dThrustBET2dr[i], dTorqueBET2dr[i], dPowerBET2dr[i], Vthetar2[i] = radial_distr_BET(multirow,currentrow,co_rotation,counter_rotation,propeller,rho,Vz1[i],
																										sigmar2[i],axlIFnew2[i],r2[i],Phirad2[i],Ca2[i],Cr2[i],(OMEGA2),nbrow)
		elif gap != 1.0: # Using Vz2
			dmdotdr[i], dThrustBET2dr[i], dTorqueBET2dr[i], dPowerBET2dr[i], Vthetar2[i] = radial_distr_BET(multirow,currentrow,co_rotation,counter_rotation,propeller,rho,Vz2[i],
																										sigmar2[i],axlIFnew2[i],r2[i],Phirad2[i],Ca2[i],Cr2[i],(OMEGA2),nbrow)
		# print  "dThrustBET2dr[i], dTorqueBET2dr[i], dPowerBET2dr[i] :",dThrustBET2dr[i], dTorqueBET2dr[i], dPowerBET2dr[i]                                                                                               
		#---------------------------------
		if (i == nsect2-1):
			Vthetar2[i] = Vthetar2[i-1] + (Vthetar2[i-1] - Vthetar2[i-2])  
		# Power in the swirl downstream of the rotor
		dPower_swirl2dr[i] = 0.5*rho*Vz2[i]*(Vthetar2[i]**2)*2.*math.pi*r2[i] 
		#-------------------------------------------------------------------------------
		# Blade Momentum theory definitions
		#------------------------------------------------------------------------------- 
		if gap == 1.0: # Using Vz1
			newterm    = np.zeros(nsect2)
			newterm[i] = angIFnew[i]*OMEGA 
			dAForceBMT2dr[i], dTForceBMT2dr[i], dPowerBMT2dr[i] = radial_distr_BMT(multirow,currentrow,co_rotation,counter_rotation,propeller,rho,Vz1[i],
																			   axlIFnew2[i],angIFnew2[i],newterm[i],r2[i],floss2[i],(OMEGA2),nbrow,gap)            
		##
		elif gap != 1.0: # Using Vz2
			dAForceBMT2dr[i], dTForceBMT2dr[i], dPowerBMT2dr[i] = radial_distr_BMT(multirow,currentrow,co_rotation,counter_rotation,propeller,rho,Vz2[i],
																			   axlIFnew2[i],angIFnew2[i],angIFnew12[i],r2[i],floss2[i],(OMEGA2),nbrow,gap)
		# print "dAForceBMT2dr[i], dTForceBMT2dr[i], dPowerBMT2dr[i] :",dAForceBMT2dr[i], dTForceBMT2dr[i], dPowerBMT2dr[i]
		# Forcing TorqueBET to be equal to TorqueBMT
		dTForceBMT2dr[i] = dTorqueBET2dr[i]
		dPowerBMT2dr[i]  = dPowerBET2dr[i]
		######
		# Equating the tip dT/dr from both theories
		if (i == nsect2-1):
			dAForceBMT2dr[i] = dThrustBET2dr[i]
		######
		print '%.3f'%dThrustBET2dr[i],'%.3f'%dAForceBMT2dr[i], '%.3f'%dTorqueBET2dr[i], '%.3f'%dTForceBMT2dr[i], '%.3f'%dPowerBET2dr[i],'%.3f'%dPowerBMT2dr[i] 
		#-------------------------------------------------------------------------------	
		# end of for loop...iteration to calculate radial distribution of BMT and BET properties.    
	###########################################
	#-------------------------------------------------------
	#-------------------------------------------------------
	# Print velocities
	print
	print "---------------------------------------------------------------------------------------------"
	print "span     Vin[m/s]     V2[m/s]     V3[m/s]   V4=Vz2[m/s]   V5[m/s]    V6[m/s]     Vexit[m/s]"
	print "---------------------------------------------------------------------------------------------"
	for i in range(nsect):
		print '%.10f'%span[i],'%.10f'%Vz1[i],'%.10f'%V2CR[i],'%.10f'%V3CR[i],'%.10f'%V4CR[i],'%.10f'%V5CR[i],'%.10f'%V6CR[i],'%.10f'%VexitCR[i]
	print
	for i in range (nsect2):
		V_inletCR = Vz1[i]
		V_exitCR  = VexitCR[i] 
		if float(V_exitCR/V_inletCR) <= 0.01 :
			negative_velocity = True
	###############################################
	if negative_velocity:
		print "Negative Exit Velocity"
	print	
	#-------------------------------------------------------
	#-------------------------------------------------------    
	# Total power and thrust calculation in kW(Trapezoidal integration)
	ThrustBET2   = trapez_integral(dThrustBET2dr,r2)
	AForceBMT2   = trapez_integral(dAForceBMT2dr,r2) 	
	TorqueBET2   = trapez_integral(dTorqueBET2dr,r2) 
	TForceBMT2   = trapez_integral(dTForceBMT2dr,r2) 	
	# PowerBMT2    = trapez_integral(dPowerBMT2dr,r) 
	Power_swirl2 = trapez_integral(dPower_swirl2dr,r2)
	PowerBET2    = trapez_integral(dPowerBET2dr,r2)#(OMEGA2*TorqueBET2)
	PowerBMT2    = trapez_integral(dPowerBMT2dr,r2)#(OMEGA2*TForceBMT2)
	#-----------------------------------------------
	#-----------------------------------------------
	print
	Vz2_avg = sum(Vz2)/float(len(Vz2))
	print
	print "Average Vz2: ", Vz2_avg
	print
	# using avg Vz for aft rotor
	Power_fluid2 = oneDvalues(Rhub2,Rtip2,axlIFnew2,r_avg,r2,Vztipv,rho,TorqueBET2,TForceBMT2,propeller)
	######################################################
	C_P2,CT_formula2,eta_Prop2,eta2,eta_HKT2 = calcEfficiencies(rho,Rhub2,Rtip2,RPM2,Lambda2,
													   Vztipv,AForceBMT2,
													   ThrustBET2,TorqueBET2,
													   PowerBET2,Power_fluid2,propeller)
	#---------------------------------------------------------
	#---------------------------------------------------------
	# Non dimensional chord
	#---------------------------------------------------------
	for i in range(nsect2):
		chrdr_nd2[i] = chordr2[i]/r2[i]
		# print r[i],chordr[i]
	#---------------------------------------------------------
	radialPrintProperties(multirow,nsect2,span2CR,floss2,axlIFnew2,
					  angIFnew2,angIFnew12,chrdr_nd2,U2,Wr2,
					  Vthetar2,inbeta2,stggr2,Phi2,lambda_r2)
	#---------------------------------------------------------
	if case == "CONTRA_AMP":
		case = "CRAMProtor2"
	elif not propeller and nbrow == 2:
		case = "CRHKTr2"
	# Creating 3DBGB input file:
	foilname2 = foilname2 + 'CR'
	xte2 = np.zeros(nsect2)
	xte2 = create_3dbgbinput(foilname2,nbrow,nblade2,nsect2,Rhub2,Rtip2,chrdr_nd2,stggr2,Wr2,case,xte1,propeller,alpha2,currentrow,Rtip2,gap)
	#---------------------------------------------------------
	#Sweeping lambda from 2 to 11 for eta_HKT
	# lambda_sweep_for_etaHKT()
	#---------------------------------------------------------
	#---------------------------------------------------------
	# print
	# R2slipstream = Rtip*math.sqrt((1.)/(1+(2.*axlIFnew[nsect-1])))#*(1 + axlIFnew2[nsect2-1])))
	# print "Radius of slip stream downstream for single Propeller: ",R2slipstream, "m"
	print
	print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	print
	print " FRONT ROTOR "
	printFinalValues(Rhub,Rtip,Vz1[nsect-1],Lambda,nblade,ThrustBET,AForceBMT,TorqueBET,TForceBMT,
					 PowerBET,PowerBMT,Power_swirl,Power_fluid,C_P,CT_formula,
					 eta_Prop,eta,eta_HKT,poundforce,hp,propeller,RPM,1)
	##########################################################################################
	print
	print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	print
	print " AFT ROTOR "
	printFinalValues(Rhub2,Rtip2,Vz2[nsect2-1],Lambda2,nblade2,ThrustBET2,AForceBMT2,TorqueBET2,TForceBMT2,
					 PowerBET2,PowerBMT2,Power_swirl2,Power_fluid2,C_P2,CT_formula2,
					 eta_Prop2,eta2,eta_HKT2,poundforce,hp,propeller,RPM2,currentrow)
	#-----------------------------------------------------------------------------------------    
	##########################################################################################
	##########################################################################################
	print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	print
	print "Total Thrust(BET) for BOTH rotors   : ",(ThrustBET + ThrustBET2), "N [",(ThrustBET + ThrustBET2)*poundforce, "lbf ]"
	print "Total Thrust(BMT) for BOTH rotors   : ",(AForceBMT + AForceBMT2), "N [",(AForceBMT + AForceBMT2)*poundforce, "lbf ]"
	print
	print "Total Power(BET)  for BOTH rotors   : ",((PowerBET) + (PowerBET2)), "W [",((PowerBET) + (PowerBET2))/hp, "HP ]"
	print "Total Power(BMT)  for BOTH rotors   : ",((PowerBMT) + (PowerBMT2)), "W [",((PowerBMT) + (PowerBMT2))/hp, "HP ]"
	print
	print "Total_Thrust         (N):",'%.5f'%(ThrustBET + ThrustBET2)
	print "Total_Thrust       (lbf):",'%.5f'%((ThrustBET + ThrustBET2)*poundforce)
	print "Total_Power         (kW):",'%.5f'%((PowerBET + PowerBET2)/1000.)
	print "Thrust_FrontRotor    (N):",'%.5f'%ThrustBET   
	print "Thrust_AftRotor      (N):",'%.5f'%ThrustBET2
	print "Thrust_FrontRotor  (lbf):",'%.5f'%(ThrustBET*poundforce)
	print "Thrust_AftRotor    (lbf):",'%.5f'%(ThrustBET2*poundforce)     
	print "Power_FrontRotor    (kW):",'%.5f'%(PowerBET/1000.)
	print "Power_AftRotor      (kW):",'%.5f'%(PowerBET2/1000.)
	print "Torque_FrontRotor   (Nm):",'%.5f'%TorqueBET
	print "Torque_AftRotor     (Nm):",'%.5f'%TorqueBET2
	print "OMEGA_FrontRotor (rad/s):",'%.5f'%OMEGA
	print "OMEGA_AftRotor   (rad/s):",'%.5f'%OMEGA2
	print "Thrust_over_Power  (s/m):",'%.5f'%((ThrustBET + ThrustBET2)/(PowerBET + PowerBET2))
	print "Thrust/Area       (N/m2): ",'%.5f'%((ThrustBET + ThrustBET2)/(math.pi*Rtip**2))  
	print "Power_over_Thrust  (m/s):",'%.5f'%((PowerBET + PowerBET2)/(ThrustBET + ThrustBET2))
	print
	if propeller:   
		etaCRprop = ((ThrustBET + ThrustBET2)*Vz1[1])/((PowerBET) + (PowerBET2)) 
		print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
		print "Counter Rotating Propeller Efficiency  :",etaCRprop
		print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	else:
		etaCRturb = (PowerBET + PowerBET2)/((ThrustBET + ThrustBET2)*Vz1[1])
		print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
		if gap != 1.0:
			print "V4 = Vz1(1 -",gap,"(a1))" 
		print "Counter Rotating Turbine Efficiency  :",etaCRturb    
		print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	#---------------------------------------------------------
	print
	avg_chrd2 = sum(chordr2)/float(len(chordr2))
	frontpitch  = 90.0 - abs(stggr[15])
	aftpitch  = 90.0 - abs(stggr2[15])
	print
	print "span		a2		a2'		a12' 	chord[m]	stagger		dTorque/dr[N]		dThrust/dr[N/m]		dPower/dr[W/m]		Vin[m/s]	V2[m/s]		V3[m/s]		V4=Vz2[m/s]   V5[m/s]    V6[m/s]     Vexit[m/s]      Vthetaslipstream[m/s]"	
	for i in range(nsect2):
		print '%.10f'%span2CR[i],'%.10f'%axlIFnew2[i],'%.10f'%angIFnew2[i],'%.10f'%angIFnew12[i],'%.10f'%chordr2[i],'%.10f'%stggr2[i],'%.10f'%dTorqueBET2dr[i],'%.10f'%dThrustBET2dr[i],'%.10f'%dPowerBET2dr[i],'%.10f'%Vz1[i],'%.10f'%V2CR[i],'%.10f'%V3CR[i],'%.10f'%V4CR[i],'%.10f'%V5CR[i],'%.10f'%V6CR[i],'%.10f'%VexitCR[i],'%.10f'%Vtheta2slipstream[i]		 
	#####
	print

	aft_geom_pitch = 2.0*math.pi*Rtip*math.tan(aftpitch*math.pi/180)
	
	print "------------------------------------------------------------"    
	print " Blade Properties "
	print "------------------------------------------------------------"
	print " Properties   | Units ||    Front    |     Aft    | "
	print " Hub Diameter |   m   || ",'%.5f'%(2.*Rhub)," | ",'%.5f'%(2.*Rhub2),  " |"
	print " Tip Diameter |   m   || ",'%.5f'%(2.*Rtip)," | ",'%.5f'%(2.*Rtip2),  " |"    
	print " Vz           |   m/s || ",'%.5f'%Vz1[1]," | ",'%.5f'%Vz2[1],  " |"
	print " Blade Count  |   -   || ",nblade," | ",nblade2 , " |"  
	print " Blade Pitch  |   deg || ",'%.5f'%frontpitch," | ",'%.5f'%aftpitch , " |"   
	print " Geom  Pitch  |   m   || ",'%.5f'%geom_pitch," | ",'%.5f'%aft_geom_pitch , " |"      
	print " Aspect Ratio |   -   || ",'%.5f'%(Rtip/avg_chrd)," | ",'%.5f'%(Rtip2/avg_chrd2),  " |"       
	print " TSR          |   -   || ",'%.5f'%Lambda," | ",'%.5f'%Lambda2,  " |"    
	print " RPM          |   rpm || ",'%.5f'%RPM," | ",'%.5f'%RPM2,  " |"
	print " Tip Mach     |   -   || ",'%.5f'%(OMEGA*Rtip/340.29)," | ",'%.5f'%(OMEGA2*Rtip2/340.29),  " |" 
	print " Thrust       |   N   || ",'%.5f'%ThrustBET," | ",'%.5f'%ThrustBET2,  " |"
	print " Thrust       |   lbf || ",'%.5f'%(ThrustBET*poundforce)," | ",'%.5f'%(ThrustBET2*poundforce),  " |"
	print " Torque       |   Nm  || ",'%.5f'%TorqueBET," | ",'%.5f'%TorqueBET2,  " |"   
	print " Power        |   kW  || ",'%.5f'%(PowerBET/1000.)," | ",'%.5f'%(PowerBET2/1000.),  " |"   
	print " Efficiency   |   -   || ",'%.5f'%eta_Prop," | ",'%.5f'%eta_Prop2,  " |"  
	print " Thrust/Power |   s/m || ",'%.5f'%(ThrustBET/PowerBET)," | ",'%.5f'%(ThrustBET2/PowerBET2),  " |"      
	print "Thrust/Area   |  N/m2 || ",'%.5f'%((ThrustBET)/(math.pi*Rtip**2))," | ",'%.5f'%((ThrustBET2)/(math.pi*Rtip2**2)),  " | "
	print

##########################################################################################    
# End of if loop # for row 2
##********************************************************************************
# Calculating different spans for velocities at different axial locations
# Single Row
if nbrow == 1:
	# print "r1SR,Vz1,span1,r1,V2,spanR1,r4SR,Vexit,span4: "
	for i in range(nsect):
		r1SR[i],r4SR[i] = spancalc(Vz1[i],V2[i],Vexit[i],r1[i])
	# span calculation
		span1[i]   = span[i]#(r1SR[i]-Rhub)/(Rtip-Rhub)
		r4SR[0] = Rhub
		span4[i]   = math.fabs(r4SR[i]-Rhub)/(Rtip-Rhub)  
		# print '%.4f'%r1SR[i],'%.4f'%Vz1[i],'%.4f'%span1[i],'%.4f'%r1[i],'%.4f'%V2[i],'%.4f'%span[i],'%.4f'%r4SR[i],'%.4f'%Vexit[i],'%.4f'%span4[i]        
elif nbrow == 2:
	# print "r1CR,Vz1,span1CR,r4CR,V4CR,span4CR,r5,V5CR,span5CR,r7,VexitCR,span7CR: "
	for i in range(nsect2):
		r1CR[i],r4CR[i],r7CR[i] = spancalc_CR(Vz1[i],V2CR[i],V4CR[i],V5CR[i],VexitCR[i],r1[i],r2[i])
	# span calculation
		span1CR[i] = span2CR[i]#(r1CR[i]-Rhub2)/(Rtip-Rhub)
		span4CR[i] = (r4CR[i]-Rhub2)/(Rtip-Rhub) 
		span5CR[i] = (r[i]-Rhub2)/(Rtip-Rhub)   
		r7CR[0] = Rhub
		span7CR[i] = math.fabs(r7CR[i]-Rhub2)/(Rtip-Rhub)    
		print '%.4f'%r1CR[i],'%.4f'%Vz1[i],'%.4f'%span1CR[i],'%.4f'%r4CR[i],'%.4f'%V4CR[i],'%.4f'%span4CR[i],'%.4f'%r2[i],'%.4f'%V5CR[i],'%.4f'%span5CR[i],'%.4f'%r7CR[i],'%.4f'%VexitCR[i],'%.4f'%span7CR[i]
print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print
stop = timeit.default_timer()
print
print "Execution Time: ", '%1.8f'%(stop - start) , "seconds"
print
print "============ END of PROGRAM ============="
print
##########################################################################################
#---------------------------------------------------------------    
#Plotting variables spanwise
#---------------------------------------------------------------   
if plot :
	if nbrow == 1:
		arguments = alpha1span,plot,plotdir,span,span1,span4,Vz1,V2,Vexit,axlIFnew,angIFnew,chordr,stggr,dThrustBETdr,dTorqueBETdr,dPowerBETdr,alpha1span,L_over_D,dAForceBMTdr,dTForceBMTdr,dPowerBMTdr
		REYNrear  = 0.
		# print alpha1span,		
		# Calling plotting routine
		plots(propeller,nbrow,arguments,REYNfront,REYNrear)
	elif nbrow == 2: 
		arguments = alpha1span,alpha2span,plot,plotdir,span,span4CR,span5CR,span7CR,Vz1,V2CR,V4CR,V5CR,VexitCR,axlIFnew,axlIFnew2,angIFnew,angIFnew2,angIFnew12,chordr,chordr2,stggr,stggr2,dThrustBETdr,dThrustBET2dr,dTorqueBETdr,dTorqueBET2dr,dPowerBETdr,dPowerBET2dr,L_over_D,L_over_D2,dAForceBMTdr,dAForceBMT2dr,dTForceBMTdr,dTForceBMT2dr,dPowerBMTdr,dPowerBMT2dr,gap
		# print alpha1span, alpha2span, L_over_D, L_over_D2
		# Calling plotting routine			
		plots(propeller,nbrow,arguments,REYNfront,REYNrear)
		
##########################################################################################		
