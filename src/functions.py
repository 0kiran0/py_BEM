# List of function definitions

import sys, os, errno, subprocess
import math, time
# import psutil, os
import socket
import getpass
import numpy as np
import pylab as py
from scipy.optimize import brentq
from scipy.interpolate import *
from scipy import ndimage
from matplotlib.pyplot import *
import subroutines
import numpy as np
from subprocess import Popen
np.seterr(divide='ignore', invalid='ignore')

# # # Constants declaration
poundforce = 0.224808943 # 1 N = poundforce
hp         = 745.699872 # W

def kill_proc_tree(pid, including_parent=True):    
	parent = psutil.Process(pid)
	for child in parent.children(recursive=True):
		child.kill()
	if including_parent:
		parent.kill()
		
#---------------------------------------------------------------------------------------
def welcome_msg():
	# Displays welcome message
	print "======================================================================="
	print "== pyBEM : Unified Blade Element Momentum Theory code                =="
	print "======================================================================="
	print "== Version: 0.0                                                      =="
	print "==                                                                   =="
	print "== This software comes with ABSOLUTELY NO WARRANTY                   =="
	print "==                                                                   =="
	print "== This is a program which generates spanwise stagger, chord and...  =="
	print "== ... velocities at different axial stations for Wind Turbines,...  =="
	print "== ... HydroKinetic Turbines (HKTs) and Propellers both in solo ...  =="
	print "== ... and contra-rotating configurations.                           =="
	print "==                                                                   =="
	print "== Inputs: Device type, Fluid Properties, Blade parameters like ...  =="
	print "==         ... flow velocity, TSR, hub & tip radii, # of blades ...  =="
	print "==         ... # of blade rows, # of airfoil sections           ...  =="
	print "==         ... Airfoil parameters like REYN, airfoil type, AoA  ...  =="
	print "==         ... Cl, Cd; spanwise distrib. of chord multipliers   ...  =="
	print "==         ... and AoA.                                              =="
	print "==                                                                   =="
	print "== Outputs: Power, Thrust, Efficiency,spanwise induction factors...  =="
	print "==         ... chord, velocites, stagger, dT/dr, dQ/dr dP/dr    ...  =="
	print "==         ... T-Blade3 (3D Blade Geometry) input file                  =="
	print "==                                                                   =="
	print "== ----------------------------by Kiran Siddappaji                   =="
	print "== -------------------------------s2kn@mail.uc.edu                   =="
	print "======================================================================="
#---------------------------------------------------------
def file_len(fname):
  with open(fname) as f:
	for i, l in enumerate(f):
	  pass
  return i + 1
#---------------------------------------------------------
#---------------------------------------------------------
def check(inputfile,string):
	datafile = file(inputfile)
	for line in datafile:
		if string in line:
			return True
	return False

#---------------------------------------------------------
#---------------------------------------------------------
#def readinput():

# Opening the input file
inputfile = sys.argv[1]
print
welcome_msg()
print
print "Input File: ",inputfile
print
file_length = file_len(inputfile)
print "Number of lines: ", file_length
print
f = open(inputfile, "r")
lineread = f.readlines()
# lines = f.read()#.split()
for i in range(file_length):
	print lineread[i]
f.close()
print
###################################################
# BEMT Blade Design Input Parameters [S.I. Units] #
###################################################
device = lineread[3].split()[0] # Splitting the line by space
# print device
case   = lineread[4].split()[0]
##########################################
# Fluid Properties
##########################################    
fluid   = lineread[9].split()[0]
rho     = float(lineread[10].split()[0])
nju     = float(lineread[11].split()[0])
a_sound = float(lineread[12].split()[0])      
##########################################
# Blade Parameters
########################################## 
inputRPM = False
Vz      = float(lineread[17].split()[0])
char1   = lineread[18].split()[0]
if str(char1) == "RPM":
	inputRPM = True
	print "inputRPM :", inputRPM
	RPMin = float(lineread[18].split()[1])
else:
	Lambda  = float(lineread[18].split()[0])
Rhub    = float(lineread[19].split()[0])
Rtip    = float(lineread[20].split()[0])
nbrow   = int(lineread[21].split()[0])
nblade  = int(lineread[22].split()[0])
nsect   = int(lineread[23].split()[0])
###
Vz1     = np.zeros(nsect)  # Variable declaration
Vz1     = [Vz]*nsect
# print Vz1
##########################################
chordr    = np.zeros(nsect) # local chord [m]
alpha1span = np.zeros(nsect) # spanwise alpha [radians]
Clspan    = np.zeros(nsect) # spanwise Cl
Cdspan    = np.zeros(nsect) # spanwise Cd
REYNspan  = np.zeros(nsect) # spanwise REYN
L_over_D  = np.zeros(nsect) # L/D spanwise
##########################################
# Airfoil Parameters
##########################################
REYN     = float(lineread[28].split()[0])
REYNvary = lineread[28].split()[1]
foilname = lineread[29].split()[0]
alpha    = lineread[30].split()[0]
### Check if the alpha, Cl and Cd are defined spanwise
spanwise = False
varyingREYN = False
if (str(alpha) == "spanwise") and (str(REYNvary) == "spanwise"):
	spanwise = True
	varyingREYN = True
	print
	print "AoA is defined spanwise..."
	print "REYN is varying along the blade..."
	print " Corresponding Cl, Cd values are used..."
	print
	Cl    = str(lineread[31].split()[0])
	Cd    = str(lineread[32].split()[0])
elif str(REYNvary) == "spanwise":
	varyingREYN = True
	print
	print "REYN is varying along the blade..."
	print " Corresponding Cl, Cd values are used..."
	print
	alpha = float(lineread[30].split()[0])	
	Cl    = float(lineread[31].split()[0])
	Cd    = float(lineread[32].split()[0])		
elif str(alpha) == "spanwise":
	spanwise = True
	print
	print "AoA is defined spanwise..."
	print
	Cl    = str(lineread[31].split()[0])
	Cd    = str(lineread[32].split()[0])
	print "Cl,Cd :", Cl,Cd
else:
	alpha = float(lineread[30].split()[0])
	Cl    = float(lineread[31].split()[0])
	Cd    = float(lineread[32].split()[0])
##
#Clc    = float(lineread[31].split()[0])
#Cdc    = float(lineread[32].split()[0])
##
# Checking if chord is defined are not in the input
if nbrow > 1:
	rotor    = str(lineread[34].split()[0])
	print "rotor: ",rotor
	print
	if rotor == "Rear":
		no_chrd_defined = True
		print "rotor, no_chrd_defined :",rotor,no_chrd_defined
	elif rotor == "#":
		no_chrd_defined = False
####################
definedchord = False
string = '#chord'
definedchord = check(inputfile,string)
print "definedchord:",definedchord
internalchord = False 
chordspline = False
#######################
if file_length > 33:
	if nbrow == 1:
		rotor    = str(lineread[34].split()[0])
	elif nbrow == 2:
		rotor    = str(lineread[53].split()[0])  
		if rotor == "#Bspline":
			chordspline = True

	print "rotor: ",rotor
	print
#######
if nbrow == 1 and file_length > 33 and rotor == "#Bspline":
	chordspline = True
	ncp           = int(lineread[36].split()[0])
	spancp        = np.ndarray(ncp)
	chrd_multipcp = np.ndarray(ncp)
	for i in range(ncp):
		spancp[i]        = float(lineread[38+i].split()[0])
		chrd_multipcp[i] = float(lineread[38+i].split()[1])	
	internalchord = True			
	cursor = 33 + 5 + ncp
elif nbrow == 1 and file_length > 33 and rotor == "#":
	for i in range(nsect):
		chordr[i] = float(lineread[33 + 2 + i])
		cursor = 35 + nsect
	internalchord = False 
elif nbrow == 2 and no_chrd_defined:
	chordr = np.zeros(nsect)
	cursor = 33    
	internalchord = True
elif nbrow == 1:
	chordr = np.zeros(nsect)
	cursor = 33    
	internalchord = True    
else:
	definedchord = False
	string = '#chord'
	definedchord = check(inputfile,string)
	print "definedchord:",definedchord
print "chordspline: ", chordspline
#######################################################
# Reading alpha spanwise
if spanwise:
	starter1 = "AoARow1"
	alphaline1 = 0
	for num, line in enumerate(lineread, 1):
		if starter1 in line:
			alphaline1 = num - 1
			print 'alphaline1', alphaline1
	num2 = alphaline1+1        
	alpha1_ncp    = int(lineread[num2].split()[0])
	print alpha1_ncp
	alpha1_cp     = np.ndarray(alpha1_ncp)
	alpha1span_cp = np.ndarray(alpha1_ncp)
	for i in range(alpha1_ncp):
		alpha1span_cp[i] = float(lineread[num2+2+i].split()[0])
		alpha1_cp[i]     = float(lineread[num2+2+i].split()[1])    
		print alpha1span_cp[i],alpha1_cp[i]        
#######################################################
#######################################################
# Reading chord spanwise
if definedchord:
	starter1 = "#chord"
	chordline1 = 0
	for num, line in enumerate(lineread, 1):
		if starter1 in line:
			chordline1 = num - 1
			print 'chordline1', chordline1
	num2 = chordline1+2        
	chord1defn     = np.ndarray(nsect)
	for i in range(nsect):
		chord1defn[i] = float(lineread[num2+i].split()[0]) 
		# print chord1defn[i]      
#######################################################
## reading the Rear rotor properties if it is counter-rotating (2 rows)
print "@ cursor_line#: ",cursor
if nbrow == 2:
	cursor = cursor + 5 # empty line + 4 lines of comment in the input file
	print "@ cursor_line#: ",cursor
	Rhub2    = float(lineread[cursor + 0].split()[0])
	Rtip2    = float(lineread[cursor + 1].split()[0])
	nblade2  = int(lineread[cursor + 2].split()[0])
	nsect2   = int(lineread[cursor + 3].split()[0])
	RPM2read = float(lineread[cursor + 4].split()[0])
	#print Rhub2,Rtip2,nblade2,nsect2
	alpha2span = np.zeros(nsect2) # spanwise alpha [radians]
	Cl2span    = np.zeros(nsect2) # spanwise Cl
	Cd2span    = np.zeros(nsect2) # spanwise Cd
	L_over_D2  = np.zeros(nsect) # L/D spanwise
	REYN2span  = np.zeros(nsect2) # spanwise REYN
	######
	cursor   = cursor + 9 # 5 data lines + empty line + 3 comment lines in the input file
	print "@ cursor_line#: ",cursor
	#######
	REYN2     = float(lineread[cursor + 0].split()[0])
	REYN2vary = lineread[cursor + 0].split()[1]
	foilname2 = lineread[cursor + 1].split()[0]
	### Check if the alpha2, Cl2 and Cd2 are defined spanwise
	spanwise2 = False
	varyingREYN2 = False	
	alpha2    = lineread[cursor + 2].split()[0]
	if (str(alpha2) == "spanwise") and (str(REYN2vary) == "spanwise"):
		spanwise2 = True
		varyingREYN2 = True
		print
		print "AoA is defined spanwise..."
		print "REYN is varying along the blade..."
		print " Corresponding Cl, Cd values are used..."
		print
		Cl2    = float(lineread[cursor + 3].split()[0])
		Cd2    = float(lineread[cursor + 4].split()[0])	    
	elif str(REYN2vary) == "spanwise":
		varyingREYN2 = True
		print
		print "REYN is varying along the blade..."
		print " Corresponding Cl, Cd values are used..."
		print
		alpha2 = float(lineread[cursor + 2].split()[0])	
		Cl2    = float(lineread[cursor + 3].split()[0])
		Cd2    = float(lineread[cursor + 4].split()[0])
	elif str(alpha2) == "spanwise":
		spanwise2 = True
		print
		print "AoA is defined spanwise..."
		print
		Cl2    = str(lineread[cursor + 3].split()[0])
		Cd2    = str(lineread[cursor + 4].split()[0])
	else:
		alpha2 = float(lineread[cursor + 2].split()[0])
		Cl2    = float(lineread[cursor + 3].split()[0])
		Cd2    = float(lineread[cursor + 4].split()[0])  
	# print REYN2,foilname2,alpha2,Cl2,Cd2
	#######
	cursor = cursor + 7 # 5 data lines + empty line + 2 comment lines
	print "@ cursor_line#: ",cursor    
	chordr2  = np.zeros(nsect2) # local chord [m]  
	print "cursor is now at: ",cursor
	print str(lineread[54].split()[0]) 
	chordspline2 = False
	if file_length > (cursor-2):
		bspline = str(lineread[53].split()[0]) 
		if bspline == "#Bspline":
			ncp           = int(lineread[55].split()[0])
			spancp        = np.ndarray(ncp)
			chrd_multipcp = np.ndarray(ncp)
			for i in range(ncp):
				spancp[i]        = float(lineread[57+i].split()[0])
				chrd_multipcp[i] = float(lineread[57+i].split()[1])
			internalchord = True
			print "cursor is now at: ",cursor
			cursor = cursor + 2 + ncp
			print "cursor is now at: ",cursor
			# Row 2
			chordspline2   = True
			print "chordspline2: ", chordspline2
			ncp2           = int(lineread[cursor+2].split()[0])
			spancp2        = np.ndarray(ncp2)
			chrd_multipcp2 = np.ndarray(ncp2)
			for i in range(ncp2):
				spancp2[i]        = float(lineread[cursor+4+i].split()[0])
				chrd_multipcp2[i] = float(lineread[cursor+4+i].split()[1])
			internalchord2 = True				
		else:
			for i in range(nsect2):
				# print i
				chordr2[i] = float(lineread[cursor + i])
			internalchord2 = True
	else:
		chordr2[i] = np.zeros(nsect2)  
		internalchord2 = False     	
		
	# print chordr2
	# End of Chord read
	#-----------------------------------------------
	# Reading alpha2 spanwise
	if spanwise:
		starter2 = "AoARow2"
		alphaline2 = 0
		for numm, line in enumerate(lineread, 1):
			if starter2 in line:
				alphaline2 = numm - 1
				print 'alphaline2', alphaline2
		numm2 = alphaline2+1        
		alpha2_ncp    = int(lineread[numm2].split()[0])
		print alpha2_ncp
		alpha2_cp     = np.ndarray(alpha2_ncp)
		alpha2span_cp = np.ndarray(alpha2_ncp)
		for i in range(alpha2_ncp):
			alpha2span_cp[i] = float(lineread[numm2+2+i].split()[0])
			alpha2_cp[i]     = float(lineread[numm2+2+i].split()[1])    
			print alpha2span_cp[i],alpha2_cp[i]     
	#------------------------------------------------      
	#######################################################
	# Reading chord spanwise
	if definedchord:
		starter2 = "#chord"
		chordline2 = 0
		for numm, line in enumerate(lineread, 1):
			if starter2 in line:
				chordline2 = numm - 1
				print 'chordline2', chordline2
		numm2 = chordline2+2        
		chord2defn     = np.ndarray(nsect2+1)
		for i in range(nsect2):
			chord2defn[i] = float(lineread[numm2+i].split()[1]) 
			# print chord2defn[i]      
	#######################################################    
# End of Rear rotor properties read
##########################################
# print device, case
# print fluid
# print rho, nju, a_sound
# print 
# print Vz, Lambda, Rhub, Rtip,  nbrow, nblade, nsect
# if device == "Propeller" or device == "propeller":
	# print rpm_given
# print
# print REYN, foilname, alpha, Cl, Cd
# print chordr
# print
##########################################
	
	# f.close() 
#    return 
#---------------------------------------------------------
# Variable declaration
IFsoln = np.zeros(2)
# Radial stations
xle      = np.zeros(nsect)
xte      = np.zeros(nsect)
r        = np.ndarray(nsect)
r1       = np.ndarray(nsect)
r1SR     = np.ndarray(nsect)
r4SR     = np.ndarray(nsect)
span     = np.ndarray(nsect)
span1    = np.ndarray(nsect)
span4    = np.ndarray(nsect)
chrdmplier = np.ndarray(nsect)
# span     = np.zeros(nsect)
# r        = np.zeros(nsect)
rnd      = np.zeros(nsect)
lambda_r = np.zeros(nsect)
U        = np.zeros(nsect)
dmdotdr  = np.zeros(nsect)
# Variable arrays Definition
axlIFin   = np.zeros(nsect)
angIFin   = np.zeros(nsect)
axlIFnew  = np.zeros(nsect)
angIFnew  = np.zeros(nsect)
axlIFcrit = 0.3
# Velocity variables
V2        = np.zeros(nsect)
V3        = np.zeros(nsect)
Vexit     = np.zeros(nsect)
#
Nr        = np.zeros(nsect)
Dr        = np.zeros(nsect)
Nr1       = np.zeros(nsect)
Dr1       = np.zeros(nsect)
Nr2       = np.zeros(nsect)
Dr2       = np.zeros(nsect)
Nr3       = np.zeros(nsect)
Dr3       = np.zeros(nsect)
#
term1     = np.zeros(nsect)
term2     = np.zeros(nsect)
termA     = np.zeros(nsect)
termB     = np.zeros(nsect)
term3     = np.zeros(nsect)
term4     = np.zeros(nsect)
term5     = np.zeros(nsect)
term6     = np.zeros(nsect)
Y1        = np.zeros(nsect)
Y2        = np.zeros(nsect)
Y3        = np.zeros(nsect)
#
b         = np.zeros(nsect)
b1        = np.zeros(nsect)
Ca        = np.zeros(nsect)
Cr        = np.zeros(nsect)
CT        = np.zeros(nsect) # Coefficient of Thrust
CP        = np.zeros(nsect) # Coefficient of Power
#
# if device == "Turbine" or device == "turbine":
# chordr    = np.zeros(nsect) # local chord [m]
##
chrdr_nd  = np.zeros(nsect) # local chord/r[i] [n.d]
chrdaxl   = np.zeros(nsect) # axial chord local
chrdaxlnd = np.zeros(nsect) # axial chord non dimensional
chrdactl  = np.zeros(nsect) # actual chord [m]
sigmar    = np.zeros(nsect) # local soliditiy
Wr        = np.zeros(nsect) # local relative velocity
Vthetar   = np.zeros(nsect) # local swirl
mrel      = np.zeros(nsect) # relative mach number at inlet
Phi       = np.zeros(nsect)
Phirad    = np.zeros(nsect)
inbeta    = np.zeros(nsect)
twist     = np.zeros(nsect)
pitch     = np.zeros(nsect)
stggr     = np.zeros(nsect)
#
fhubloss  = np.zeros(nsect)
ftiploss  = np.zeros(nsect)  
ftiploss_corrected  = np.zeros(nsect)
floss     = np.zeros(nsect)
# Power and Axial Force from Blade Element theory
PowerBE   = np.zeros(nsect)
AxForceBE = np.zeros(nsect)
# Power and Axial Force from Blade Momentum theory
PowerBM   = np.zeros(nsect)
AxForceBM = np.zeros(nsect)
##
dCP        = np.zeros(nsect)
dTorqueBETdr = np.zeros(nsect)
dThrustBETdr = np.zeros(nsect)
dTForceBMTdr = np.zeros(nsect)
dAForceBMTdr = np.zeros(nsect)
##
dPowerBETdr = np.zeros(nsect)
dPowerBMTdr = np.zeros(nsect)
dPower_swirldr = np.zeros(nsect)
PowerBMT  = PowerBET = 0.
AForceBMT = ThrustBET = 0.
TorqueBET = TForceBMT = 0.
Power_swirl = 0.
##
L       = np.zeros(10)
eta_HKTrange = np.zeros(10)
##
f1 = np.zeros(nsect)
minREYN = 0.
maxREYN = 0.
avgREYN = 0.
REYNfront  = np.zeros(nsect) 
#Re_lookup = np.zeros(nsect)
#---------------------------------------
#---------------------------------------
# REAR ROTOR VARIABLES
#---------------------------------------
#---------------------------------------
# Radial stations
if nbrow == 2 :
	xle2      = np.zeros(nsect2)
	xte2      = np.zeros(nsect2)
	span1CR   = np.zeros(nsect2)
	span2CR   = np.zeros(nsect2)  
	span4CR   = np.zeros(nsect2)    
	span5CR   = np.zeros(nsect2)    
	span7CR   = np.zeros(nsect2)    
	r         = np.zeros(nsect2)
	r2         = np.zeros(nsect2)    
	r1CR      = np.zeros(nsect2)    
	r4CR      = np.zeros(nsect2)    
	r7CR      = np.zeros(nsect2)      
	rnd       = np.zeros(nsect2)
	lambda_r2 = np.zeros(nsect2)
	Vz2       = np.zeros(nsect2)
	U2        = np.zeros(nsect2)
	dmdotdr   = np.zeros(nsect2)
	# Variable arrays Definition
	axlIFin2   = np.zeros(nsect2)
	angIFin2   = np.zeros(nsect2)
	angIFin12  = np.zeros(nsect2)
	axlIFnew2  = np.zeros(nsect2)
	angIFnew2  = np.zeros(nsect2)
	angIFnew12 = np.zeros(nsect2)
	axlIFcrit  = 0.3
	# Velocity variables
	V2CR      = np.zeros(nsect2)
	V3CR      = np.zeros(nsect2)
	V4CR      = np.zeros(nsect2)
	V5CR      = np.zeros(nsect2)
	V6CR      = np.zeros(nsect2)
	VexitCR   = np.zeros(nsect2)    
	#
	Nrr        = np.zeros(nsect2)
	Drr        = np.zeros(nsect2)
	Nrr1       = np.zeros(nsect2)
	Drr1       = np.zeros(nsect2)
	Nrr2       = np.zeros(nsect2)
	Drr2       = np.zeros(nsect2)
	Nrr3       = np.zeros(nsect2)
	Drr3       = np.zeros(nsect2)
	#
	termm1     = np.zeros(nsect2)
	termm2     = np.zeros(nsect2)
	termmA     = np.zeros(nsect2)
	termmB     = np.zeros(nsect2)
	termm3     = np.zeros(nsect2)
	termm4     = np.zeros(nsect2)
	termm5     = np.zeros(nsect2)
	termm6     = np.zeros(nsect2)
	#
	b2         = np.zeros(nsect2)
	Ca2        = np.zeros(nsect2)
	Cr2        = np.zeros(nsect2)
	CT2        = np.zeros(nsect2) # Coefficient of Thrust
	CP2        = np.zeros(nsect2) # Coefficient of Power
	#
	# if device == "Turbine" or device == "turbine":
	# chordr    = np.zeros(nsect2) # local chord [m]
	##
	chrdr_nd2  = np.zeros(nsect2) # local chord/r[i] [n.d]
	chrdaxl2   = np.zeros(nsect2) # axial chord local
	chrdaxlnd2 = np.zeros(nsect2) # axial chord non dimensional
	chrdactl2  = np.zeros(nsect2) # actual chord [m]
	sigmar2    = np.zeros(nsect2) # local soliditiy
	Wr2        = np.zeros(nsect2) # local relative velocity
	Vthetar2   = np.zeros(nsect2) # local swirl
	mrel2      = np.zeros(nsect2) # relative mach number at inlet
	Phi2       = np.zeros(nsect2)
	Phirad2    = np.zeros(nsect2)
	inbeta2    = np.zeros(nsect2)
	twist2     = np.zeros(nsect2)
	stggr2     = np.zeros(nsect2)
	pitch2     = np.zeros(nsect2)
	#
	fhubloss2  = np.zeros(nsect2)
	ftiploss2  = np.zeros(nsect2)  
	ftiploss_corrected2  = np.zeros(nsect2)
	floss2     = np.zeros(nsect2)
	# Power and Axial Force from Blade Element theory
	PowerBE2   = np.zeros(nsect2)
	AxForceBE2 = np.zeros(nsect2)
	# Power and Axial Force from Blade Momentum theory
	PowerBM2   = np.zeros(nsect2)
	AxForceBM2 = np.zeros(nsect2)
	##
	dCP2          = np.zeros(nsect2)
	dTorqueBET2dr = np.zeros(nsect2)
	dThrustBET2dr = np.zeros(nsect2)
	dTForceBMT2dr = np.zeros(nsect2)
	dAForceBMT2dr = np.zeros(nsect2)
	##
	dPowerBET2dr = np.zeros(nsect2)
	dPowerBMT2dr = np.zeros(nsect2)
	dPower_swirl2dr = np.zeros(nsect2)
	PowerBMT2  = PowerBET2  = 0.
	AForceBMT2 = ThrustBET2 = 0.
	TorqueBET2 = TForceBMT2 = 0.
	Power_swirl2 = 0.
	REYNrear  = np.zeros(nsect2) 
	term3 = np.zeros(nsect2)

#---------------------------------------------------------------------------------------
def make_sure_path_exists(path):
	try:
		os.makedirs(path)
	except OSError as exception:
		if exception.errno != errno.EEXIST:
			raise
#---------------------------------------------------------------------------------------
# Function to calculate TSR
def calcTSR(RPM,Vz,Rtip):
	# Vz [m/s]; Rtip [m]
	TSR = (math.pi*RPM*Rtip)/(30.*Vz)
	print
	print "The TSR is :",TSR
	print
	return TSR
#---------------------------------------------------------------------------------------
# Function to calculate RPM
def calcRPM(TSR,Vz,Rtip):
	# Vz [m/s]; Rtip [m]
	RPM = (30.*Vz*TSR)/(math.pi*Rtip)
	print
	print "The RPM is :",RPM
	print
	return RPM

#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------

# Calculates radius for a specific power requirement
def radius_for_Power():

	# declare_Constants();
	# # #Power at radius R
	# # Power1 = 0.5*rho*(math.pi*Radius**2)*(Vz**3)*eta_turbine*eta_drivetrain
	# # # radius for 1MW Power
	Power  = 1.0 # [W]
	multiplier     = 1.0e6
	coeff_Power    = 0.35*0.8
	eta_drivetrain = 0.80
	radius = math.sqrt(2.*Power*multiplier/(rho*math.pi*(Vz**3)*coeff_Power*eta_drivetrain))
	print
	Power_2mradius = 0.5*rho*(math.pi*4.)*(Vz**3)*coeff_Power*eta_drivetrain
	print "Power_2mradius",Power_2mradius/multiplier, "MW"
	print
	print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	print
	print 'rho,Vz:',rho,Vz
	print
	print "Radius for", Power, "MW Power: ", radius, "meters."
	print "Radius for", Power, "MW Power: ", radius*100., "centimeters."
	print
	print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	print
	#---------------------------------------------------------
	return radius
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------	
def trapez_integral(y,x):
	#integral (a,b):ydx
	trapez = 0.
	for i in range(1,y.size):
		trapez = trapez + 0.5*(x[i] - x[i-1])*(y[i-1] + y[i])

	return trapez
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------	
def summation(y):
	
	sum = 0.
	for i in range(1,y.size):
		sum += (y[i-1] + y[i])

	return sum
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
	
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Reynolds Number calculation
def Reyn_Calc(Vz,chrdactl,nju_fluid,propeller,alpha,airfoil,plot,span,currentrow,k):

	# Re = Vz*chrdactl/nju_fluid
	nrange   = 1#0
	#Vz       = np.zeros((chrdactl.size,nrange))
	Re_range  = np.zeros((chrdactl.size,nrange))
	Re_lookup = np.zeros(chrdactl.size)
	REYN      = np.zeros(chrdactl.size)
	# Re calculation
	print
	print "Reynolds number calculation: "
	for i in range(nrange):
		print
		#print "Vz: ",Vz_range, " m/s"
		print "actual_chrd [m]  REYN    REYN     REYN"
		for j in range(chrdactl.size):
			Re_range[j][i] =  abs(Vz[i]*chrdactl[j]/nju_fluid)           
			if Re_range[j][i] > 10e5:
				base = -5
				Re_lookup[j]   = float(round(int(Re_range[j][i]),base))/10e5     
			elif Re_range[j][i] < 5000.0:
				base = -3
				Re_lookup[j]   = float(round(int(Re_range[j][i]),base))/10e5                 
			else:
				base = -4
				Re_lookup[j]   = float(round(int(Re_range[j][i]),base))/10e5                
			#
			print '%2.8f'%chrdactl[j],'%7.4f'%Re_range[j][i],int(round(int(Re_range[j][i]),base)),'%.3E'%(int(round(int(Re_range[j][i]),base))),Re_lookup[j]
			REYN[j] = abs(Re_range[j][i])
	print
	minREYN = float(min(Re_range))
	maxREYN = float(max(Re_range))
	avgREYN = float(np.mean(Re_range))
	print "Minimum REYN: ", minREYN
	print "Maximum REYN: ", maxREYN
	print "Average REYN: ", avgREYN 
	print
	# print REYN
	# Plotting REYN vs span
	# if plot and k == 1:
		# figname = 'REYNrow'+ str(currentrow) +'.png'
		# fignum = 1
		# py.figure(fignum, figsize=(10, 6))    
		# py.plot(Re_range,span,'--o',color = 'red',markersize=7,mfc ='none',mew=2,lw=2)#
		# py.title("REYN vs span for row "+ str(currentrow))
		# py.xlabel("REYN")
		# py.ylabel("span")
		# leg = py.legend(['REYN row '+ str(currentrow)], loc=2, bbox_to_anchor=(1.05, 1),
				  # borderaxespad=0.)
		# py.savefig(figname,bbox_extra_artists=(leg,), bbox_inches='tight')	
	# Cl_Cd_rangeREYN(airfoil,alpha,propeller,Re_lookup)
	
	return Re_lookup,REYN
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Running Xfoil for a range of REYN and interpolating them to get Cl and Cd at reqd REYN
def Cl_Cd_rangeREYN(airfoil,alphaspan,propeller,Re_lookup):

	#--------------------------------------------------------
	# Look up Cl and Cd at the specified alpha for the integer part of the REYN from the table for the specified airfoil
	# print Re_lookup
	username = getpass.getuser()	
	# print
	Clspan = np.zeros(len(Re_lookup))
	Cdspan = np.zeros(len(Re_lookup))
	for i in range(len(Re_lookup)):
		alpha = alphaspan[i]
		###########
		if airfoil == "clarky":
			if os.name == 'nt':
				REYNdir = "C:/Users/" + username + "/Dropbox/Hydroturbine/GitHydroBEMT/Propeller/UC/clarky_REYNrange/clarky_10k_1000k_0.01alpha_increment/"        
			elif os.name == 'posix':
				REYNdir = "/home/kiran/Counter_Rotating_Propellers/clarky_REYNrange/clarky_10k_1000k_0.01alpha_increment/"        
			REYNfilename = 'CLARK Y AIRFOIL_T1_Re' + '%1.3f'%Re_lookup[i] + '_M0.00_N9.0.txt'
		elif airfoil == "naca4415":
			if os.name == 'nt':
				REYNdir = "C:/Users/" + username + "/Dropbox/Hydroturbine/GitHydroBEMT/NACA4415/"       	
			elif os.name == 'posix':
				REYNdir = "/home/kiran/HydroKineticTurbines/NACA4415/"
			REYNfilename = 'NACA 4415_T1_Re' + '%1.3f'%Re_lookup[i] + '_M0.00_N9.0.txt'  		
		elif airfoil == "s809m":
			if os.name == 'nt':
				REYNdir = "C:/Users/" + username + "/Dropbox/Hydroturbine/GitHydroBEMT/S809/"       	
			elif os.name == 'posix':
				REYNdir = "/home/" + username + "/Dropbox/Hydroturbine/GitHydroBEMT/S809/"
			REYNfilename = "NREL's S809 Airfoil_T1_Re" + '%1.3f'%Re_lookup[i] + '_M0.00_N9.0.txt'  						
		elif airfoil == "ep857m":
			if os.name == 'nt':
				REYNdir = "C:/Users/" + username + "/Dropbox/Hydroturbine/GitHydroBEMT/HKT/case_runs/e857m_REYNrange/"       
				# REYNdir  = "C:/Users/user/Dropbox/Hydroturbine/GitHydroBEMT/HKT/case_runs/e857m_REYNrange	/"			
			elif os.name == 'posix':
				REYNdir = "/home/kiran/HydroKineticTurbines/e857m_REYNrange/"
			REYNfilename = 'EPPLER 857 AIRFOIL_T1_Re' + '%1.3f'%Re_lookup[i] + '_M0.00_N9.0.txt'    
		else:
			print "Need an airfoil name to be set in functions.py, def Cl_Cd_rangeREYN, line 190..."
		#######
		REYNfile = REYNdir + REYNfilename   
		#print "-------------------------------------------"
		file_length = file_len(REYNfile)
		f = open(REYNfile, "r")
		lineread = f.readlines()
		# print REYNfilename
		# print "-------------------------------------------"		
		#print "Number of lines: ", file_length      
		if propeller:
			keyword  = str('%1.2f'%(math.fabs(alpha)))        
			keyword2 = str('%1.2f'%(math.fabs(alpha)+0.01))
			keyword3 = str('%1.2f'%(math.fabs(alpha)-0.01))
		else:
			keyword  = str('%1.2f'%((alpha)))        
			keyword2 = str('%1.2f'%((alpha)+0.01))
			keyword3 = str('%1.2f'%((alpha)-0.01))       
		#print keyword, keyword2
		for num, line in enumerate(lineread, 1):
			if keyword in line:
				# print line
				linelist  = line.split()
				# break
			elif keyword2 in line:
				# print line
				linelist  = line.split()  
			elif keyword3 in line:
				# print line
				linelist  = line.split()                 
			#else:
				#print " The keyword ", keyword, " not found."
				
		f.close()
		# print       
		if propeller:
			Clspan[i] = abs(float(linelist[1]))      
		else:
			Clspan[i] = (float(linelist[1]))
		Cdspan[i] = linelist[2]
		# print 'alpha, Cl, Cd :',alpha,Clspan[i],Cdspan[i]		
	return Clspan,Cdspan
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
def myround(x, base):
	return int(base * round(float(x)/base))

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

# Defining a function which writes the 3dbgb input file taking calculated values from
# induction_factors.py code
def create_3dbgbinput(casename,nbrow,nblades,nsect,Rhub,Rtip,chrdr_nd,stagger,Wr,case,prevXTE,propeller,alpha,currentrow,bsf,axialgap):
	print
	print
	print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	print
	print "Creating 3dbgbinput file..."
	print
	print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	#
	# Variable declaration
	r    = np.zeros(nsect)
	span = np.zeros(nsect)
	#
	#Blade scale
	# Linear interpolation of radius from HUB to TIP
	r    = np.linspace(Rhub, Rtip, num=nsect)
	rnd  = r/Rtip
	#xle  = np.linspace(2.0, 2.25, num=nsect)
	xle = np.zeros(nsect)
	print "xte before:",prevXTE
#	axial_offset = (prevXTE[0] + 0.1*chrdr_nd[0])
	if nbrow == 2 and currentrow == 2:
		axial_offset = (prevXTE[0] + (0.5*axialgap*Rtip))
	# constant radius construction lines/streamlines
	xsl  = np.zeros(20)
	#
	ctrlspan = np.zeros(5)
	ctrlspan = np.linspace(0.0, 1.00, num=5)
	#
	for i in range(nsect):
	# 
		span[i]      = (r[i]-Rhub)/(Rtip-Rhub) 
		chrdaxl[i]   = math.fabs(math.cos(np.radians(stagger[i]))*(chrdr_nd[i]*r[i]))
		chrdaxlnd[i] = chrdaxl[i]/Rtip
		chrdactl[i]  = chrdr_nd[i]*r[i]
		if currentrow == 1 and casename == 'clarky':
			xle[i] = 2.00 - (0.42)*chrdaxlnd[i] # quick fix for CLARKY AIRFOILS to be stacked at CG
		elif currentrow == 1:
			xle[i] = 2.00		
		elif currentrow == 2:
			xle[i] = axial_offset - (0.42)*chrdaxlnd[i] # quick fix for CLARKY AIRFOILS to be stacked at CG       
		xte[i]       = xle[i] + chrdaxlnd[i]
	#
	# Reynolds_Number Calculation
	# Reyn_Calc(Wr,chrdactl,nju,propeller,alpha,casename)
	#
	xsl  = np.linspace((xle[1] - 1),(xle[1] + 1),num=20)
	#
	if casename == 'clarky' :
		casename = 'clark'
	print
	print "casename,nblades,nsect,Rtip,Rhub"
	print
	print casename,nblades,nsect,Rtip,Rhub
	print
	print "Geometric properties :"
	print
	print "------------------------------------------------------------"
	print "  r    axchrd_nd    chrdr_nd      chrdactl[m]    stagger"
	print "------------------------------------------------------------"
	print 
	# for i in range(nsect-1,1,-1):
	for i in range(nsect):
		print '%.10f'%r[i],'%.10f'%chrdaxlnd[i],'%.10f'%chrdr_nd[i],'%.10f'%chrdactl[i],'%.10f'%stagger[i]    
	# print '%.10f'%span[0],'%.10f'%chrdr_nd[0],'%.10f'%stagger[0]
	
	#-------------------------------------------------------------------
	#f.write("" + '\n')
	#f.write(str() + '\n')
	#
	f = open("3dbgbinput." + str(currentrow) + "." + str(case) + ".dat", "w")
	f.write("Input parameters (version 1.1)" + '\n')
	f.write("    " + str(casename) + '\n')
	f.write(" Blade row #:" + '\n')
	f.write("    " + str(currentrow) + '\n')
	f.write(" Number of blades in this row:" + '\n')
	f.write("    " + str(nblades) + '\n')
	f.write(" Blade Scaling factor (mm):" + '\n')
	f.write("    " + str(bsf*1000) + '\n')
	f.write(" Number of streamlines:" + '\n')
	f.write("    " + str(nsect) + '\n')
	f.write(" Angles in the input file (0=Beta_z (default),1=Beta_r):" + '\n')
	f.write("    " + str(0) + '\n')
	f.write(" Airfoil camber defined by curvature control (0=no,1=yes):" + '\n')
	f.write("    " + str(0) + '\n')
	f.write(" Airfoil Thickness distribution (0=Wennerstrom,1=Spline):" + '\n')
	f.write("    " + str(0) + '\n')
	f.write(" Airfoil Thickness multiplier (0=no,1=yes):" + '\n')
	f.write("    " + str(0) + '\n')
	f.write(" Airfoil LE defined by spline (0=no,1=yes):" + '\n')
	f.write("    " + str(0) + '\n')
	f.write(" Non-dimensional Actual chord (0=no,1=yes,2=spline):" + '\n')
	f.write("    " + str(1) + '\n')
	f.write(" Sectionwise properties:" + '\n')
	f.write(" J      in_Beta     out_Beta     mrel_in      chord      t/c_max     Incidence     Deviation    Sec. Flow Angle" + '\n')
	for i in range(nsect):
		# stagger is negative in T-Blade3 for Turbines (convention)
		if not propeller:
			stagger[i] = -stagger[i]
		#		
		mrel[i] = Wr[i]/a_sound
		#f.write(str(i + stagger[i] + stagger[i] + mrel[i] + chrdr_nd[i] + "0.2000" + "0.0000" + "0.0000" + "0.0000") + '\n')
		f.write(('%02d'%(i+1)) + "   ")
		f.write(('%2.8f'%stagger[i]) + "  ")
		f.write(('%2.8f'%stagger[i]) + "  ")
		f.write(('%2.8f'%mrel[i]) + "  ")
		f.write(('%2.8f'%chrdr_nd[i]) + "  ")
		f.write(('%2.8f'%0.2000) + "  ")
		f.write(('%2.8f'%0.0000) + "  ")
		f.write(('%2.8f'%0.0000) + "  ")
		f.write(('%2.8f'%0.0000) + "\n")    
	f.write('\n')
	f.write(" LE / TE curve (x,r) definition :" + '\n')
	f.write(" Number of Curve points :" + '\n')
	f.write("    " + str(nsect) + '\n')
	f.write("   xLE          rLE           xTE          rTE" + '\n')
	for i in range(nsect):
		f.write("    " + ('%2.8f'%xle[i]) + "  ")
		f.write(('%2.8f'%rnd[i]) + "  ")
		f.write(('%2.8f'%xte[i]) + "  ")
		f.write(('%2.8f'%rnd[i]) + "\n")
	f.write('\n') 
	f.write(" # Airfoil type and Variable Radial Stacking information.         #" + '\n')
	f.write(" # stack_u: % chord stack (0.00 to 100.00).                       #" + '\n')
	f.write(" # stack_v: % below or above meanline stack (-100.00 to +100.00). #" + '\n')
	f.write(" # Use +200 for stacking on airfoil area centroid.                #" + '\n')
	f.write(" Variable Radial stacking (0=no,1=yes):" + '\n')
	f.write("    " + ('%01d'%0) + '\n')    
	f.write(" J   type |stk_u |stk_v |umxthk |lethk |tethk  |Jcells(Grid:4n+1) |eta_ofst(<=10){%thkc/Jmax}  |BGgrid(0=no,1=yes) |" + '\n')
	for i in range(nsect):
		f.write('%02d'%(i+1) + "   ")
		f.write(str(casename) + "  ")
		f.write(('%2.2f'%25.000) + "  ")
		f.write(('%2.2f'%00.300) + "  ")
		f.write(('%2.2f'%00.000) + "  ")
		f.write(('%2.2f'%00.020) + "  ")
		f.write(('%2.2f'%00.020) + "  ")
		f.write(('%02d'%15) + "  ")
		f.write(('%02d'%10) + "  ")
		f.write(('%01d'%0)  + '\n')
	f.write('\n')
	f.write(" Control table for blending section variable:" + '\n')
	f.write("           5           0           0" + '\n')    
	f.write("       span                       bf1         bf2" + '\n')
	for i in range(5):
		f.write("  " + ('%2.15f'%ctrlspan[i]) + "           ")
		f.write('%01d'%1  + "           ")
		f.write('%01d'%0  + '\n')
	f.write('\n')
	f.write(" Stacking axis as a fraction of chord(2.=centroid):" + '\n')
	f.write("   " + str('030000')  + '\n')
	f.write('\n')
	f.write(" Control points for delta_m:" + '\n')
	f.write("           " + str(5)  + '\n')
	f.write(str('        span                   delta_m')  + '\n')
	for i in range(5):
		f.write("  " + ('%2.15f'%ctrlspan[i]) + "     ")
		f.write("  " + str('0.000000000000000') + '\n')
	f.write('\n')
	f.write(" Control points for delta_theta:" + '\n')
	f.write("           " + str(5)  + '\n')
	f.write(str('        span                   delta_theta')  + '\n')
	for i in range(5):
		f.write("  " + ('%2.15f'%ctrlspan[i]) + "     ")
		f.write("  " + str('0.000000000000000') + '\n')
	f.write('\n')
	f.write(" Control points for in_beta*:" + '\n')
	f.write("           " + str(5)  + '\n')
	f.write(str('        span                   in_beta*')  + '\n')
	for i in range(5):
		f.write("  " + ('%2.15f'%ctrlspan[i]) + "     ")
		f.write("  " + str('0.000000000000000') + '\n')
	f.write('\n')
	f.write(" Control points for out_beta*:" + '\n')
	f.write("           " + str(5)  + '\n')
	f.write(str('        span                   out_beta*')  + '\n')
	for i in range(5):
		f.write("  " + ('%2.15f'%ctrlspan[i]) + "     ")
		f.write("  " + str('0.000000000000000') + '\n')
	f.write('\n')
	f.write(" Control points for chord:" + '\n')
	f.write("           " + str(5)  + '\n')
	f.write(str('        span                   chord')  + '\n')
	for i in range(5):
		f.write("  " + ('%2.15f'%ctrlspan[i]) + "     ")
		f.write("  " + str('0.000000000000000') + '\n')  
	f.write('\n')
	f.write(" Control points for tm/c:" + '\n')
	f.write("           " + str(5)  + '\n')
	f.write(str('        span                   tm/c')  + '\n')
	for i in range(5):
		f.write("  " + ('%2.15f'%ctrlspan[i]) + "     ")
		f.write("  " + str('0.000000000000000') + '\n')
	f.write('\n')
	f.write(" Hub offset" + '\n')
	f.write(" " + str('0.000000000000000') + '\n')
	f.write(" Tip offset" + '\n')   
	f.write(" " + str('0.000000000000000') + '\n')
	f.write('\n')
	f.write("  Streamline Data" + '\n')    
	f.write("  x_s      r_s" + '\n')     
	for i in range(nsect):
		for j in range(20):
			f.write(('%2.8f'%xsl[j]) + "  ")
			f.write(('%2.8f'%rnd[i]) + '\n')
		f.write('0 0' + '\n')
	# f.write('\n')
	# End of File
	
	f.close
	print "xte:",xte
	return xte
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Function for varying eta_HKT with lambda
def lambda_sweep_for_etaHKT():

	##
	L       = np.linspace(2,11,num=10)
	print L
	eta_HKTrange = [43.966,51.724,55.476,57.436,58.493,59.049,59.305,59.368,59.302,59.144]
	
	fignum = 1
	py.figure(fignum, figsize=(10, 6))    
	py.plot(L,eta_HKTrange,'-o',color = 'blue')
	py.title("eta_HKT vs Lambda (TSR) ")
	py.xlabel("TSR (Tip Speed Ratio)")
	py.ylabel(" eta_HKT ")
	# py.axis('equal')
	# py.axis([0.0,1.0,0.0,10.0])
	py.legend(['eta_HKT'], loc='best', bbox_to_anchor=(0.9, 1.10),
			  fancybox=True, shadow=True, ncol=2)  
			  
	return L,eta_HKTrange

#---------------------------------------------------------------------------------------	
#---------------------------------------------------------------------------------------
# Function which writes Xfoil run file, runs Xfoil and obtains Cl, Cd and L/D values
# Inputs: AoA values, number of AoA values, REYN and airfoil data filename
# Outputs: Cl, Cd and L_over_D values corresponding to AoA
def xfoilrun(alpha,alphacount,REYN,airfoil):

	airfoildatafile = str(airfoil) + '.dat'
	airfoilname     = str(airfoil)
	f = open(airfoilname + ".run", "w")
	f.write("PLOP" + '\n')
	f.write("G" + '\n')
	f.write("C" + '\n')
	f.write('\n')
	f.write("LOAD" + '\n')
	f.write(str(airfoil) + ".dat" + '\n')
	f.write("OPER" + '\n')
	f.write("VISC" + '\n')
	f.write(str(int(REYN)) + '\n')
	f.write("PACC" + '\n')
	f.write(str(airfoil) + "_Re_" + str(int(REYN)) + ".dat" + '\n')
	f.write("x" + str(airfoil) + "_Re_" + str(int(REYN)) + ".dat" + '\n')
	#for i in range(alphacount):
	#	f.write("ALFA" + " " + ('%2.4f'%alpha[i]) + '\n')
	f.write("ASEQ" + '\n')
	f.write("-8" + '\n')
	f.write("8" + '\n')
	f.write("0.25")
	f.write('\n')
	f.write('\n')

	f.write("QUIT" + '\n')
	f.write("EOF" + '\n')
	f.write('\n')	
	# f.write('\n')	


	f.close
	##########################
	# Creating xfoilrun bash script for linux machines
	f = open("xfoilrun", "w")
	f.write("xfoil << EOF" + '\n')
	f.write("PLOP" + '\n')
	f.write("G" + '\n')
	f.write("C" + '\n')
	f.write('\n')
	f.write("LOAD" + '\n')
	f.write(str(airfoil) + ".dat" + '\n')
	f.write("OPER" + '\n')
	f.write("VISC" + '\n')
	f.write(str(int(REYN)) + '\n')
	f.write("PACC" + '\n')
	f.write(str(airfoil) + "_Re_" + str(int(REYN)) + ".dat" + '\n')
	f.write("x" + str(airfoil) + "_Re_" + str(int(REYN)) + ".dat" + '\n')
	for i in range(alphacount):
		f.write("ALFA" + " " + ('%2.4f'%alpha[i]) + '\n')
	f.write('\n')
	f.write('\n')
	f.write("QUIT" + '\n')	
	f.write("EOF" + '\n')
	f.write('\n')
	f.close
	##########################
	# Running Xfoil and creating the polar data
	polarfile = str(airfoil) + "_Re_" + str(int(REYN)) + ".dat"
	xpolarfile = "x" + str(airfoil) + "_Re_" + str(int(REYN)) + ".dat"
	# Removing the files if they exist before the current run	
	# if os.path.isfile(polarfile):
		# os.remove(polarfile)
	# if os.path.isfile(xpolarfile):
		# os.remove(xpolarfile)
	###########################	
	workdir = os.getcwd()
	print
	print 'Running Xfoil to obtain Cl, Cd and L/D for AoA values...'
	print
	xfoilrunfile = airfoilname + ".run"
	# xfoilrun     = "xfoil <" + xfoilrunfile + " >xfoilrun.log" 
	#############
	#xfoilrun = "xfoil < e857.run" #+ xfoilrunfile # + ' >xfoilrun.log'
	#subprocess.Popen('chmod u+x xfoilrun', shell=True)	
	#subprocess.Popen('./xfoilrun', shell = True)
	# subprocess.Popen(xfoilrun, shell = True)
	# Reading the Xfoil polar output file
	time.sleep(0.25)
	xfoilpolar = str(airfoil) + "_Re_" + str(int(REYN)) + ".dat"
	file_length = file_len(xfoilpolar)
	f = open(xfoilpolar, "r")
	lineread = f.readlines()
	# lines = f.read()#.split()
	#for i in range(file_length):
		#print lineread[i]
	f.close()
	print "Number of lines: ", file_length
	print	
	# Calculating how many successful polars obtained
	alpha_polars = file_length - int(12) #116# 12 lines before value list
	print
	starter = "alpha"
	alphaline = 0
	for numm, line in enumerate(lineread, 1):
		if starter in line:
			alphaline = numm - 1
			print 'alphaline', alphaline
	numm2 = alphaline + 2    
	print "alpha_polars:", alpha_polars
	alfapolars = np.ndarray(alpha_polars)
	Clpolars   = np.ndarray(alpha_polars)
	Cdpolars   = np.ndarray(alpha_polars)
	Clspan     = np.ndarray(alphacount)
	Cdspan     = np.ndarray(alphacount)	
	L_over_D   = np.ndarray(alphacount)
	## Interpolating the polars for spanwise alpha values
	for i in range(alpha_polars):
		alfapolars[i] = float(lineread[numm2+i].split()[0])
		Clpolars[i]   = float(lineread[numm2+i].split()[1])
		Cdpolars[i]   = float(lineread[numm2+i].split()[2])  
		# print '%.4f'%alfapolars[i], '%.4f'%Clpolars[i], '%.4f'%Cdpolars[i]
	######
	### PLotting the alpha vs Cl and Cd
	# fignum = 1 
			 
	# py.figure(fignum, figsize=(10, 6))    
	# py.plot(alfapolars,Clpolars,color = 'red')
	# py.plot(alfapolars,Cdpolars)
	# py.title("Cl and Cd vs AoA")
	# py.ylabel("Cl and Cd")
	# py.xlabel("AoA")
	# # py.axis([0.0,1.0,0.0,1.4])
	# # py.axis('equal')
	# leg = py.legend(['Cl','Cd'], loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0.)
	# py.savefig('Xfoil_polars.png', bbox_extra_artists=(leg,), bbox_inches='tight')	

	### Interpolating the values to get Cl, Cd and L/D for
	### a specific alpha (spanwise)
	Clspan = subroutines.curv_line_inters(alpha,alphacount,Clpolars,alfapolars,alpha_polars)
	Cdspan = subroutines.curv_line_inters(alpha,alphacount,Cdpolars,alfapolars,alpha_polars)
	print
	print "alpha    Cl      Cd      L/D"	
	for i in range(alphacount):
		L_over_D[i] = Clspan[i]/Cdspan[i]
		print '%.5f'%alpha[i], '%.5f'%Clspan[i], '%.5f'%Cdspan[i], '%.5f'%L_over_D[i] 
	
	return Clspan, Cdspan, L_over_D
#---------------------------------------------------------------------------------------	
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------	
#---------------------------------------------------------------------------------------
# Function which writes Xfoil run file, runs Xfoil and obtains Cl, Cd and L/D values
# Inputs: AoA values, number of AoA values, REYN and airfoil data filename
# Outputs: Cl, Cd and L_over_D values corresponding to AoA
def xfoilrun_REYN(REYN,airfoil,alpha,propeller):

	airfoildatafile = str(airfoil) +str(int(REYN))+ '.dat'
	airfoilname     = str(airfoil) + str(int(REYN))
	f = open(airfoilname + ".run", "w")
	f.write("PLOP" + '\n')
	f.write("G" + '\n')
	f.write("C" + '\n')
	f.write('\n')
	f.write("LOAD" + '\n')
	f.write(str(airfoil) + ".dat" + '\n')
	f.write("OPER" + '\n')
	f.write("VISC" + '\n')
	f.write(str(int(REYN)) + '\n')
	f.write("PACC" + '\n')
	f.write(str(airfoil) + "_Re_" + str(int(REYN)) + ".dat" + '\n')
	f.write('\n')
	f.write("ASEQ" + '\n')
	if propeller:
		f.write("-7" + '\n')
		f.write("0" + '\n')
	else: # for turbines 
		f.write("0" + '\n')
		f.write("9" + '\n')    
	f.write("0.25" + '\n')
	# f.write("ALFA " + str(float(alpha)) + '\n')
	f.write('\n')
	f.write("QUIT" + '\n')
	f.write("EOF")
	f.close
	##########################
	##########################
	# # Running Xfoil and creating the polar data
	# polarfile = str(airfoil) + "_Re_" + str(int(REYN)) + ".dat"
	xpolarfile = "x" + str(airfoil) + "_Re_" + str(int(REYN)) + ".dat"
	# # Removing the files if they exist before the current run	
	# # if os.path.isfile(polarfile):
		# # os.remove(polarfile)
	# ###########################	
	# workdir = os.getcwd()
	# print
	# print 'Running Xfoil to obtain Cl, Cd and L/D for AoA values...'
	# print
	# xfoilrunfile = airfoilname + ".run"
	# # creating a batchfile
	# xfoilbatch = "batchxfoil.bat"
	# f = open(xfoilbatch,"w")
	# f.write("C:\\Users\\kolimanja\\Dropbox\\Hydroturbine\\GitHydroBEMT\\XFoil_Data\\xfoil<"+str(airfoil)+str(int(REYN))+".run>xfoilrun.log")
	# f.close
	# #############
	# # p = Popen("batchxfoil.bat", cwd=r"C:\Users\kolimanja\Dropbox\Hydroturbine\GitHydroBEMT\Propeller\UC\soloProp\1+a\OptimizationBEMT\Run01\optimum_run\temp")
	# # stdout, stderr = p.communicate()
	# dir="C:\\Users\\kolimanja\\Dropbox\\Hydroturbine\\GitHydroBEMT\\Propeller\\UC\\soloProp\\1+a\\OptimizationBEMT\\Run01\\optimum_run\\temp\\"
	# # subprocess.call([dir+"batchxfoil.bat"])
	# # os.system('xfoil <'+xfoilrunfile+ '> xfoil.'+airfoil+'.'+REYN+'.log')
	# process = Popen("batchxfoil.bat",shell=True)
	# #time.sleep(3)	
	# # os.remove(xfoilbatch)	
	# # subprocess.call(['taskkill', '/F', '/T', '/PID', str(p.pid)])	
	# # os.system("pause")    
	# xfoilpolar = str(airfoil) + "_Re_" + str(int(REYN)) + ".dat"
	# file_length = file_len(xfoilpolar)
	# f = open(xfoilpolar, "r")
	# lineread = f.readlines()
	# # lines = f.read()#.split()
	# #for i in range(file_length):
		# #print lineread[i]
	# f.close()
	# # me = os.getpid()
	# # print me
	# # os.system("pause")
	# # kill_proc_tree(me)
	if os.path.isfile(xpolarfile):
		os.remove(xpolarfile)    
	# print "Number of lines: ", file_length
	# print	
	# # Calculating how many successful polars obtained
	# alpha_polars = file_length - int(12) #116# 12 lines before value list
	# print
	# starter = "alpha"
	# alphaline = 0
	# for numm, line in enumerate(lineread, 1):
		# if starter in line:
			# alphaline = numm - 1
			# print 'alphaline', alphaline
	# numm2 = alphaline + 2    
	# print "alpha_polars:", alpha_polars
	# alfapolars = np.ndarray(alpha_polars)
	# Clpolars   = np.ndarray(alpha_polars)
	# Cdpolars   = np.ndarray(alpha_polars)
	# ## Interpolating the polars for spanwise alpha values
	# for i in range(alpha_polars):
		# alfapolars[i] = float(lineread[numm2+i].split()[0])
		# Clpolars[i]   = float(lineread[numm2+i].split()[1])
		# Cdpolars[i]   = float(lineread[numm2+i].split()[2])  
		# # print '%.4f'%alfapolars[i], '%.4f'%Clpolars[i], '%.4f'%Cdpolars[i]
	# ######
	# ### PLotting the alpha vs Cl and Cd
	# fignum = 1 
			 
	# py.figure(fignum, figsize=(10, 6))    
	# py.plot(alfapolars,Clpolars,color = 'red')
	# py.plot(alfapolars,Cdpolars)
	# py.title("Cl and Cd vs AoA"+" for REYN "+str(int(REYN)))
	# py.ylabel("Cl and Cd")
	# py.xlabel("AoA")
	# py.axis([-8.0,0.0,-0.8,0.8])
	# # py.axis('equal')
	# leg = py.legend(['Cl','Cd'], loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0.)
	# py.savefig(('Xfoil_polars.'+str(int(REYN))+'.png'), bbox_extra_artists=(leg,), bbox_inches='tight')	
	# Popen("TASKKILL /F /PID {pid} /T".format(pid=process.pid))

	
	return #Clspan, Cdspan, L_over_D
#---------------------------------------------------------------------------------------	
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------	
#---------------------------------------------------------------------------------------
# Writes output for analysis mode and contains
# 1) Airfoil properties, Vz, RPM, dimensions, nsections, Design Thrust, Torque and Power
# 2) spanwise properties of the design point: 
#     span, radius, a, a', chord (m), stagger, phi, Cl, Cd, alpha, Ca, Cr, (dT/dr)_BMT, (dT/dr)_BET, (dQ/dr)_BMT, (dQ/dr)_BET,
#
# Can vary Vz or RPM or both
# OUTPUTS:
#    spanwise dT/dr, dQ/dr, Velocities at various stations
#    Thrust, Torque, Power
#
def write_analysis_input(arg1,arg2,case):
	Rhub,Rtip,Thrust,Torque,Power,RPM,Vz,Prop_eta = arg1
	radius, axialIF, anglrIF, chordr, REYN, stagger, phi, Clarray, Cdarray, alpha, dTdrBMT, dTdrBET, dQdrBMT, dQdrBET,V2,V3,Vexit = arg2
	
	alpha_array = np.zeros(nsect)
	for i in range(nsect):
		alpha_array[i] = alpha
	print " Writing analysis input file for this design..."
	# open the file to write to
	f = open("BEMTanalysis" + "." + str(case) + ".dat", "w")
	f.write("#####################################################" + '\n') 
	f.write("# BEMT Blade ANALYSIS Input Parameters [S.I. Units] #" + '\n')    
	f.write("#####################################################" + '\n')   
	f.write(str(device) + "    # Device type" + '\n')
	f.write(str(case)   + "    # Casename # Solo Propellers" + '\n')    
	f.write('\n')    
	f.write("##########################################" + '\n')    
	f.write("# Fluid Properties"+ '\n')      
	f.write("##########################################" + '\n') 
	f.write(str(fluid)   + "    # Fluid Type" + '\n')
	f.write(str(rho)     + "    # Density [kg/m^3]" + '\n')   
	f.write(str(nju)     + "    # Kinematic Viscosity @ 20C [m^2/s]" + '\n')
	f.write(str(a_sound) + "    # Speed of Sound [m/s^2]" + '\n')    
	f.write('\n')
	f.write("##########################################" + '\n') 
	f.write("# Blade Parameters" + '\n')  
	f.write("##########################################" + '\n')     
	f.write('%2.3f'%(Vz1[1]) + "   #Vz" + '\n') 
	f.write('%2.3f'%(RPM)    + "   # RPM" + '\n') 
	f.write('%2.3f'%(Rhub)   + "   # Hub Radius [m]" + '\n') 
	f.write('%2.3f'%(Rtip)   + "   # Tip Radius [m]" + '\n') 
	f.write('%02d'%(nbrow)  + "   # Number of blade rows" + '\n') 
	f.write('%02d'%(nblade) + "   # Number of blades" + '\n') 
	f.write('%02d'%(nsect)  + "   # Number of airfoil sections" + '\n')  
	f.write('\n')    
	f.write("##########################################" + '\n') 
	f.write( "#Analysis_variables" + '\n')
	f.write("##########################################" + '\n')        
	f.write( "#parameter LL UL step" + '\n')  
	f.write( "Vz 3 10 10" + '\n') 
	f.write('\n')
	f.write("##########################################" + '\n') 
	f.write("# Airfoil Parameters" + '\n')  
	f.write("##########################################" + '\n')  
	f.write(str(foilname) + "   # Type of Airfoil" + '\n') 
	f.write('\n') 
	f.write("##########################################" + '\n') 
	f.write("#Design Point properties" + '\n')  
	f.write("##########################################" + '\n')      
	f.write(('%4.4f'%Thrust) + " # Thrust (N)" + '\n')
	f.write(('%4.4f'%Torque) + " # Torque (Nm)" + '\n')
	f.write(('%4.4f'%Power) + "  # Power (W) " + '\n')
	f.write(('%2.4f'%Prop_eta) + " # Propeller Efficiency " + '\n')   
	f.write('\n')
	f.write("radius  	chord[m]     Stagger	 Phi    REYN      AoA     Cl      Cd   L/D     a    a'	 dThrust/dr[N]		dTorque/dr[N/m]		V2[m/s]		V3[m/s]		Vexit[m/s]	" + '\n')
	for i in range(nsect):
		f.write(('%3.8f'%radius[i]) + "  ")
		f.write(('%3.8f'%chordr[i] ) + "  ")
		f.write(('%3.8f'%stagger[i] ) + "  ")
		f.write(('%3.8f'%phi[i] ) + "  ")
		f.write(('%3.8f'%REYN[i]  ) + "  ")
		f.write(('%3.8f'%alpha_array[i] ) + "  ")        
		f.write(('%3.8f'%Clarray[i] ) + "  ")
		f.write(('%3.8f'%Cdarray[i] ) + "  ")
		f.write(('%3.8f'%(Clarray[i]/Cdarray[i]) ) + "  ")
		f.write(('%3.8f'%axialIF[i]  ) + "  ")  
		f.write(('%3.8f'%anglrIF[i]  ) + "  ")  		
		f.write(('%3.8f'%dTdrBET[i]  ) + "  ") 
		f.write(('%3.8f'%dQdrBET[i]  ) + "  ")       
		f.write(('%3.8f'%V2[i]  ) + "  ")   
		f.write(('%3.8f'%V3[i]  ) + "  ")   
		f.write(('%3.8f'%Vexit[i]  ) + "  " + '\n')   
	f.close()        
#---------------------------------------------------------------------------------------	
#---------------------------------------------------------------------------------------