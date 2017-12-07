import sys
import math
import numpy as np
import pylab as py
from scipy.optimize import brentq
from matplotlib.pyplot import *
from functions import * 
import subroutines
import numpy as np
np.seterr(divide='ignore', invalid='ignore')


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

def calcChord(rsv,Phiradsv,nblade,Cl,lambda_rsv):

	# Calculates chord. (single value)
	chordrsv  = (8*(math.pi)*rsv*math.sin(Phiradsv))/(3*nblade*Cl*lambda_rsv) # optimal chord (obtained for max Power)

	return chordrsv
#---------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# hub_loss_factor
def hublossfactor(nblade,lambda_rsingle,Phiradsingle,rsingle,Rhub):

	term         = 2./math.pi
	# term1        = -0.125*((nblade*lambda_rsingle) - 21)
	# b1           = math.exp(term1) + 0.1
	b            = 1.
	
	#### hub loss factor
	term2        = (nblade/2.)*((rsingle - Rhub)/(rsingle*abs(math.sin(Phiradsingle))))
	termA        = np.exp(-b*term2)       
	fhubloss     = term*math.acos(termA)	
	
	return fhubloss
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# tip_loss_factor
def tiplossfactor(nblade,lambda_rsingle,Phiradsingle,rsingle,Rtip):

	term         = 2./math.pi
	# term1        = -0.125*((nblade*lambda_rsingle) - 21)
	# b1           = math.exp(term1) + 0.1
	b            = 1.
	
	#### tip loss factor
	term2        = (nblade/2.)*((Rtip - rsingle)/(rsingle*abs(math.sin(Phiradsingle))))
	# print "Phirad,sinPhirad,term2:",Phiradsingle,math.sin(Phiradsingle),term2
	termA        = np.exp(-b*term2)     
	# print "termA:",termA
	ftiploss     = term*math.acos(termA)
	#corrected ftiploss near the tip
	#ftiploss     = term*math.acos(math.exp(-b1*term3))
	# ftiploss[i]     = 1.	
	
	return ftiploss
#---------------------------------------------------------------------------------------
def total_loss(nblade,lambda_rsv,Phiradsv,rsv,Rhub,Rtip):

	#--------------------------------------------------------------   
	# Prandtl's hub and tip loss factors
	#--------------------------------------------------------------        
	#### hub loss factor
	# fhublosssv   = hublossfactor(nblade,lambda_rsv,Phiradsv,rsv,Rhub)
	fhublosssv     = 1.
	#        
	#### tip loss factor        
	ftiplosssv     = tiplossfactor(nblade,lambda_rsv,Phiradsv,rsv,Rtip)
	#      
	#### total loss factor
	flosssv        = fhublosssv*ftiplosssv 
	
	return flosssv
	
#---------------------------------------------------------------------------------------    
#---------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

def Ca_Cr(Cl,Cd,Phiradsv,propeller):

	# Using the same vector triangle for both Propellers and Turbines
	if propeller:
		Casv       =  Cl*math.cos(Phiradsv) - Cd*math.sin(Phiradsv) # also called Ct
		Crsv       =  Cl*math.sin(Phiradsv) + Cd*math.cos(Phiradsv) # also called Cn
	else: # turbine
		Casv       =  Cl*math.cos(Phiradsv) + Cd*math.sin(Phiradsv) # also called Ct
		Crsv       =  Cl*math.sin(Phiradsv) - Cd*math.cos(Phiradsv) # also called Cn
		
	# print "Ca,Cr: ",Casv,Crsv

	return Casv,Crsv
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
def calcAxlIFnew(currentrow,co_rotation,counter_rotation,propeller,Casv,Crsv,Phiradsv,flosssv,sigmarsv,gap):

	# Calculating the new axlIF. Single value.
	sphi  = math.sin(Phiradsv)
	kappa = (4.*flosssv*(sphi)**2)/(sigmarsv*Casv)
	#
	if propeller:
		if counter_rotation or co_rotation:
			if currentrow == 1: #front rotor   
				#gap = 1.0
				print "gap in rotor 1:",gap
				A = 1 - kappa*gap
				B = kappa*gap - 1
				C = (kappa*gap**2) + 1
				axlIFnewsv = (A - math.sqrt((B**2) - (C)))/C
			elif currentrow == 2: # aft rotor
				term4sv    = ((4*flosssv*(math.sin(Phiradsv))**2)/(sigmarsv*Casv)) - 1.     
				axlIFnewsv = 1./(term4sv)            			
		else: # single row 
			term4sv    = ((4*flosssv*(math.sin(Phiradsv))**2)/(sigmarsv*Casv)) - 1.     
			axlIFnewsv = 1./(term4sv)            
	else: # turbine
		if counter_rotation or co_rotation:
			if currentrow == 1: # front rotor
				A = 1 + kappa*gap
				B = kappa*gap + 1                
				C = (kappa*gap**2) + 1
				axlIFnewsv = (A - math.sqrt((B**2) - (C)))/C
			elif currentrow == 2: # aft rotor
				term4sv    = ((4*flosssv*(math.sin(Phiradsv))**2)/(sigmarsv*Casv)) + 1.
				axlIFnewsv = 1./(term4sv)			
		else: # single row 
			term4sv    = ((4*flosssv*(math.sin(Phiradsv))**2)/(sigmarsv*Casv)) + 1.
			axlIFnewsv = 1./(term4sv)  
	
	return axlIFnewsv

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# function to calculate axial induction factor if a > 0.4 using Glauert's correction
def calc_axlIF_Glauert(propeller,axlIF,flosssv,phiradsv,sigmarsv,Casv,CThrust):

	#
	# root_term = math.sqrt((CThrust*(50 - (36*flosssv))) + ((12*flosssv)*(3*flosssv - 4)))
	# Nrr = (18*flosssv) - 20 - 3*(root_term)
	# Drr = 36*flosssv - 50
	# axlIFGlrt = Nrr/Drr
	
	# Nrr2 = 0.143 + math.sqrt(0.0203 - 0.60427*(0.889 - CThrust))
	# axlIFGlrt = Nrr2/flosssv
	#
	#if flosssv == 0.0:
	#	flosssv = 0.1
	sphi  = math.sin(phiradsv)
	if flosssv == 0.0:
		flosssv = 0.0001	
	kappa = (sigmarsv*Casv)/(4.*flosssv*(sphi)**2)
	# print "flosssv: ",flosssv
	# print "kappa: ", kappa
	N1   = 2.*flosssv*kappa - ((10./9.) - flosssv)
	A    = (2.*flosssv*kappa)
	B    = (flosssv*((4./3.) - flosssv))
	# if (A-B) < 0. :
		# print 'A,B:',A,B
	N2   = math.sqrt(abs(A - B))
	D1   = 2.*flosssv*kappa - ((25./9.) - (2.*flosssv))
	######
	# print "D1 :",D1    
	if D1 == 0. :
		D1 = 0.000001
	axlIFGlrt = (N1 - N2)/D1 
	
	return axlIFGlrt

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
def calcAngIFnew(propeller,Casv,Crsv,Phiradsv,flosssv,sigmarsv):

	# Calculates angIFnew (single value)
	Nr2sv      = 4.0*flosssv*math.sin(Phiradsv)*math.cos(Phiradsv)
	Dr2sv      = sigmarsv*Crsv 
	if propeller:
		term5sv    = ((Nr2sv)/(Dr2sv)) + 1.    
	else:
		term5sv    = ((Nr2sv)/(Dr2sv)) - 1.
	angIFnewsv = 1./(term5sv)	
		
	return angIFnewsv
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
def calcAngIF12(axlIF,angIF,lambda_rsv,propeller):
	# Calculates the angIF12 caused by the wake rotation of rotor 1 a the inlet of aft rotor

	t1  =  axlIF*(1 + axlIF)  
	t2  =  angIF*(1 + angIF)
	t3  =  -4.*((t1/(lambda_rsv)**2) + t2)
	t4  =  1. - (t3)
	if t4 < 0.:
		print "Imaginary roots found for angIFnew12."
		t4 = 0.
	#
	angIF12 = 0.5*(-1. + math.sqrt(t4)) # taking only the (-,+) root.   

	return angIF12
	
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Fixed point iteration
def fixedPointIteration(currentrow,co_rotation,counter_rotation,propeller,file_length,chordrsv,chrdmpliersv,lambda_rsv,nbrow,rsv,Rhub,Rtip,nblade,Cl,Cd,gap):
	deftol = 1.0e-5
	itr    = 100

	axlIF = 0.0
	angIF = 0.0
	#
	error1 = 1.0e-1
	error2 = 1.0e-1
	itr = 0
	#
	while (error1 > deftol and error2 > deftol):
		Nr = 1 - axlIF
		Dr = (1 + angIF)*lambda_rsv
		if Dr == 0.0:
		 Dr = 0.0001*lambda_rsv
		#
		phirad = math.atan2(Nr,Dr)
		if phirad == 0.0:
		 print " Phi is zero... "
		 print " New Phi is 0.00001"
		 phirad = 0.00001
		 print
		#
		if internalchord and not definedchord:
			chordrsv  = calcChord(rsv,phirad,nblade,Cl,lambda_rsv)# optimal chord (obtained for max Power)
			chordrsv  = math.fabs(chordrsv)
			
		# Using the chord multiplier to the chord
		chordrsv   = chordrsv*chrdmpliersv
		sigmarsv   = (nblade*chordrsv)/(2.*math.pi*rsv)
		Cafinal,Crfinal = Ca_Cr(Cl,Cd,phirad,propeller)      
		flossfinal      = total_loss(nblade,lambda_rsv,phirad,rsv,Rhub,Rtip)    
		#
		axlIFnew  = calcAxlIFnew(currentrow,co_rotation,counter_rotation,propeller,Cafinal,Crfinal,phirad,flossfinal,sigmarsv,gap)
		angIFnew  = calcAngIFnew(propeller,Cafinal,Crfinal,phirad,flossfinal,sigmarsv)
		#
		# Error calculation
		error1 = math.fabs(axlIFnew - axlIF)
		error2 = math.fabs(angIFnew - angIF)
		#
		# Current solution as initial solution for next iteration
		axlIF = axlIFnew
		angIF = angIFnew
		
		# next iteration count
		itr = itr + 1
		#
		print
		print "Iterations:",itr
		print        

	return phirad,axlIF,angIF,flossfinal,Cafinal,Crfinal,sigmarsv,chordrsv
	
#---------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------- 
# calculate a and a'

def calc_a_aprime(currentrow,co_rotation,counter_rotation,propeller,file_length,chordrsv,chrdmpliersv,lambda_rsv,nbrow,rsv,Rhub,Rtip,nblade,Cl,Cd,gap):
	#
	print "Cl,Cd in the frontrotor a' routine :", Cl, Cd
	# Find phi using Brent's method
	phiradstar = findphi(propeller,chordrsv,lambda_rsv,nblade,rsv,Rhub,Rtip,Cl,Cd,currentrow,co_rotation,counter_rotation,gap)
	#
	# Final values
	# print 'internalchord:',internalchord
	if internalchord and not definedchord:
		chordrsv  = calcChord(rsv,phiradstar,nblade,Cl,lambda_rsv)# optimal chord (obtained for max Power)
		chordrsv  = math.fabs(chordrsv)
		# print 'chordrsv:',chordrsv
	#
	# print chordrsv
	# Using the chord multiplier to the chord
	# print "Chord multiplier: ",chrdmpliersv
	chordrsv   = chordrsv*chrdmpliersv
	print "New chord: ",chordrsv
	sigmarsv   = (nblade*chordrsv)/(2.*math.pi*rsv)
	Cafinal,Crfinal = Ca_Cr(Cl,Cd,phiradstar,propeller)      
	flossfinal      = total_loss(nblade,lambda_rsv,phiradstar,rsv,Rhub,Rtip)    
	#
	axlIF  = calcAxlIFnew(currentrow,co_rotation,counter_rotation,propeller,Cafinal,Crfinal,phiradstar,flossfinal,sigmarsv,gap)
	# Calculating the coefficient of thrust
	CThrust = 4.*flossfinal*axlIF*(1 - axlIF)    
	#########################################
	# Use this correction for the front rotor only in a Counter Rotation Device
	# Check if axlIF is greater than 0.45 and apply Gauert's correction
	if not propeller:
		if math.fabs(axlIF) > 0.45:
			print "Applying Glauert and Spera correction..."
			print
			print "abolute value Old AxlIF = ",math.fabs(axlIF), " > 0.45 "
			axlIFglrt = calc_axlIF_Glauert(propeller,axlIF,flossfinal,phiradstar,sigmarsv,Cafinal,CThrust)
			print "New AXlIF = ",axlIFglrt  
			axlIF = axlIFglrt    
	###############################################        
	angIF  = calcAngIFnew(propeller,Cafinal,Crfinal,phiradstar,flossfinal,sigmarsv)

	return phiradstar,axlIF,angIF,flossfinal,Cafinal,Crfinal,sigmarsv,chordrsv

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Brent's method to find phi from fzero
def findphi(propeller,chordrsv,lambda_rsv,nblade,rsv,Rhub,Rtip,Cl,Cd,currentrow,co_rotation,counter_rotation,gap):

	eps = 1e-6
	pi  = math.pi
	#
	phi_lwr = eps
	phi_upr = pi/2.    
	#
	# arguments based on the index for fphi
	args = (lambda_rsv,chordrsv,propeller,Cl,Cd,rsv,Rhub,Rtip,nblade,currentrow,co_rotation,counter_rotation,gap)
	#
	fzero_lwr = fphi(phi_lwr,*args)
	fzero_upr = fphi(phi_upr,*args)
	#
	if fzero_lwr*fzero_upr > 0:
		print "fzero_lwr*fzero_upr > 0"
		#
		fzero_1 = fphi((3.*pi/4.),*args)
		#
		fzero_2 = fphi(eps,*args)
		if fzero_1 < 0 and fzero_2 > 0:
			phi_lwr = (3.*pi/4.)
			phi_upr = eps
		else:
			if propeller:
				phi_lwr = eps
				phi_upr = pi/2.
			else:
				phi_lwr = pi/2.
				phi_upr = pi - eps            
	#       
	# try:
	radphistar = brentq(fphi,phi_lwr,phi_upr,args = args)
	# print radphistar
	# except:
		# radphistar = 0.0
		# print "Phi is zero. Check your input values."
		
	return radphistar

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# function which calculates fzero for a 1D equation of phi 
def fphi(inphirad,lambda_rsv,chordrsv,propeller,Cl,Cd,rsv,Rhub,Rtip,nblade,currentrow,co_rotation,counter_rotation,gap):
	# calculates the equation f(phi)
	# f_prop(phi) = [sin(phi)/(1 + axlIFnewsv(phi))] - [cos(phi)/((1 - angIFnewsv(phi))*lambda_rsv)]
	# f_WT(phi)   = [sin(phi)/(1 - axlIFnewsv(phi))] - [cos(phi)/((1 + angIFnewsv(phi))*lambda_rsv)]

	#
	sphi = math.sin(inphirad)
	cphi = math.cos(inphirad)
	#
	if internalchord and not definedchord:
		chordrsv  = calcChord(rsv,inphirad,nblade,Cl,lambda_rsv)# optimal chord (obtained for max Power)
	#
	sigmarsv   = (nblade*chordrsv)/(2.*math.pi*rsv)
	Casv,Crsv  = Ca_Cr(Cl,Cd,inphirad,propeller)      
	flossfinal = total_loss(nblade,lambda_rsv,inphirad,rsv,Rhub,Rtip)
	axlIFsv = calcAxlIFnew(currentrow,co_rotation,counter_rotation,propeller,Casv,Crsv,inphirad,flossfinal,sigmarsv,gap)  
	angIFsv = calcAngIFnew(propeller,Casv,Crsv,inphirad,flossfinal,sigmarsv)

	A1 = (1 - axlIFsv)
	A2 = (1 + angIFsv)
	if  A1 == 0.:
		A1 = 0.000001
	if A2 == 0.:
		A2 = 0.000001
	fzero = (sphi/(1 - axlIFsv)) - (cphi/((1 + angIFsv)*lambda_rsv)) 
	
	return fzero

#---------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------
# AFT ROTOR
#---------------------------------------------------------------------------------------
# calculate a and a' for aft rotor

def calc_a_aprime_aftrotor(currentrow,co_rotation,counter_rotation,propeller,file_length,chordrsv,chrdmpliersv,lambda_rsv,nbrow,rsv,Rhub,Rtip,nblade,Cl,Cd,gap):
	#
	print "Cl,Cd in the aftrotor a' routine :", Cl, Cd
	# Find phi using Brent's method
	phiradstar = findphi2(propeller,chordrsv,lambda_rsv,nblade,rsv,Rhub,Rtip,Cl,Cd,currentrow,co_rotation,counter_rotation,gap)
	if propeller:
		phiradstar = (math.pi - phiradstar)
	# Final values
	if internalchord and not definedchord:
		chordrsv  = calcChord(rsv,phiradstar,nblade,Cl,lambda_rsv)# optimal chord (obtained for max Power)
	#
	# Using the chord multiplier to the chord
	# print "Chord multiplier: ",chrdmpliersv
	chordrsv   = chordrsv*chrdmpliersv
	# print "New chord: ",chordrsv
	sigmarsv   = (nblade*chordrsv)/(2.*math.pi*rsv)
	Cafinal,Crfinal = Ca_Cr(Cl,Cd,phiradstar,propeller)      
	flossfinal      = total_loss(nblade,lambda_rsv,phiradstar,rsv,Rhub,Rtip)    
	#   
	axlIF2  = calcAxlIFnew(currentrow,co_rotation,counter_rotation,propeller,Cafinal,Crfinal,phiradstar,flossfinal,sigmarsv,gap)
	# Calculating the coefficient of thrust
	CThrust = 4.*flossfinal*axlIF2*(1 - axlIF2)    
	# # # Check if axlIF2 is greater than 0.45 and apply Gauert's correction
	if not propeller:
		if math.fabs(axlIF2) > 0.45:
			 print "Applying Glauert and Spera correction..."
			 print
			 print "abolute value Old AxlIF = ",math.fabs(axlIF2), " > 0.45 "
			 axlIFglrt = calc_axlIF_Glauert(propeller,axlIF2,flossfinal,phiradstar,sigmarsv,Cafinal,CThrust)
			 print "New AXlIF = ",axlIFglrt 
			 axlIF2 = axlIFglrt    
	######
	sphi = math.sin(phiradstar)
	cphi = math.cos(phiradstar)
	kappa1 = (4.*flossfinal*sphi**2)/(sigmarsv*Cafinal)
	kappa2 = (4.*flossfinal*sphi*cphi)/(sigmarsv*Crfinal)
	k_invr = 1./(kappa2 - 1)
	# angIF2 = 0.5*(k_invr + ((axlIF2*(1 + axlIF2))/(k_invr*lambda_rsv**2)) - 1)
	angIF2 = calcAngIFnew(propeller,Cafinal,Crfinal,phiradstar,flossfinal,sigmarsv)    
		
	angIF12 = calcAngIF12(axlIF2,angIF2,lambda_rsv,propeller)  
	
	# if propeller:
		# axlIF2 = -axlIF2
		# angIF2 = -angIF2
		# angIF12 = -angIF12
		
	return phiradstar,axlIF2,angIF2,angIF12,flossfinal,Cafinal,Crfinal,sigmarsv,chordrsv
#---------------------------------------------------------------------------------------   

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Brent's method to find phi from fzero
def findphi2(propeller,chordrsv,lambda_rsv,nblade,rsv,Rhub,Rtip,Cl,Cd,currentrow,co_rotation,counter_rotation,gap):

	eps = 1e-6
	pi  = math.pi
	#
	phi_lwr = eps
	phi_upr = pi    
	#
	# arguments based on the index for fphi2
	args = (lambda_rsv,chordrsv,propeller,Cl,Cd,rsv,Rhub,Rtip,nblade,currentrow,co_rotation,counter_rotation,gap)
	#
	fzero_lwr = fphi2(phi_lwr,*args)
	fzero_upr = fphi2(phi_upr,*args)
	#
	if fzero_lwr*fzero_upr > 0:
		print "fzero_lwr*fzero_upr > 0"
		#
		fzero_1 = fphi2((3.*pi/4.),*args)
		#
		fzero_2 = fphi2(eps,*args)
		if fzero_1 < 0 and fzero_2 > 0:
			print "fzero_1 and fzero_2 have opp sign"
			phi_lwr = (3.*pi/4.)
			phi_upr = pi
		else:
			print "other option"
			phi_lwr = eps
			phi_upr = pi
	#       
	# try:
	radphistar = brentq(fphi2,phi_lwr,phi_upr,args = args)
	# print "Phi: ", np.degrees(radphistar)
	# except:
		# radphistar = 0.0
		# print "Phi is zero. Check your input values."
		
	return radphistar
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# function which calculates fzero for a 1D equation of phi2 
def fphi2(inphirad,lambda_rsv,chordrsv,propeller,Cl,Cd,rsv,Rhub,Rtip,nblade,currentrow,co_rotation,counter_rotation,gap):
	# calculates the equation f(phi)
	# f_prop(phi) = [sin(phi)/(1 + axlIFnewsv(phi))] - [cos(phi)/((1 - kinvrs)*lambda_rsv)]
	# f_WT(phi)   = [sin(phi)/(1 - axlIFnewsv(phi))] - [cos(phi)/((1 - kinvrs)*lambda_rsv)]
	# kinvrs   = 1/(kappa2+-1) , + for propeller and - for turbines
	# kappa2    = 4*floss*sinphi2*cosphi2/(sigma2*Crsv)

	sphi = math.sin(inphirad)
	cphi = math.cos(inphirad)
	#
	if internalchord and not definedchord:
		chordrsv  = calcChord(rsv,inphirad,nblade,Cl,lambda_rsv)# optimal chord (obtained for max Power)
		chordrsv  = math.fabs(chordrsv)
	#
	sigmarsv   = (nblade*chordrsv)/(2.*math.pi*rsv)
	Casv,Crsv  = Ca_Cr(Cl,Cd,inphirad,propeller)      
	flossfinal = total_loss(nblade,lambda_rsv,inphirad,rsv,Rhub,Rtip)
	# axlIFnew2 = calcAxlIFnew(currentrow,co_rotation,counter_rotation,propeller,Casv,Crsv,inphirad,flossfinal,sigmarsv,gap)  
	# angIFnew2 = calcAngIFnew(propeller,Casv,Crsv,inphirad,flossfinal,sigmarsv)
	# angIFnew12 = calcAngIF12(axlIFnew2,angIFnew2,lambda_rsv,propeller)
	####
	kappa1 = (4.*flossfinal*sphi**2)/(sigmarsv*Casv)
	kappa2 = (4.*flossfinal*sphi*cphi)/(sigmarsv*Crsv)
	k_invr = 1./(kappa2 - 1)
	# fzero = (sphi*(1 - (1./kappa1))) - (cphi/((1 - k_invr)*lambda_rsv))
	fzero = (sphi*(1 - (1./kappa1))) - (cphi*(1 - (1./kappa2))/lambda_rsv)    
	# fzero = (sphi/(1 - axlIFnew2)) - (cphi/((1 + angIFnew2 - angIFnew12)*lambda_rsv)) 
	
	return fzero

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# Calculating the radial distribution of mass flow rate,Thrust,Torque,Power and Vtheta using BET
def radial_distr_BET(multirow,currentrow,co_rotation,counter_rotation,propeller,rho,Vzsv,sigmarsv,axlIFsv,rsv,Phiradsv,Casv,Crsv,OMEGA,nbrow):

	if propeller:
		k = -1
	else: # turbine
		k = +1
	# k = -1
	#
	# mass flow rate
	dmdotdrsv = rho*(1 - (k*axlIFsv))*Vzsv*2*math.pi*rsv
	#-------------------------------------------------------------------------------
	# Blade Element theory definitions
	#------------------------------------------------------------------------------- 
	# Thrust 
	dThrustBETdrsv = sigmarsv*math.pi*rho*((Vzsv*(1 - (k*axlIFsv))/(math.sin(Phiradsv)))**2)*Casv*rsv
	#print " BET : sigma(1+a)*Ca/(sinphi)^2 :", sigmarsv*(1-axlIFsv)*Casv/(math.sin(Phiradsv))**2
	# Torque
	# print "sigmarsv,Crsv in BET radial fn : ", sigmarsv,Crsv
	dTorqueBETdrsv = sigmarsv*math.pi*rho*((Vzsv*(1 - (k*axlIFsv))/(math.sin(Phiradsv)))**2)*Crsv*(rsv**2)    
	# Power
	# if propeller and currentrow == 2:
		# if dTorqueBETdrsv < 0.:
			# dTorqueBETdrsv = -dTorqueBETdrsv
	dPowerBETdrsv  = (OMEGA*dTorqueBETdrsv) # W/m
	#Tangential_Velocity
	Vthetarsv      = dTorqueBETdrsv/(2.*math.pi*rho*(1. - (k*axlIFsv))*Vzsv*(rsv)**2)

	return dmdotdrsv, dThrustBETdrsv, dTorqueBETdrsv, dPowerBETdrsv, Vthetarsv

#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Calculating the radial distribution of mass flow rate,Thrust,Torque,Power and Vtheta using BMT
def radial_distr_BMT(multirow,currentrow,co_rotation,counter_rotation,propeller,rho,Vzsv,axlIFsv,angIFsv,angIF12sv,rsv,flosssv,OMEGA,nbrow,gap):
	#-------------------------------------------------------------------------------
	# Blade Momentum theory definitions
	#------------------------------------------------------------------------------- 
	if propeller:
		k = -1
	else: # turbine
		k = +1
	# k = -1
	# print 'k is ',k, 'in BMT'
	##
	if multirow:
		Aft = (angIFsv - angIF12sv)
	elif multirow and gap == 1.0: #a'2xOMEGA2 - a'1xOMEGA1
		Aft = angIFsv - (angIF12sv/float(OMEGA))  # Had to use the same variable. Not a good idea but temporary fix.
	elif int(currentrow) == 1:
		Aft = angIFsv
	else: # solo row
		Aft = angIFsv        
	# print "multirow,Aft: ",multirow,Aft
	####
	# print " Counter_rotation,nbrow,currentrow,multirow : ", counter_rotation,nbrow,currentrow,multirow
	if nbrow == 2: #Including the axial gap value for the rotors in a multirow arrangement
		if int(currentrow) == 1: #front rotor
			#gap = 1.0
			# print "gap in rotor 1:",gap
			# Axial Force
			dAForceBMTdrsv = flosssv*math.pi*rho*axlIFsv*gap*(2. - (k*axlIFsv*gap))*(Vzsv**2)*rsv
			# Tangential Force            
			dTForceBMTdrsv = flosssv*math.pi*rho*4*Aft*(1 - (k*axlIFsv))*(Vzsv)*OMEGA*(rsv**3)			
		elif int(currentrow) == 2: # aft rotor, using V4 as Vz here and if gap = 1, V4 = V2
			if gap != 1.0:
				# Axial Force
				dAForceBMTdrsv = flosssv*math.pi*rho*4*axlIFsv*(1 - (k*axlIFsv))*(Vzsv**2)*rsv
				# Tangential Force
				dTForceBMTdrsv = flosssv*math.pi*rho*4*Aft*(1 - (k*axlIFsv))*(Vzsv)*OMEGA*(rsv**3)                
			elif gap == 1.0:
				# Axial Force
				dAForceBMTdrsv = flosssv*math.pi*rho*axlIFsv*gap*(2. - (k*axlIFsv*gap))*(Vzsv**2)*rsv  
				# Tangential Force
				dTForceBMTdrsv = flosssv*math.pi*rho*4*Aft*(1 - (k*axlIFsv))*(Vzsv)*OMEGA*(rsv**3)                 
	else: #single row
		# Axial Force
		# print " single row "
		dAForceBMTdrsv = flosssv*math.pi*rho*4*axlIFsv*(1 - (k*axlIFsv))*(Vzsv**2)*rsv
		# Tangential Force
		dTForceBMTdrsv = flosssv*math.pi*rho*4*angIFsv*(1 - (k*axlIFsv))*(Vzsv)*OMEGA*(rsv**3)          
	####    

	# Power
	dPowerBMTdrsv  = (OMEGA*dTForceBMTdrsv) #W/m   
	

	return dAForceBMTdrsv, dTForceBMTdrsv, dPowerBMTdrsv

#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#Calculating final geometry parameters for front rotor
def geometry_param(propeller,Phiradsv,alpha,Vzsv,axlIFsv,angIFsv,lambda_rsv,currentrow,co_rotation,counter_rotation,gap,angular_term):

	Phisv     = np.degrees(Phiradsv)
	if (currentrow == 2) and counter_rotation:   
		if propeller:
			Wrsv      = Vzsv*math.sqrt((1 + axlIFsv)**2 + (lambda_rsv*(1 - angIFsv))**2) 
			### gap = 1
			if gap == 1.0:
				Wrsv  = math.sqrt((Vzsv*(1 + axlIFsv))**2 + (angular_term)**2) 
			###
			inbetasv  = 90.0 + Phisv            
			pitchsv   = -Phisv - alpha
			stggrsv   = -90.0  + pitchsv           
		else: # turbine
			Wrsv      = Vzsv*math.sqrt((1 - axlIFsv)**2 + (lambda_rsv*(1 + angIFsv))**2)
			###
			inbetasv  = 90.0 - Phisv
			pitchsv   = -Phisv + alpha
			stggrsv   = -90.0 + pitchsv 
	else: #solo
		if propeller:
			Wrsv      = Vzsv*math.sqrt((1 + axlIFsv)**2 + (lambda_rsv*(1 - angIFsv))**2) 
			inbetasv  = 90.0 + Phisv            
			pitchsv   = Phisv + alpha
			stggrsv   = 90.0 - pitchsv
		else: # turbine
			Wrsv      = Vzsv*math.sqrt((1 - axlIFsv)**2 + (lambda_rsv*(1 + angIFsv))**2)
			inbetasv  = 90.0 - Phisv    
			pitchsv   = Phisv - abs(alpha)
			stggrsv   = 90.0 - pitchsv        
	#pitchsv   = Phisv - alpha
	#stggrsv   = 90.0 - pitchsv
   
	return Phisv,inbetasv,pitchsv,stggrsv,Wrsv
#---------------------------------------------------------------------------------------		
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
# defining the print out statement for the radial distribution of blade section properties
def radialPrintProperties(multirow,nsect,span_all,floss_all,axlIF_all,angIF_all,angIF12_all,chrdnd_all,U_all,Wr_all,Vtheta_all,inbeta_all,stggr_all,Phi_all,lambda_r_all):
	print
	print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	print
	# Screen printouts:
	if multirow:
		print "-------------------------------------------------------------------------------"
		print "    span      floss     axlIFnew      angIfnew      angIFnew12 "
		print "-------------------------------------------------------------------------------"
		for i in range(nsect):
			print '%.10f'%span_all[i],'%.10f'%floss_all[i],'%.10f'%axlIF_all[i],'%.10f'%angIF_all[i],'%.10f'%angIF12_all[i]
		#print '%.10f'%span[0],'%.10f'%fhubloss[0],'%.10f'%ftiploss[0],'%.10f'%floss[0],'%.10f'%axlIFnew[0],'%.10f'%angIFnew[0]
		print
		#---------------------------------------------------------
		print
		print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
		print
		print "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
		print "    span   |    chordr_nd     |      U      |      W    |      Vtheta    |    inBeta      |    stagger       |     Phi      |     TSR       |     a2      |     a2'       |     a12'    |"    
		print "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
		print
		# for i in range(nsect-1,2,-1):
		for i in range(nsect):
			print '%.10f'%span[i],"|",'%.10f'%chrdnd_all[i],"|",'%.10f'%U_all[i],"|",'%.10f'%Wr_all[i],"|",'%.10f'%Vtheta_all[i],"|",'%.10f'%inbeta_all[i],"|",'%.10f'%stggr_all[i],"|",'%.10f'%Phi_all[i],"|",'%.10f'%lambda_r_all[i],"|",'%.10f'%axlIF_all[i],"|",'%.10f'%angIF_all[i],"|",'%.10f'%angIF12_all[i],"|"
		#print '%.10f'%span[0],'%.10f'%chrdr_nd[0],'%.10f'%U[0],'%.10f'%Wr[0],'%.10f'%inbeta[0],'%.10f'%stggr[0],'%.10f'%Phi[0],'%.10f'%lambda_r[0]
		print
		print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"        
	else: #single row
		print "-------------------------------------------------------------------------------"
		print "    span      floss     axlIFnew      angIfnew "
		print "-------------------------------------------------------------------------------"    
		#for i in range(nsect-1,2,-1):
		for i in range(nsect):
			print '%.10f'%span_all[i],'%.10f'%floss_all[i],'%.10f'%axlIF_all[i],'%.10f'%angIF_all[i]
		#print '%.10f'%span[0],'%.10f'%fhubloss[0],'%.10f'%ftiploss[0],'%.10f'%floss[0],'%.10f'%axlIFnew[0],'%.10f'%angIFnew[0]
		print
		#---------------------------------------------------------
		print
		print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
		print
		print "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
		print "    span   |    chordr_nd     |      U      |      W    |      Vtheta    |    inBeta      |    stagger       |     Phi      |     TSR       |     a      |     a'       |"    
		print "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
		print
		# for i in range(nsect-1,2,-1):
		for i in range(nsect):
			print '%.10f'%span[i],"|",'%.10f'%chrdnd_all[i],"|",'%.10f'%U_all[i],"|",'%.10f'%Wr_all[i],"|",'%.10f'%Vtheta_all[i],"|",'%.10f'%inbeta_all[i],"|",'%.10f'%stggr_all[i],"|",'%.10f'%Phi_all[i],"|",'%.10f'%lambda_r_all[i],"|",'%.10f'%axlIF_all[i],"|",'%.10f'%angIF_all[i],"|"
		#print '%.10f'%span[0],'%.10f'%chrdr_nd[0],'%.10f'%U[0],'%.10f'%Wr[0],'%.10f'%inbeta[0],'%.10f'%stggr[0],'%.10f'%Phi[0],'%.10f'%lambda_r[0]
		print
		print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	##########################3
	return
#--------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------  
def printFinalValues(Rhubsv,Rtipsv,Vztipsv,Lambdasv,nbladesv,ThrustBETsv,AForceBMTsv,TorqueBETsv,TForceBMTsv,
					 PowerBETsv,PowerBMTsv,Power_swirlsv,Power_fluidsv,C_Psv,CT_formulasv,
					 eta_Propsv,etasv,eta_HKTsv,poundforce,hp,propeller,RPM,currentrow):
					 
	# Defining the print statements to print final values of Thrust, Power and Etas
	#---------------------------------------------------------
	print
	print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	print
	print "Tip Radius, Hub radius : ", Rtipsv, Rhubsv, "m"
	# calcTSR(RPM,Vz,Rtip)
	# calcRPM(Lambdasv,Vztipsv,Rtip)
	print "RPM   : ", RPM
	print "Vztip : ", Vztipsv, " m/s"
	print "TSR   : ", Lambdasv
	print
	print "Total Torque  (BET) for ",nbladesv," rotors : ",TorqueBETsv, "Nm "  
	print
	print "Total Torque  (BMT) for ",nbladesv," rotors : ",TForceBMTsv, "Nm "  
	print
	# print "Total Thrust (BET) each rotor      : ", ThrustBET, "N"  
	print "Total Thrust  (BET) for ",nbladesv," rotors : ",ThrustBETsv, "N [",ThrustBETsv*poundforce, "lbf ]"  
	#---------------------------------------------------------
	#---------------------------------------------------------
	# print "Total AxForce (BMT) each rotor     : ",AForceBMT, "N"  
	print "Total AxForce (BMT) for ",nbladesv," rotors : ",AForceBMTsv, "N [",AForceBMTsv*poundforce, "lbf ]"   
	print
	#---------------------------------------------------------
	#---------------------------------------------------------
	print
	# print "Total Power (BET) each rotor       : ",PowerBET, "W"  
	print "Total Power (BET) for ",nbladesv," rotors   : ",(PowerBETsv), "W  [", (PowerBETsv)/hp, "HP ]"
	#---------------------------------------------------------
	#---------------------------------------------------------
	# print "Total Power (BMT) each rotor       : ",PowerBMT, "W"  
	print "Total Power (BMT) for ",nbladesv," rotors   : ",(PowerBMTsv), "W  [", (PowerBMTsv)/hp, "HP ]"
	print
	print
	print "Total Power_swirl for ",nbladesv," rotors   : ",Power_swirlsv, "W  [", Power_swirlsv/hp, "HP ]"
	print
	#---------------------------------------------------------
	#---------------------------------------------------------
	print
	print "Power for ",Rtipsv," m radius  : ",Power_fluidsv, "W  [", Power_fluidsv/hp, "HP  ]"
	print
	if RPM < 0.:
		print "Row",currentrow,"Coefficient of Power (C_P)    : ",-C_Psv     
	else:
		print "Row",currentrow,"Coefficient of Power (C_P)    : ",C_Psv 
	print
	print "Row",currentrow,"Coefficient of Thrust (C_T)   : ",CT_formulasv
	print
	#
	if propeller:
		print "Row",currentrow,"Propeller_efficiency          :",abs(eta_Propsv)
		print
	else:
		print "Row",currentrow,"Turbine_efficiency            :",etasv
		print
		print "Row",currentrow,"USV_efficiency                :",eta_HKTsv
		print
	#----------------------------------------------------------------------------------
	return
	
	
#---------------------------------------------------------------------------------------      
#---------------------------------------------------------------------------------------
# Calculating 1D numbers of various properties
def oneDvalues(Rhubsv,Rtipsv,axlIF_all,r_avgsv,r_all,Vztipsv,rho,TorqueBETsv,TorqueBMTsv,propeller):

	print
	print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	print
	axlIF_int = trapez_integral(axlIF_all*r_all,r_all)
	r_int     = trapez_integral(r_all,r_all)
	r_intt    = 0.5*(Rtipsv**2 - Rhubsv**2)
	axlIF_area_avg = axlIF_int/r_int
	print
	print "axlIF_int      : ",axlIF_int
	print "r_int          : ",r_int
	print "r_intt         : ",r_intt
	print "axlIF_area_avg : ",axlIF_area_avg
	CP_from_axlIF = 4.*axlIF_area_avg*(1 - axlIF_area_avg)**2
	CT_from_axlIF = 4.*axlIF_area_avg*(1 - axlIF_area_avg)
	print
	print "CP_from_axlIF  : ",CP_from_axlIF
	print "CT_from_axlIF  : ",CT_from_axlIF
	#--------------------------------------------------------------
	#
	# mdot = trapez_integral(dmdotdr,r)
	print
	print "TorqueBET:", TorqueBETsv, "Nm"
	print "TorqueBMT:", TorqueBMTsv, "Nm"
	print   
	# print "mdot_integration: ", mdot, "kg/s"
	mdot = rho*math.pi*(Rtipsv**2 - Rhubsv**2)*Vztipsv*(1 - axlIF_area_avg)
	print
	print "mdot_formula:", mdot, "kg/s"
	vtheta_out = TorqueBETsv/(r_avgsv*mdot)
	print
	print "Vtheta_out:",vtheta_out,"m/s"
	swirl_power = 0.5*mdot*(vtheta_out**2)
	print
	print "Swirl_Power 1D:", swirl_power/1000.,"kW"
	# Power Coefficient
	print
	# print "CP from integration:",CP
	Power_fluid = 0.5*rho*math.pi*(Rtipsv**2)*(Vztipsv**3)
	print "Power_fluid    :",Power_fluid/1000., "kW"
	
	return Power_fluid
#---------------------------------------------------------------------------------------    

#---------------------------------------------------------------------------------------
# Calculating efficiencies
def calcEfficiencies(rho,Rhubsv,Rtipsv,RPMsv,Lambdasv,Vztipsv,AForceBMTsv,ThrustBETsv,TorqueBETsv,PowerBETsv,Power_fluidsv,propeller):

	if propeller:
		Power       = 2*math.pi*(RPMsv/60.)*TorqueBETsv
		print
		print "Power_propeller: ",Power," W"
		C_P         = (PowerBETsv)/(rho*(((RPMsv)/60.)**3)*(2.*Rtipsv)**5)
		print
		print "C_P:",C_P
		CT_formula = (ThrustBETsv)/(rho*(((RPMsv)/60.)**2)*(2.*Rtipsv)**4)
		print
		print "CT_formula:",CT_formula
		#
		eta_Prop    = ((Vztipsv/(((RPMsv)/60.)*(2.*Rtipsv)))*(CT_formula/C_P)) #(ThrustBETsv*Vztipsv)/PowerBETsv#
		eta         = 0.
		eta_HKT     = 0.
	else:
		C_P         = (PowerBETsv)/(Power_fluidsv)
		print
		print "C_P:",C_P
		##
		CT_formula = ThrustBETsv/(0.5*rho*math.pi*(Rtipsv**2)*(Vztipsv**2))
		print
		print "CT_formula:",CT_formula
		#
		eta_HKT     = PowerBETsv/(ThrustBETsv*Vztipsv)
		print "eta_HKT",eta_HKT
		eta         = C_P/CT_formula
		eta_Prop    = 0.

	return C_P,CT_formula,eta_Prop,eta,eta_HKT

#---------------------------------------------------------------------------------------    
#---------------------------------------------------------------------------------------   
#---------------------------------------------------------------------------------------
#Calculate velocities at different axial stations
def velocitycalc(propeller,Vinsv,axlIFsv):
# Calculates velocities for a Single Row Machine 
	if propeller:
		k = -1
	else: # turbine
		k = +1
	####
	V2sv    = Vinsv*(1. - (k*axlIFsv))
	V3sv    = V2sv
	Vexitsv = Vinsv*(1. - (k*2.*axlIFsv))
		   
	return V2sv,V3sv,Vexitsv
#---------------------------------------------------------------------------------------   
#---------------------------------------------------------------------------------------
#Calculate velocities at different axial stations for a Counter Rotating Machine
def velocitycalc_CR(propeller,Vinsv,axlIF1sv,axlIF2sv,gap):
# Calculates velocities for a Single Row Machine 
	if propeller:
		k = -1
	else: # turbine
		k = +1
	####
	V2CRsv = Vinsv*(1. - (k*axlIF1sv))
	V3CRsv = V2CRsv
	V4CRsv = Vinsv*(1. - (k*gap*axlIF1sv))
	V5CRsv = V4CRsv*(1. - (k*axlIF2sv))
	V6CRsv = V5CRsv
	Vexitsv = V4CRsv*(1. - (k*2.*axlIF2sv))
	if gap == 1.0:
		V3CRsv=V2CRsv
		V4CRsv=V2CRsv
		V5CRsv=V2CRsv
		V6CRsv=V2CRsv
		Vexitsv = Vinsv*(1. - (k*2.*axlIF1sv))
		
	return V2CRsv,V3CRsv,V4CRsv,V5CRsv,V6CRsv,Vexitsv
#---------------------------------------------------------------------------------------   
#---------------------------------------------------------------------------------------

# Calculating spans for various axial locations accounting for streamtube expansion/contraction
def spancalc(V1,V2,V4,rR1):
	# Using continuity principle radius at various axial locations
	# can be calculated.
	# V1A1 = V2A2 = V3A3 = V4A4
	# V2 = V3 
	# V1A1 = V2A2 = V4A4
	# V1(r1^2) = V2(r2^2)  = V4(r4^2) 
	# rR1 is known which is rR1 = r2 
	
	r1 = math.sqrt(math.fabs(V2/V1))*rR1
	r4 = math.sqrt(math.fabs(V2/V4))*rR1
	
	return r1,r4
#---------------------------------------------------------------------------------------   
#---------------------------------------------------------------------------------------

# Calculating spans for various axial locations accounting for streamtube expansion/contraction
def spancalc_CR(V1,V2,V4,V5,V7,rR1,rR2):
	# Using continuity principle radius at various axial locations
	# can be calculated.
	# V1A1 = V2A2 = V3A3 = V4A4 = V5A5 = V6A6 = V7A7
	# V2 = V3 and V5 = V6
	# V1A1 = V2A2 = V4A4 = V5A5 = V7A7
	# V1(r1^2) = V2(r2^2)  = V4(r4^2)  = V5(r5^2)  = V7(r7^2) 
	# rR1 and rR2 are known which is rR1 = r2 and rR2 = r5
	
	r1 = math.sqrt(math.fabs(V2/V1))*rR1
	r4 = math.sqrt(math.fabs(V2/V4))*rR1
	r7 = math.sqrt(math.fabs(V5/V7))*rR2
	
	return r1,r4,r7
#---------------------------------------------------------------------------------------   
#---------------------------------------------------------------------------------------

def aPhi(p,q):
	phirad = p
	propeller,Z,TSR,Rtip,Vz,omega,chrd_current,rcurrent,Cl,Cd = q   
	if propeller:
		C_axial = Cl*math.cos(phirad) - Cd*math.sin(phirad)	
		C_tngnt = Cl*math.sin(phirad) + Cd*math.cos(phirad)
	else: # turbine
		C_axial =  Cl*math.cos(phirad) + Cd*math.sin(phirad)	
		C_tngnt = -Cl*math.sin(phirad) + Cd*math.cos(phirad)
	# ####
	if phirad == 0.:
		phirad = 0.1
	f = Z*(Rtip-rcurrent)/(2.*rcurrent*abs(math.sin(phirad)))
	Ftotal = (2./math.pi)*math.acos(np.exp(-f))
	lambdar = TSR*(rcurrent/Rtip)	
	sigmar  = Z*chrd_current/float(2*math.pi*rcurrent)
	Term1   = 1./(lambdar*math.tan(phirad))
	Nr      = float(sigmar)*(Vz)*(C_tngnt)
	Dr      = 4*Ftotal*omega*(rcurrent)*(math.sin(phirad))**2
	if Dr == 0. :
		Dr = 0.000001
	Term2   = Nr/Dr

	###
	if propeller:
		b  = 1.0/(Term1 + Term2)  
		f1 = float(4.*Ftotal*(b-1) - (sigmar*b*C_axial)/(math.sin(phirad))**2)    
	elif not propeller: # turbine WT/HKT
		# print 'Term1, Term2 in b aPhi:', Term1, Term2
		b  = 1.0/(Term1 - Term2)
		a  = 1.0 - b
		if abs(a) > 0.40: # Using correction to calculate a from Glauert and Spera correction
			a = axlIF_Glauert(a,Ftotal,phirad,sigmar,C_axial)
		f1 = float((4.*Ftotal*a) - (sigmar*(1.0-a)*C_axial)/(math.sin(phirad))**2) 

	return f1
	
#---------------------------------------------------------------------------------------   
#---------------------------------------------------------------------------------------	
def bisection(a1,a3,q,tol):
	propeller,Z,TSR,Rtip,Vz,omega,chrd_current,rcurrent,Cl,Cd = q 
	a2 = float((a1+a3)/2.0)
	while (a2-a1)/2.0 > tol:
		if aPhi(a3,q) == 0:
			return a3
		elif aPhi(a1,q)*aPhi(a3,q) < 0:
			a2 = a3
		else :
			a1 = a3
		a3 = (a1+a2)/2.0
		#print a1,a2,a3
	return a3
	
#---------------------------------------------------------------------------------------   
#---------------------------------------------------------------------------------------
# function for phi including the axial gap between rotors for Contra rotating rotors
def func_phi(argout,argin):
	phirad = argout
	Z,TSR,Rtip,Vz,omega,chrd_current,rcurrent,Cl,Cd,axial_gap = argin
	# print 'phirad:',phirad
	if phirad == 0.:
		phirad = 0.000001  
		print "phirad is zero...so taking it as 0.000001"
	f = Z*(Rtip-rcurrent)/(2.*rcurrent*abs(math.sin(phirad)))
	Ftotal  = (2./math.pi)*math.acos(np.exp(-f))
	lambdar = TSR*(rcurrent/Rtip)	
	sigmar  = Z*chrd_current/float(2*math.pi*rcurrent)
	Term1   = 1./(lambdar*math.tan(phirad))
	Nr      = float(sigmar)*(Vz)*(Cl*math.sin(phirad) - Cd*math.cos(phirad))
	Dr      = 4*Ftotal*omega*(rcurrent)*(math.sin(phirad))**2
	if Dr == 0. :
		Dr = 0.000001
	Term2   = Nr/Dr
	b       = 1.0/(Term1 + Term2)  # axlIF = b - 1
	#
	c  = (b-1)*axial_gap
	f1 = float((Ftotal*c*(2+c)) - (sigmar*(b**2)*(Cl*math.cos(phirad) + Cd*math.sin(phirad)))/(math.sin(phirad))**2)
	
	return f1
#---------------------------------------------------------------------------------------   
#---------------------------------------------------------------------------------------
def root_bisection(a1,a3,argin,tol):
	Z,TSR,Rtip,Vz,omega,chrd_current,rcurrent,Cl,Cd,axial_gap = argin 
	a2 = float((a1+a3)/2.0)
	while (a2-a1)/2.0 > tol:
		if func_phi(a3,argin) == 0:
			return a3
		elif func_phi(a1,argin)*func_phi(a3,argin) < 0:
			a2 = a3
		else :
			a1 = a3
		a3 = (a1+a2)/2.0
		# print a1,a2,a3
	return a3
	
#---------------------------------------------------------------------------------------   
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------   
#---------------------------------------------------------------------------------------
# function for phi including the axial gap between rotors for Contra rotating rotors
def func_phi_closely_placed_rotors(argout,argin):
	phirad = argout
	propeller,Z,TSR,Rtip,Vz,angIF1,axlIF1,omega1,omega2,chrd_current,rcurrent,Cl,Cd,gap = argin
	# print 'phirad:',phirad
	if phirad == 0.:
		phirad = 0.000001  
		print "phirad is zero...so taking it as 0.000001"
	f = Z*(Rtip-rcurrent)/(2.*rcurrent*abs(math.sin(phirad)))
	Ftotal  = (2./math.pi)*math.acos(np.exp(-f))
	#lambdar = TSR*(rcurrent/Rtip)	
	sigmar  = Z*chrd_current/float(2*math.pi*rcurrent)
	# Term1   = Vz/(rcurrent*math.tan(phirad))
	# Nr      = float(sigmar)*(Vz)*(Cl*math.sin(phirad) - Cd*math.cos(phirad))
	# Dr      = 4*Ftotal*(rcurrent)*(math.sin(phirad))**2
	# if Dr == 0. :
		# Dr = 0.000001
	# Term2   = Nr/Dr
	b       = 1.0 + axlIF1 # a1=a2 as assumption(omega2 - (3.0*angIF1*omega1))/float(Term1 + Term2)  # 
	#
	c  = (b-1)*gap
	f1 = float((Ftotal*c*(2+c)) - (sigmar*(b**2)*(Cl*math.cos(phirad) + Cd*math.sin(phirad)))/(math.sin(phirad))**2)
	
	return f1
#---------------------------------------------------------------------------------------   
#---------------------------------------------------------------------------------------
def bisection_closely_placed_rotors(a1,a3,argin,tol):
	propeller,Z,TSR,Rtip,Vz,angIF1,axlIF1,omega1,omega2,chrd_current,rcurrent,Cl,Cd,gap = argin 
	a2 = float((a1+a3)/2.0)
	while (a2-a1)/2.0 > tol:
		if func_phi_closely_placed_rotors(a3,argin) == 0:
			return a3
		elif func_phi_closely_placed_rotors(a1,argin)*func_phi_closely_placed_rotors(a3,argin) < 0:
			a2 = a3
		else :
			a1 = a3
		a3 = (a1+a2)/2.0
		# print a1,a2,a3
	return a3
	
#---------------------------------------------------------------------------------------   
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# function to calculate axial induction factor if a > 0.4 using Glauert's correction
def axlIF_Glauert(axlIF,flosssv,phiradsv,sigmarsv,Casv):

	# print " axial IF  > 0.45 so using Glauert-Spera correction..."
	# print "Current axial IF :", axlIF
	# print
	sphi  = math.sin(phiradsv)
	if flosssv == 0.0:
		flosssv = 0.0001	
	kappa = (sigmarsv*Casv)/(4.*flosssv*(sphi)**2)
	# print "flosssv: ",flosssv
	# print "kappa: ", kappa
	N1   = 2.*flosssv*kappa - ((10./9.) - flosssv)
	A    = (2.*flosssv*kappa)
	B    = (flosssv*((4./3.) - flosssv))
	# if (A-B) < 0. :
		# print 'A,B:',A,B
	N2   = math.sqrt(abs(A - B))
	D1   = 2.*flosssv*kappa - ((25./9.) - (2.*flosssv))
	######
	# print "D1 :",D1    
	if D1 == 0. :
		D1 = 0.000001
	axlIFGlrt = (N1 - N2)/D1 
	# print "Corrected axial IF :", axlIFGlrt
	
	return axlIFGlrt

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------