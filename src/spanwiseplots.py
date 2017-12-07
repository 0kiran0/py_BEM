## Plotting routine
import sys,os, errno
import math, timeit
import numpy as np
import pylab as py
from functions import * 
from matplotlib.pyplot import *
import mpl_toolkits.mplot3d.axes3d as p3

workdir = os.getcwd()
if plot:
	if os.name == 'nt':
		plotdir = workdir + '\plots'
	elif os.name == 'posix':
		plotdir = workdir + '/plots'
	## Check if the plot directory exist and if not then make one
	make_sure_path_exists(plotdir)
	####
	print nbrow
##
#########################################################################
def plots(propeller,nbrow,arguments,REYNfront,REYNrear):

	print "Plotting radial plots..."
	print   
	# print arguments[0]
	###################################
	
	# ## Axial Velocity vs span ###
	fignum = 1
	py.figure(fignum, figsize=(10, 6))
	if nbrow == 1:
		py.plot(Vz1,span1,'-o',color = 'black',markersize=7,lw=2)    
		py.plot(V2,span,'-d',color = 'black',markersize=7,lw=2)    
		py.plot(Vexit,span4,'-s',color = 'black',markersize=7,lw=2)
		py.title("Velocities at various axial stations [m/s] vs span")
		if propeller:
			py.xlabel("Axial Velocities for Prop")
		else:
			py.xlabel("Axial Velocities for Turbine")
		py.ylabel("span")
		leg = py.legend(["Vz1","V_LE = V_TE","Vexit"], loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0.)
	elif nbrow == 2:
		py.plot(Vz1,span,'-o',color = 'black')    
		py.plot(V2CR,span,'-o',color = 'black',markersize=7,mfc ='none',mew=2,lw=2)    
		py.plot(V4CR,span4CR,'-+',color = 'black',markersize=7,lw=2)    
		py.plot(V5CR,span5CR,'-s',color = 'black',markersize=7,lw=2)    
		py.plot(VexitCR,span7CR,'-^',color = 'black',markersize=7,lw=2)
		py.title("Velocities at various axial stations [m/s] vs span")
		if propeller:
			py.xlabel("Axial Velocities for CR-Prop")
		else:
			py.xlabel("Axial Velocities for CR-Turbine")
		py.ylabel("span")
		leg = py.legend(["Vz_R1","V_LE_R1 = V_TE_R1","Vz_R2","V_LE_R2 = V_TE_R2","Vexit"], loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0.)        
	if propeller:
		py.ylim([0.0,1.0])    
	else:
		py.ylim([0.0,1.0]) #py.xlim([0.0,10.0])
	py.savefig(os.path.join(plotdir, 'Vz.png'), bbox_extra_artists=(leg,),bbox_inches='tight')      

	# ## Induction Factors vs span ###    
	fignum += 1
	py.figure(fignum, figsize=(10, 6))
	py.plot(axlIFnew,span,'-o',color = 'blue',markersize=7,lw=2)   
	py.plot(angIFnew,span,'-s',color = 'red',markersize=7,lw=2)
	if nbrow == 2:
		py.plot(axlIFnew2,span5CR,'-d',color = 'black',markersize=7,lw=2)   
		py.plot(angIFnew2,span5CR,'-o',color = 'magenta',markersize=7,lw=2)
		py.plot(angIFnew12,span5CR,'-+',color = 'green',markersize=7,lw=2)
		py.title("a, a' vs span for front and aft Rotor")
		py.xlabel("a,a',a2,a2',a12'")
		leg = py.legend(["a1","a1'","a2","a2'","a12'"], loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0.)
	elif nbrow == 1:
		py.title("a, a' vs span ")
		py.xlabel("a,a'")
		leg = py.legend(["a","a'"], loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0.)
	py.ylabel("span")
	py.savefig(os.path.join(plotdir, 'IFs.png'),bbox_extra_artists=(leg,),bbox_inches='tight')

	# ## Actual Chord (m) vs span ###
	fignum += 1		  
	py.figure(fignum, figsize=(10, 6))    
	py.plot(chordr,span,'-o',color = 'blue',markersize=7,lw=2)
	if nbrow == 2:
		py.plot(chordr2,span5CR,'-o',color = 'red',markersize=7,lw=2)
		py.title("chordr [m] vs span for front and Aft Rotor")
		py.xlabel("chordr, chordr2 [m]")
		leg = py.legend(['chordr Front Rotor [m]','chordr2 Aft Rotor [m]'], loc=2, bbox_to_anchor=(0.6, 1), borderaxespad=0.)
	else:
		py.title("chordr [m] vs span")
		py.xlabel("chordr [m]")
		leg = py.legend(['chordr[m]'], loc=2, bbox_to_anchor=(0.6, 1), borderaxespad=0.)
	py.ylabel("span")
	if propeller:
		py.axis('equal')
		py.axis([0.0,1.0,0.0,1.2])
	else:
		py.axis('equal')
		py.axis([0.0,1.0,0.0,1.2])
		#py.xlim([0.0,0.1])
	py.savefig(os.path.join(plotdir, 'chord.png'), bbox_extra_artists=(leg,), bbox_inches='tight')

	# ## Stagger vs span ###			  
	fignum += 1	  
	if nbrow == 1:          
		py.figure(fignum, figsize=(10, 6))
		# py.plot(inbeta,span,color = 'red')    
		py.plot(stggr,span,'-o',color = 'blue',markersize=7,lw=2)
		py.title("Stagger vs span ")
		py.xlabel("Stagger")
		py.ylabel("span")
		# py.axis('equal')
		leg = py.legend(['Stagger'], loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0.)
		py.savefig(os.path.join(plotdir, 'stagger.png'), bbox_extra_artists=(leg,), bbox_inches='tight')			  
	if nbrow == 2:
		fignum += 1
		py.figure(fignum, figsize=(10, 6))
		py.plot(stggr,span,'-o',color = 'blue',markersize=7,lw=2)    
		py.plot(stggr2,span5CR,'-s',color = 'red',markersize=7,lw=2)
		py.title("Stagger vs span for Front and Aft Rotors")
		py.xlabel("Stagger for Front and Aft Rotor")
		py.ylabel("span")
		# py.axis('equal')
		leg = py.legend(['Stagger Front','Stagger Aft'], loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0.)
		py.savefig(os.path.join(plotdir, 'stagger_CR.png'), bbox_extra_artists=(leg,), bbox_inches='tight')

	# ## dThrustdr BET and BMT vs span ###        
	fignum += 1		  
	py.figure(fignum, figsize=(10, 6))    
	py.plot(dThrustBETdr,span,'-o',color = 'blue',markersize=7,lw=2)
	py.plot(dAForceBMTdr,span,'-^',color = 'blue',markersize=7,lw=2)
	if nbrow == 2:
		py.plot(dThrustBET2dr,span5CR,'-o',color = 'red',markersize=7,lw=2)
		py.plot(dAForceBMT2dr,span5CR,'-^',color = 'red',markersize=7,lw=2)        
		py.title("Counter-Rotating: Axial-Force/m (N/m) vs span (Blade Element theory and BMT)")
		py.xlabel("dThrustBETdr , dThrustBET2dr, dAForceBMTdr , dAForceBMT2dr in N/m")
		leg = py.legend(['dThrustBETdr','dAForceBMTdr','dThrustBET2dr','dAForceBMT2dr'], loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0.)
	else:
		py.title("Axial-Force/m (N/m) vs span (Blade Element theory and BMT)")
		py.xlabel("dThrustBETdr, dAForceBMTdr in N/m")
		leg = py.legend(['dThrustBETdr','dAForceBMTdr'], loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0.) 
	py.ylabel("span")
	py.savefig(os.path.join(plotdir, 'dThrust_dr.png'), bbox_extra_artists=(leg,), bbox_inches='tight')
		  
	# ## dTorquedr BET and BMT vs span ###   			  
	fignum += 1
	py.figure(fignum, figsize=(10, 6))    
	py.plot(dTorqueBETdr,span,'-o',color = 'blue',markersize=7,lw=2)
	py.plot(dTForceBMTdr,span,'-^',color = 'blue',markersize=7,lw=2)    
	if nbrow == 2:
		py.plot(dTorqueBET2dr,span5CR,'-o',color = 'red',markersize=7,lw=2)
		py.plot(dTForceBMT2dr,span5CR,'-^',color = 'red',markersize=7,lw=2)        
		py.title("Counter-Rotating: dTorque/dr (Nm/m) vs span (Blade Element theory and BMT)")
		py.xlabel("dTorqueBETdr, dTorqueBET2dr,dTForceBMTdr, dTForceBMT2dr in Nm/m")
		leg = py.legend(['dTorqueBETdr','dTForceBMTdr','dTorqueBET2dr','dTForceBMT2dr'], loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0.)
	else:
		py.title("dTorque/dr (Nm/m) vs span (Blade Element theory and BMT)")
		py.xlabel("dTorqueBETdr,dTForceBMTdr in Nm/m")
		leg = py.legend(['dTorqueBETdr','dTorqueBMTdr'], loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0.)         
	py.ylabel("span")
	py.savefig(os.path.join(plotdir, 'dTorque_dr.png'), bbox_extra_artists=(leg,), bbox_inches='tight') 

	# ## dPowerdr BET and BMT vs span ###     
	fignum += 1
	py.figure(fignum, figsize=(10, 6))    
	py.plot(dPowerBETdr,span,'-o',color = 'blue',markersize=7,lw=2)
	py.plot(dPowerBMTdr,span,'-^',color = 'blue',markersize=7,lw=2)
	if nbrow == 2:
		py.plot(dPowerBET2dr,span5CR,'-o',color = 'red',markersize=7,lw=2)
		py.plot(dPowerBMT2dr,span5CR,'-^',color = 'red',markersize=7,lw=2)        
		py.title("Counter Rotating: dPower/dr (W/m) vs span (Blade Element theory and BMT)")
		py.xlabel("dPowerBETdr, dPowerBMTdr, dPowerBET2dr, dPowerBMT2dr in W/m")
		leg = py.legend(['dPowerBETdr','dPowerBMTdr','dPowerBET2dr','dPowerBMT2dr'], loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0.)         
	else:
		py.title("dPower/dr (W/m) vs span (Blade Element theory and BMT)")
		py.xlabel("dPowerBETdr,dPowerBMTdr in W/m")
		leg = py.legend(['dPowerBETdr','dPowerBMTdr'], loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0.)
	py.ylabel("span")
	py.savefig(os.path.join(plotdir, 'dPower_dr.png'), bbox_extra_artists=(leg,), bbox_inches='tight')

	# ## REYN vs span ###     
	fignum += 1
	py.figure(fignum, figsize=(10, 6))    
	py.plot(REYNfront,span,'-o',color = 'blue',markersize=7,mfc ='none',mew=2,lw=2)
	if nbrow == 2:
		py.plot(REYNrear,span5CR,'-^',color = 'red',markersize=7,mfc ='none',mew=2,lw=2)
		py.title("Counter Rotating: REYN vs span ")
		py.xlabel("REYN for front and aft rotors")
		leg = py.legend(['REYN front rotor','REYN aft rotor'], loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0.)         
	else:
		py.title("REYN vs span")
		py.xlabel("REYN for the rotor")
		leg = py.legend(['REYN'], loc=2, bbox_to_anchor=(0.8, 1), borderaxespad=0.)
	# py.axis('equal')	
	py.ylabel("span")
	py.ticklabel_format(style='sci', axis='x', scilimits=(0,0))    
	py.savefig(os.path.join(plotdir, 'REYN.png'), bbox_extra_artists=(leg,), bbox_inches='tight')
	
	if varyingREYN or spanwise:

		# ## AoA vs span ### 
		fignum += 1	
		py.figure(fignum, figsize=(10, 6))    
		py.plot(arguments[0],span,'-o',color = 'blue',markersize=7,mfc ='none',mew=2,lw=2)
		if nbrow == 2:
			py.plot(arguments[1],span5CR,'-^',color = 'red',markersize=7,mfc ='none',mew=2,lw=2)
			py.title("Counter Rotating: Angle of Attack (AoA) vs span ")
			py.xlabel("AoA for front and aft rotors")
			leg = py.legend(['AoA front rotor','AoA aft rotor'], loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0.)         
		elif nbrow == 1:
			py.title("Angle of Attack (AoA) vs span")
			py.xlabel("AoA for the rotor")
			leg = py.legend(['AoA'], loc=2, bbox_to_anchor=(0.8, 1), borderaxespad=0.)
		py.ylabel("span")
	    #py.axis('equal')			
		py.ticklabel_format(style='sci', axis='x', scilimits=(0,0))    
		py.savefig(os.path.join(plotdir, 'AoA.png'), bbox_extra_artists=(leg,), bbox_inches='tight')

		# ## L/D vs span ### 
		fignum += 1		
		py.figure(fignum, figsize=(10, 6))    
		py.plot(abs(L_over_D),span,'-o',color = 'blue',markersize=7,mfc ='none',mew=2,lw=2)
		if nbrow == 2:
			py.plot(abs(L_over_D2),span5CR,'-^',color = 'red',markersize=7,mfc ='none',mew=2,lw=2)
			py.title("Counter Rotating: L/D vs span ")
			py.xlabel("L/D for front and aft rotors")
			leg = py.legend(['L/D front rotor','L/D aft rotor'], loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0.)         
		else:
			py.title("L/D vs span")
			py.xlabel("L/D for the rotor")
			leg = py.legend(['L/D'], loc=2, bbox_to_anchor=(0.8, 1), borderaxespad=0.)
		py.ylabel("span")
		# py.axis('equal')
		py.ticklabel_format(style='sci', axis='x', scilimits=(0,0))    
		py.savefig(os.path.join(plotdir, 'L_over_D.png'), bbox_extra_artists=(leg,), bbox_inches='tight')

	return