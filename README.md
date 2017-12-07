# py_BEM
Python based design analysis tool for all horizontal axis unducted rotors including propellers and turbines using Blade Element Momentum Theory.

A robust design and analysis tool using the blade element momentum theory is created to obtain spanwise stagger and chord definitions as an initial design of 3D blade generation for unducted rotors both solo and contra-rotating configurations. Designers of Propellers, Wind Turbines, Hydrokinetic and Tidal Turbines, Helicopter rotors and UAV/Quad/Octo-copter rotors can greatly benefit from pyBEM. Trade off studies with power, thrust, weight, rotor disc area, number of blades, hover performance can be assessed very quickly for any mission.

• Written in Python27
• Design and analysis tool for Propellers, Wind Turbines, Hydrokinetic and Tidal Turbines, Helicopter rotors in hover
• Solo and Contra rotating configuration design capability
• Airfoil characteristics using XFOIL http://web.mit.edu/drela/Public/web/xfoil/
• Local effect of Reynolds Number is incorporated which varies with angle of attack and provides a choice of lift to drag ratio
• Spanwise variation of lift and drag coefficient, lift to drag ratio expands the design space to obtain a better optimum
• Spline routines written for T-Blade3 [1] are also utilized here for spanwise definition of chord and other properties
• Airfoil lift and drag extrapolation at higher angles of attack is also implemented to capture performance in stall conditions
• Stall Delay, Rotational effect on lift coefficient, compressibility effect, Tip loss model, turbulent wake axial induction factor model are also implemented
• Three types of contra rotating configurations are included which affects the rotor-rotor interaction
• Spanwise twist and stagger are the outputs if chord is defined. Spanwise lift, drag and glide ratio can also be defined
• 3D blade input file is also created to connect to high fidelity design and analysis tools

Any desiginer of unducted rotors solo or contra rotating configuration can use this tool to design and analyze their
unducted rotors. Quickest way to get into 3D high fidelity design space using py_BEM and T-Blade3 https://github.com/GTSL-UC/T-Blade3 written by the author as part of Masters research and improved it in collaboration with other co-researchers.
• Air and Water Propellers
• UAVs, Multirotor applications
• Wind Turbines
• Tidal and Hydrokinetic Turbines
• Helicopter Rotors, mainly in hover/vertical ascent
