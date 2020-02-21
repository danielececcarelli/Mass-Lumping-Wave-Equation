# Mass Lumped triangular elements for the Wave Equation - MATLAB code

- by Daniele Ceccarelli & Tommaso Missoni - NAPDE project

## How to run the simulation:

1) Choose one of the MAIN.m:
	- a) MAIN_rectangle for a simulation on a 2D rectangle
	- b) MAIN_pdetool for a simulation on a 2D mesh coming from pde_tool draw (you can also build your own mesh - see later *)
	- c) MAIN_convergence_h for a convergence analysis in space (rectangle mesh)
	- d) MAIN_convergence_dt for a convergence analysis in space (rectangle mesh)

2) Choose your data (for example Time, dt, degree, k(x,y), f(x,y,t)); 
	if you are in case 1.a), to change the the mesh size h you need to change n in "T = RectangleMeshD1(n)"
	while, if you are in case 1.b), you have to change the mesh size with a refining in pde_tool directly or using the function "Refine.m"

3) if you are in case 1.b), to change the mesh you have to change path = 'narrow_channel' with a new string
	indicating the folder where you have saved the mesh from pde_tool

4) Run the code

*: if you want to use a new mesh (with Dirichlet Boundary), you can draw it using pde_tool, export it 
(in p.xlsx, t.xlsx, e.xlsx), and put it in a new folder "my_mesh" inside "meshes" folder. 
To call it, follow point 3), you have to change path = 'narrow_channel' with path = 'my_mesh' and automatically 
the code will read it.


## CODE ORGANIZATION:

1) masslumping_classes : the part of the project concerning the mass lumping technique; it is subdivided in another 6 folders:
	1. errors (H1 and L2, with Dunavant quadrature or Mass Lumping quadrature)
	2. mass_lumping_mesh (Generate the ML lagrange mesh of degree d)
	3. matrix_assembly (assembly of K, M and the rhs + Dirichlet and Neumann b.c.)
	4. quadratures (Mass Lumping quadrature(qpts_qwts) and Dunavant quadrature)
	5. shape_functions (ML FE space)
	6. triangle_ref (Ref.Triangle, get_nodes and transformation from a generic triangle to reference triangle)

2) meshes : some meshes that we draw using pde_tool and that we used in our numerical analysis

3) meshgenerator : code to link the output of pde_tool export (p,t,e) to our mesh structure, in order to use it in our algorithm

4) non_ML_wave : code that we used to compute Finite Element comparison (wave equation). It is basically the "elliptic_code" with 
		 a leap-frog scheme to manage time-step.

5) elliptic_code: the code that we use as starting point. It's a classical FE method for Elliptic Problem
		  with higher-order polynomial degree 	["Understanding and Implementing the Finite
		  Element Method" by Mark S. Gockenbach (copyright SIAM 2006).
		  
## Report and Slides

If you want to see more detail of the project, you can see the extended [report](https://github.com/danielececcarelli/Mass-Lumping-Wave-Equation/blob/master/Report%20and%20presentation/Ceccarelli_Missoni_Report.pdf)  or the final [slides](https://github.com/danielececcarelli/Mass-Lumping-Wave-Equation/blob/master/Report%20and%20presentation/Ceccarelli_Missoni_Slides.pdf).
