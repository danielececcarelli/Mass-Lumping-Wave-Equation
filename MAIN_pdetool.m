clear all, close all, clc

addpath masslumping_classes\errors
addpath masslumping_classes\mass_lumping_mesh
addpath masslumping_classes\matrix_assembly
addpath masslumping_classes\quadratures
addpath masslumping_classes\shape_functions
addpath masslumping_classes\triangle_ref
addpath meshgenerator
addpath meshes

%                   Wave equation
%
%          (d^2)u/(dt^2) - div(k*grad(u)) = f in Omega
%                       u=g on Gamma
%                       u=u_0 in Omega, for t=0
%                   du/dt=du_0 in Omega, for t=0
%
%   Omega is set by the user, chosen among some already defined domains,
%   and Gamma is the whole boundary.
%
%   We can use mass lumped elements of degree:
%     - 1 (3 nodes);
%     - 2 (7 nodes);
%     - 3 (12 nodes).

%   Daniele Ceccarelli & Tommaso Missoni - NAPDE project

%% PROBLEM DATA

% Define the problem data

degree = 3;

Time = 20;      % final time
dt = 0.0001;    % time step

u_bc = @(x,y,t) 0.*x.*y; % we don't know the exact solution, but we want homogeneous D b.c.

r = 0.005;
%initialcond = @(x,y) (x.^2+(y-1.25).^2<2*r).*exp((x.^2+(y-1.25).^2)/r^2); % u_0
initialcond = @(x,y) (x.^2+(y-1.25).^2<2*r).*exp(r*(x.^2+(y-1.25).^2)); % u_0
initialcond_time = @(x,y) 0.*x.*y; % du_0

forcetime = @(t) 0;
forcespace = @(x,y) 0;

k = @(x,y) 1; % k

ux = @(x,y,t) 0.*x.*y;
uy = @(x,y,t) 0.*x.*y;

%% MESH

% Establish a mesh of quadratic Lagrange triangles on the domain:
path = 'square';
[p,e,t] = get_mesh_data(path);
T = mesh_pdetool_dirichletBC(p,e,t); %fully Dirichlet BC
T = Refine1(T);
 
T = GenLagrangeMesh2_masslumping(T,degree); % generate the mass lumping mesh
figure(1), ShowMesh_masslumping(T,'r')
savefig('CHAN_1_mesh.fig')

%% RESOLUTION OF THE PROBLEM

x = T.Nodes(T.FNodePtrs,1);          %coordinates of the nodes
y = T.Nodes(T.FNodePtrs,2);

K = Stiffness2_masslumping(T,k);     % Assemble stiffness matrix with mass lumping nodes
M = MassMatrix2_masslumping(T);      % Assemble mass matrix with mass lumping nodes
                                             % M should be diagonal!

t = 0;
% Create the boundary data:
g = getDirichletData_time(T,u_bc,t);
h = getNeumannData2_time(T,ux,uy,k);

% Assemble the load vector:
F = Load2_masslumping(T,forcespace,k,g,h);

% Solve the finite element equations:

u0 = initialcond(x,y); % initial condition

figure(2)
clf
ShowPWPolyFcn2_masslumping(T,u0,g)
title('initial condition')
savefig('CHAN_1_ic.fig')

u0_dt = initialcond_time(x,y); % initial condition for time derivative

% first step
b = (M-0.5*dt^2*K)*u0 + dt*M*u0_dt  + 0.5*dt^2*F.*forcetime(t);
u1 = M\b;

% g = getDirichletData_time(T,u_bc,t);

figure(3)
clf
ShowPWPolyFcn2_masslumping(T,u1,g)
title('evolution in time')
savefig('CHAN_1_start.fig')

for t = dt : dt : Time      % solve the problem for each time step
    
    %g = getDirichletData_time(T,u_bc,t);
    f = forcetime(t); %F = Load2_masslumping_tempo(T,force,k,g,h,t);
    
    b = (dt)^2*(F*f - K*u1) + 2*M*u1 - M*u0;    

    u2 = b./diag(M);        % enough, since M is diag

    % update
    u0 = u1;
    u1 = u2;
    
    % plot of the evolution    
    if(mod(t,0.25)==0)
        figure(3)
        clf
        ShowPWPolyFcn2_masslumping(T,u1,g)
        title(['evolution in time, t = ',num2str(t)])
        savefig(['CHAN_1_',num2str(t),'.fig'])
    end
    %zlim([-0.015, 0.015])
    pause(0.00001)
end


