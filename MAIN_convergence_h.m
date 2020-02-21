clear all, close all, clc

addpath masslumping_classes\errors
addpath masslumping_classes\mass_lumping_mesh
addpath masslumping_classes\matrix_assembly
addpath masslumping_classes\quadratures
addpath masslumping_classes\shape_functions
addpath masslumping_classes\triangle_ref
addpath meshgenerator
addpath meshes

% CONVERGENCE IN SPACE

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

Time = 0.2;   % final time
dt = 0.00001; % time step

initialcond = @(x,y) x.*(1-x).*y.*(1-y); % u_0
initialcond_time = @(x,y) -x.*(1-x).*y.*(1-y); % du_0

% force = @(x,y,t) (x.*(1-x).*y.*(1-y)+2*(y-y.^2+x-x.^2)).*exp(-t); % f
forcespace = @(x,y) 2*(-x.^2+x)+2*(-y.^2+y)+(-x.^2+x).*(-y.^2+y);
forcetime = @(t) exp(-t);

k = @(x,y) 1; % k

% u is also used to compute the boundary values;
% also, these functions are used to check the accuracy of the computed
% solution in L^2 and H^1 norm: you can set the partial derivatives to
% whatever value you want if you are not interested in computing the H^1
% error, and you can set u in a way that it satisfies the b.c. you want,
% but not being the actual solution, if you are not interested in computing
% the L^2 norm

u = @(x,y,t) x.*y.*(1-x).*(1-y).*exp(-t);
u_end = @(x,y) x.*y.*(1-x).*(1-y).*exp(-Time); % u evaluated in t=Time
ux = @(x,y,t) (1-2*x).*y.*(1-y).*exp(-Time); % ux evaluated in t=Time
uy = @(x,y,t) x.*(1-x).*(1-2*y).*exp(-Time); % uy evaluated in t=Time

N = 7;
L2_err = zeros(1,N);
H1_err = zeros(1,N);
hh = 1./2.^[1:N];

for n = 1:N
    
    n
    
    %% MESH
    
    % Establish a mesh of Lagrange triangles on the domain:
    
    T0 = RectangleMeshD1(2^n); %fully Dirichlet BC
    
    T = GenLagrangeMesh2_masslumping(T0,degree); % generate the mass lumping mesh
   
%     figure(1), ShowMesh_masslumping(T,'r'), title('mesh')
%     figure(2), ShowMesh(TH,'r'), title('mesh')
    
    %% RESOLUTION OF THE PROBLEM
    
    x = T.Nodes(T.FNodePtrs,1);          %coordinates of the nodes
    y = T.Nodes(T.FNodePtrs,2);
    
    K = Stiffness2_masslumping(T,k);     % Assemble stiffness matrix with mass lumping nodes
    M = MassMatrix2_masslumping(T);      % Assemble mass matrix with mass lumping nodes
    % M should be diagonal!
    
    t = 0;
    % Create the boundary data:
    g = getDirichletData_time(T,u,t);
    h = getNeumannData2_time(T,ux,uy,t,k);
    
    % Assemble the load vector:
    F = Load2_masslumping(T,forcespace,k,g,h); %F = Load2_masslumping_tempo(T,force,k,g,h,t);
    
    % Solve the finite element equations:
    
    u0 = initialcond(x,y); % initial condition
    
%     figure(2)
%     clf
%     ShowPWPolyFcn2_masslumping(T,0.00001*u0,g)
%     title('initial condition')
    
    u0_dt = initialcond_time(x,y); % initial condition for time derivative
    
    % first step
    b = (M-0.5*dt^2*K)*u0 + dt*M*u0_dt  + 0.5*dt^2*F.*forcetime(t);
    u1 = M\b;
%     figure(3)
%     clf
%     ShowPWPolyFcn2_masslumping(T,u1,g)
%     title('evolution in time')
    
    for t = dt : dt : Time      % solve the problem for each time step
        
        g = getDirichletData_time(T,u,t);
        %h = getNeumannData2_time(T,ux,uy,t,k);
        f = forcetime(t);
        
        b = (dt)^2*(F*f - K*u1) + 2*M*u1 - M*u0;
        
        u2 = b./diag(M);        % enough, since M is diag
        
        % update
        u0 = u1;
        u1 = u2;
        
%         %     % plot of the evolution
%         figure(3)
%         clf
%         ShowPWPolyFcn2_masslumping(T,u1,g)
%         title(['evolution in time, t = ', num2str(t)])
%         %zlim([-0.1, 0.1])
%         pause(0.0000000001)
    end
    
    % compute L2 err at final step
    L2_err(n) = L2NormErr2_dunavant(T,u_end,u2,g);
    H1_err(n) = H1EnergyNormErr_dunavant(T,1,ux,uy,u2,g);
end

figure(1)
loglog(hh,L2_err,'.-','markersize',14,'Color','g'), hold on, grid on
%loglog(hh,hh,'o-','Color','r')
%loglog(hh,hh.^2,'*-','Color','b')
loglog(hh,hh.^3,'^-','Color','k')
loglog(hh,hh.^4,'d-','Color','m')
%legend('L^2 error','h','h^2','h^3','h^4','location','best')
legend('L^2 error','h^3','h^4','location','best')
xlabel('h'), ylabel('|| u-u_h ||_{L^2}')

figure(2)
loglog(hh,H1_err,'.-','markersize',14,'Color','g'), hold on, grid on
%loglog(hh,hh,'o-','Color','r')
loglog(hh,hh.^2,'*-','Color','b')
loglog(hh,hh.^3,'^-','Color','k')
%loglog(hh,hh.^4,'d-','Color','m')
%legend('H^1 error','h','h^2','h^3','h^4','location','best')
legend('H^1 error','h^2','h^3','location','best')
xlabel('h'), ylabel('|| u-u_h ||_{H^1}')

format short
rates_L2 = log(L2_err(1:end-1)./L2_err(2:end))./log(hh(1:end-1)./hh(2:end))
rates_H1 = log(H1_err(1:end-1)./H1_err(2:end))./log(hh(1:end-1)./hh(2:end))

format long e
L2_err
H1_err


