function g=getDirichletData_time(T,u,t)

% g=getDirichletData_time(T,u,t)
%
%   This function assembles the Dirichlet data from
%   the mesh T and function u at time t; that is, it evaluates
%   u(x,y,t) at time t at the constrained nodes of T.
%
%   The output is the Nc by 1 array g, where Nc is the
%   number of constrained nodes.
%
%   For a description of the data structure T, see
%   "help Mesh1".

%   Modified starting from "Understanding and Implementing the Finite
%   Element Method" by Mark S. Gockenbach (copyright SIAM 2006).

%   Daniele Ceccarelli & Tommaso Missoni - NAPDE project

x=T.Nodes(T.CNodePtrs,1);
y=T.Nodes(T.CNodePtrs,2);
g=feval(u,x,y,t);
