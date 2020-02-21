function [T] = mesh_generator(p,edges,t,boundary_nodes,near_element_edges,edges_of_a_element)
% [T] = mesh_generator(p,edges,t,boundary_nodes,near_element_edges,edges_of_a_element)
% In this code, we use the input we compute during mesh_pdetool_dirichletBC
% to build a mesh structure for our code

% Daniele Ceccarelli & Tommaso Missoni - NAPDE project

T.Degree = 1;

[~,Nv] = size(p);
T.Nodes = zeros(Nv,2);
T.NodePtrs = zeros(Nv,1);

% Compute the number of constrained nodes and allocate space for CNodePtrs:
[~,Nc] = size(boundary_nodes);
T.CNodePtrs = zeros(Nc,1);

% Compute the number of free nodes and allocate space for FNodePtrs:
Nf = Nv-Nc;
T.FNodePtrs = zeros(Nf,1);

% Compute the number of triangles and allocate space for Elements:
[~,Nt] = size(t);
T.Elements = zeros(Nt,3);

% Compute the number of edges and allocate space for Edges, EdgeEls,
% and EdgeCFlags:
Ne = max(edges(1,:));
T.Edges = zeros(Ne,2);
T.EdgeEls = zeros(Ne,2); %puntatori agli elementi adiacenti al lato (zero se bordo)
T.EdgeCFlags = zeros(Ne,1);

% Compute the number of boundary edges and allocate space for FBndyEdges:
Nb = 0;
T.FBndyEdges = zeros(Nb,1);

% Loop over the rows and columns of the mesh, defining the nodes, node
% pointers, free node pointers, and constrained node pointers.
%
% k is the number of the node.
% kf is the number of the free node.
% kc is the number of the constrained node.

k = 0;
kf = 0;
kc = 0;

% Loop over the rows of the grid
for i = 1:Nv
      y = p(2,i); 
      x = p(1,i);

      % Insert the coordinates of the node
      T.Nodes(i,:) = [x , y];

      % Interior nodes are free
      if (~boundary_nodes_control(i,boundary_nodes))
         kf = kf+1;
         T.NodePtrs(i) = kf;
         T.FNodePtrs(kf) = i;
      else
         kc = kc+1;
         T.NodePtrs(i) = -kc;
         T.CNodePtrs(kc) = i;
      end
end

for i = 1:Ne
   T.Edges(i,:) = [edges(2,i), edges(3,i)];
   T.EdgeEls(i,:) = [near_element_edges(1,i), near_element_edges(2,i)];
end

T.Elements = edges_of_a_element';

end
