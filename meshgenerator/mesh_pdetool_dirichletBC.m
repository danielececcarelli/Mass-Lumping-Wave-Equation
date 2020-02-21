function [T] = mesh_pdetool_dirichletBC(p,e,t)
% [T] = mesh_pdetool_dirichletBC(p,e,t)
% p, e and t define a mesh which has been made with pdetool: t contains the
% informations about the triangles, e about the (boundary) edges and p
% about the nodes
%
% this function, starting from the data from pdetool, produces a mesh which
% is equivalent to the original one, but through structures which are
% suitable for our code to be used

% Daniele Ceccarelli & Tommaso Missoni - NAPDE project

addpath masslumping_classes
addpath meshgenerator

% mesh_pdetool_prova % per visualizzare la mesh (e volendo anche per estrarla)

boundary_nodes = [e(1,:) , e(2,:)];
boundary_nodes = unique(boundary_nodes);
% boundary_node is now a vector which contains the 
% labels of all the boundary nodes

[~,nt] = size(t); % nt = number of triangles in the mesh

edges = zeros(3,nt*3);
curr_edge = 1;
for i = 1:nt
    curr_triang = t(1:3,i);
    for j = 1:3
        if(j~=3)
            if(edge_control(curr_triang(j),curr_triang(j+1),edges))
                edges(1,(i-1)*3 + j) = curr_edge;
                edges(2:3,(i-1)*3 + j) = [curr_triang(j); curr_triang(j+1)];
                curr_edge = curr_edge+1;
            end
        else
            if(edge_control(curr_triang(j),curr_triang(1),edges))
                edges(1,(i-1)*3 + j) = curr_edge;
                edges(2:3, (i-1)*3 + j) = [curr_triang(j); curr_triang(1)];
                curr_edge = curr_edge+1;
            end
        end
    end
end

nl = max(edges(1,:)); % nl = number of all the edges in the mesh
edges_new = zeros(3,nl);
for i = 1:3*nt
    if(edges(1,i)~=0)
        edges_new(1:3,edges(1,i)) = edges(1:3,i);
    end
end
% edges_new is now a 3xnl matrix; first row: label of the edge, second and
% third rows: labels of the two vertices

near_element_edges = zeros(2,nl);
for i = 1:nl
    nodes = [edges_new(2,i) , edges_new(3,i)];
    [a,b] = search_triangle(nodes,t);
    near_element_edges(:,i) = [a;b];
end
% near_element_edges is a 2xnl matrix: the i-th column contains the labels
% of the two triangles adjacent to the i-th edge

edges_of_a_element = zeros(3,nt);
for i = 1:nt
    current_nodes = t(1:3,i);
    edges_of_a_element(1,i) = find_edge(current_nodes(1),current_nodes(2),edges_new);
    edges_of_a_element(2,i) = find_edge(current_nodes(2),current_nodes(3),edges_new);
    edges_of_a_element(3,i) = find_edge(current_nodes(3),current_nodes(1),edges_new);
end
% edges_of_a_element is a 3xnt matrix: the i-th column contains the labels
% (with "sign") of the 3 edges of the i-th triangle

T = mesh_generator(p,edges_new,t,boundary_nodes,near_element_edges,edges_of_a_element);


