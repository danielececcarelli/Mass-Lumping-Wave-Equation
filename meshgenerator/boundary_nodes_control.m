function result = boundary_nodes_control(node,boundary_nodes)
% result = boundary_nodes_control(node,boundary_nodes)

% node is the label of a node, boundary_nodes is the set of the boundary nodes
% of the mesh
%
% the function "boundary_nodes_control" gives 1 if the node is a
% boundary nodes, 0 otherwise

% Daniele Ceccarelli & Tommaso Missoni - NAPDE project


result = 0;

[~,n] = size(boundary_nodes);

for i = 1:n
    if(boundary_nodes(i) == node)
        result = 1;
    end
end

end

