function [result] = edge_control(nodes1,nodes2,edges)
% [result] = edge_control(nodes1,nodes2,edges)
% nodes1 and nodes2 are the labels of two nodes, edges is a certain set of
% edges of the mesh
%
% the function "control" gives 0 if the edge defined by nodes1 and nodes2
% IS ALREADY PRESENT in edges, 1 otherwise

% Daniele Ceccarelli & Tommaso Missoni - NAPDE project

[~,n] = size(edges);

result = 1;

for i = 1:n
    if(nodes1 == edges(2,i) && nodes2 == edges(3,i))
        result = 0;
    end
    if(nodes2 == edges(2,i) && nodes1 == edges(3,i))
        result = 0;
    end
end

end

