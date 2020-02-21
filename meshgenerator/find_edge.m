function edge = find_edge(node1,node2,edges)
% edge = find_edge(node1,node2,edges)
% node1 and node2 are the labels of two nodes, edges is a set of edges of
% the mesh
% 
% function "find_edge" gives the label of the edge defined by node1 and
% node2 (with the right sign, in order to mantain a certain order)

% Daniele Ceccarelli & Tommaso Missoni - NAPDE project

edge = 0;

[~,n] = size(edges);

for i = 1:n
    if(edges(2,i) == node1 && edges(3,i) == node2)
        edge = edges(1,i);
    end
    if(edges(2,i) == node2 && edges(3,i) == node1)
        edge = -edges(1,i);
    end
end

end

