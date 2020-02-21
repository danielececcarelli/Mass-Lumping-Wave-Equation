function [a,b] = search_triangle(nodes,t)
% t is the set of all the triangles of the mesh, nodes is a vector of two
% node labels
% 
% the function "search_triangoli" gives as output the labels of the two
% triangles which are adjacent to the edge defined by the nodes in nodes
% set

% Daniele Ceccarelli & Tommaso Missoni - NAPDE project


    [~,n] = size(t);
    a = 0;
    b = 0;
    for i = 1:n
        if(nodes(1)==t(1,i) || nodes(1)==t(2,i) || nodes(1)==t(3,i))
            if(nodes(2)==t(1,i) || nodes(2)==t(2,i) || nodes(2)==t(3,i))
                if(a==0)
                    a = i;
                else
                    b = i;
                end
            end
        end
    end

end

