function [M] = eval_matrix_masslumping(nodes, qpts)

%   [M] = eval_matrix_masslumping(nodes, qpts)
%   evaluation of the x,y part of shape functions (the coeff will be considered)
%   thanks to the phi matrix (filled with fun coeff_basis) on qpts.
%   We use the input nodes to understand which degree we have:
%            size nodes = 3 -> degree = 1
%            size nodes = 7 -> degree = 7
%            size nodes = 12 -> degree = 3

%   Daniele Ceccarelli & Tommaso Missoni - NAPDE project

if(nargin < 2)
    qpts = nodes;
end

[n, ~] = size(nodes);
[m, ~] = size(qpts);

M = zeros(m,n);

if(n == 3)
    for i=1:m
        M(i,1) = 1;
        M(i,2) = qpts(i,1); %x
        M(i,3) = qpts(i,2); %y
    end
    
elseif(n == 7)
    for i=1:m
        M(i,1) = 1;
        M(i,2) = qpts(i,1); %x
        M(i,3) = qpts(i,2); %y
        M(i,4) = qpts(i,1)^2;   %x^2
        M(i,5) = qpts(i,2)^2;   %y^2
        M(i,6) = qpts(i,1)*qpts(i,2);  %x*y
        M(i,7) = qpts(i,1)*qpts(i,2)*(1-qpts(i,1)-qpts(i,2)); %bubble = x*y*(1-x-y)
    end
else %%% n == 12
    for i=1:m
        M(i,1) = 1;
        M(i,2) = qpts(i,1); %x
        M(i,3) = qpts(i,2); %y
        M(i,4) = qpts(i,1)^2;   %x^2
        M(i,5) = qpts(i,2)^2;   %y^2
        M(i,6) = qpts(i,1)*qpts(i,2);  %x*y
        M(i,7) = qpts(i,1)^3;   %x^3
        M(i,8) = qpts(i,2)^3;   %y^3
        M(i,9) = qpts(i,1)^2*qpts(i,2);  %x^2*y
        M(i,10) = qpts(i,1)*qpts(i,2)^2; %x*y^2
        M(i,11) = qpts(i,1)^2*qpts(i,2)*(1-qpts(i,1)-qpts(i,2));   %bubble1 x^2*y*(1-x-y)
        M(i,12) = qpts(i,1)*qpts(i,2)^2*(1-qpts(i,1)-qpts(i,2));   %bubble2 x*y^2*(1-x-y)
    end

end

