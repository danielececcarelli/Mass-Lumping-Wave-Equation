function [grad_phi_x, grad_phi_y] = eval_grad_basis(nodes, qpts, phi)

%   [grad_phi_x, grad_phi_y] = eval_grad_basis(nodes, qpts, phi)
%   evaluation of Gradient of basis functions on qpts
%   we have as imput the matrix phi of coeff (a_ij)
%   the nodes of the triangle and the points of quadrature qpts

%   Daniele Ceccarelli & Tommaso Missoni - NAPDE project

if(nargin < 2)
    qpts = nodes;
end

% for example, if d=2 => 7 nodes
% grad phi (x,y) : 
%(a2 + 2*a4*x + a6*y + a7*y*(1-x-y) - a7*x*y)
%(a3 + 2*a5*y + a6*x + a7*x*(1-x-y) - a7*x*y)

[n, ~] = size(nodes);
[q, ~] = size(qpts);

grad_phi_x = zeros(q,n);
grad_phi_y = zeros(q,n);

if(n == 3) % degree d = 1
    for i=1:q
        for j=1:n
            grad_phi_x(i,j) = phi(2,j);
            grad_phi_y(i,j) = phi(3,j);
        end
    end
    
elseif(n == 7) % degree d = 2
    for i=1:q
    xx = qpts(i,1);
    yy = qpts(i,2);
        for j=1:n
            grad_phi_x(i,j) = phi(2,j) + 2*phi(4,j)*xx + phi(6,j)*yy + ...
                phi(7,j)*(yy - yy^2 - 2*xx*yy);
            grad_phi_y(i,j) = phi(3,j) + 2*phi(5,j)*yy + phi(6,j)*xx + ...
                phi(7,j)*(xx - xx^2 - 2*xx*yy);
        end
    end
    
else % n == 12, degree d = 3
    for i=1:q
    xx = qpts(i,1);
    yy = qpts(i,2);
        for j=1:n
            grad_phi_x(i,j) = phi(2,j) + 2*phi(4,j)*xx + phi(6,j)*yy + ...
                3*phi(7,j)*xx^2 + 2*phi(9,j)*xx*yy + phi(10,j)*yy^2 + ...
                phi(11,j)*(-xx^2*yy - 2*xx*yy*(xx+yy-1))+ ...
                phi(12,j)*(-yy^2*(xx+yy-1) - xx*yy^2);
            grad_phi_y(i,j) = phi(3,j) + 2*phi(5,j)*yy + phi(6,j)*xx + ...
                3*phi(8,j)*yy^2 + phi(9,j)*xx^2 + 2*phi(10,j)*xx*yy + ...
                phi(11,j)*(-xx^2*(xx+yy-1) - xx^2*yy)+ ...
                phi(12,j)*(-xx*yy^2 - 2*xx*yy*(xx+yy-1));
        end
    end
end

end

