function [eval_phi, grad_phi_x,grad_phi_y] = shape_functions_masslumping(nodes, qpts)

%   [eval_phi, grad_phi_x,grad_phi_y] = shape_functions_masslumping(nodes, qpts)
%   evaluation of shape functions and their gradient with lagrangian nodes
%   in nodes and we evaluate it in qpts 

%   Daniele Ceccarelli & Tommaso Missoni - NAPDE project

if(nargin < 2)
    qpts = nodes;   % if we want to evaluate in nodes -> Mass Lumping
end

%example: reference triangle with 6+1 nodes
%   5
%   |  \
%   |    \
%   6      4
%   |  7     \
%   1----2----4
%
%   For every i in 1:7
%   phi_i(x,y) = a_i1 + a_i2*x + a_i3*y + a_i4*x^2 + a_i5*y^2 +
%                   a_i6*x*y + a_i7*bubble(x,y)
%                        where bubble(x,y) = 27*x*y*(1-x-y)

[n, ~] = size(nodes); 

phi = zeros(n,n); % matrix of coefficients

% We split the evaluation of coeff a_ij and then the evaluation of x,y part
% of the shape function (in S) -> then we multiply S*phi to obtain the
% evaluation of the shape functions on qpts

for i=1:n
    phi(:,i) = coeff_basis(nodes,i); % compute the coeff a_ij (phi(:,i) is the vector a_i:)
end

S = eval_matrix_masslumping(nodes, qpts);     %for example, if d = 2, -> 7 nodes
      %here we have the evaluations of 1,x,y,x^2,y^2,xy,bubble on the nodes

eval_phi = S*phi;   %mult. of S and coefficient -> evaluation of shape functions 
                    %should be diagonal if we choose qpts = nodes!

[grad_phi_x,grad_phi_y] = eval_grad_basis(nodes, qpts, phi);  %eval of grad of basis functions

end
