function [phi] = coeff_basis(nodes,k)

% Computation of the coefficients for k-basis
% in order to be = 1 on the k-node, 0 in all the other 
% (Lagrangian Basis)

[n, ~] = size(nodes);
phi = zeros(n,1);

% we set the rhs equal to 1 for k, 0 otherwise
b = zeros(n,1);
b(k) = 1;

% Evaluation of shape functions on the nodes
M = eval_matrix_masslumping(nodes); 

% Backslash to compute the coefficients
phi = M\b;

end

