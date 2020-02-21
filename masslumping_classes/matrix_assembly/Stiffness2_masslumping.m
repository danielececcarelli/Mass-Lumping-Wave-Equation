function K=Stiffness2_masslumping(T,fnk)

%   K = Stiffness2_masslumping(T,fnk)
%   Generate stiffness matrix in our mass lumping space 
%   T is our mass lumping mesh of degree d
%   We use our quadrature rule [qpts, qwts] = qpts_qwts(d);
%   and shape_functions_masslumping for shape functions

%   Daniele Ceccarelli & Tommaso Missoni - NAPDE project

if nargin<2 || isempty(fnk)
   fnk=1.0;
end
if isnumeric(fnk)
   nkflag=1;
else
   nkflag=0;
end

% Get the number of free nodes and initialize the stiffness
% matrix to zero:

Nf=length(T.FNodePtrs);
K=sparse(Nf,Nf);

% Get the number of triangles:

Nt=size(T.Elements,1);

% d is the degree of the elements 

d=T.Degree;

% id is the number of nodes per triangle.

id=round((d+2)*(d+1)/2 + d-1);

% Create the reference triangle and the quadrature weights and nodes
% on it.

TR=RefTri_masslumping(d);
[qpts, qwts] = qpts_qwts(d);    %quadrature nodes and weights for masslumping of degree d
npts = length(qwts);    %number of nodes

% Evaluate the gradients of all the basis functions at all
% of the quadrature nodes:

% we can use qpts or TR.Nodes, they are equal
[~,Vs,Vt]= shape_functions_masslumping(qpts, qpts);   %shape_funct for mass lumping

% Add the contributions from each element

for i=1:Nt

   % Get the coordinates and pointers of the nodes:

   [coords,ll]=getNodes_masslumping(T,i);

   % Extract the coordinates of the vertices of the triangle:

   c=coords(1:d:2*d+1,1:2);

   % Transform the triangle to the reference triangle:
   % (The object trans is a struct that describes the
   % transformation (matrix J, etc.).  See TransToRefTri
   % for details.)

   trans=TransToRefTri(c);

   % Transform the gradients to the reference triangle:

   Grads1=zeros(2,npts,id);
   for ii=1:id
      Grads1(:,:,ii)=trans.J'\[Vs(:,ii)';Vt(:,ii)'];
   end

   % Compute all the quadrature nodes on T:

   z=trans.z1*ones(1,npts)+trans.J*qpts';

   % Compute quantities common to the integrals:

   if nkflag
      scale=fnk*trans.j;
      ghat=qwts;
   else
      scale=trans.j;
      ghat=feval(fnk,z(1,:)',z(2,:)').*qwts;
   end

   % Loop over all possible combinations of (global) indices related
   % to this triangle.

   for r=1:id
      llr=ll(r);
      if llr>0
         for s=r:id
            lls=ll(s);
            if lls>0

               % Estimate the integral:

               I=scale*((Grads1(1,:,r).*Grads1(1,:,s)+...
                         Grads1(2,:,r).*Grads1(2,:,s))*ghat);
               ii=min(llr,lls);
               jj=max(llr,lls);
               K(ii,jj)=K(ii,jj)+I;
            end
         end
      end
   end

end

% Fill in the lower triangle of K.

K=K + triu(K,1)';
