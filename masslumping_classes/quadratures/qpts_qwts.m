function [qpts, qwts] = qpts_qwts(degree)

%   [qpts, qwts] = qpts_qwts(degree)
%   Mass Lumping quadrature nodes 

%   Daniele Ceccarelli & Tommaso Missoni - NAPDE project

if(degree == 1)
    qpts = [0,0; 1,0; 0,1];
    qwts = [1/3; 1/3; 1/3];
elseif(degree == 2)
    qpts = [0,0; 0.5,0; 1,0; 0.5,0.5; 0,1; 0, 0.5; 1/3,1/3];
    qwts = [1/40; 1/15; 1/40; 1/15; 1/40; 1/15; 9/40];
elseif(degree == 3)
    alpha = 1/2 - sqrt(441-84*(7-sqrt(7)))/42;
    beta = (1 - 1/sqrt(7))/3;
    qpts = [0,0; alpha,0; 1-alpha,0; 1,0; 1-alpha,alpha; alpha,1-alpha; 0,1; 0,1-alpha; 0,alpha;
        beta,beta; 1-2*beta,beta; beta,1-2*beta];
    w1 = 1/90 - sqrt(7)/720;
    w2 = 7/720 + sqrt(7)/180;
    w3 = 49/360 - 7*sqrt(7)/720;
    qwts = [w1; w2; w2; w1; w2; w2; w1; w2; w2; w3; w3; w3];
else
    error("Error, degree not implemented yet. Choose degree = 1,2,3")
end


end

