function [V,H,D] = trivoller(X,K)
################################################
## Evaluates the properties of a triangle element
## using Voller formulation,
################################################

## the area matrix
dx1 = -X(1,1)/6. +X(1,2)/3. -X(1,3)/6.;
dy1 = -X(2,1)/6. +X(2,2)/3. -X(2,3)/6.;
dx2 = -X(1,1)/6. -X(1,2)/6. +X(1,3)/3.;
dy2 = -X(2,1)/6. -X(2,2)/6. +X(2,3)/3.;
dx3 = +X(1,1)/3. -X(1,2)/6. -X(1,3)/6.;
dy3 = +X(2,1)/3. -X(2,2)/6. -X(2,3)/6.;
da = [-dy1 dx1; -dy2 dx2; -dy3 dx3];

## volume
v = .5 * det( [1 X(1,1) X(2,1); 1 X(1,2) X(2,2); 1 X(1,3) X(2,3) ] );

## derivative
n1x = ( X(2,2) - X(2,3) ) / 2 / v;
n2x = ( X(2,3) - X(2,1) ) / 2 / v;
n3x = ( X(2,1) - X(2,2) ) / 2 / v;
n1y = ( X(1,3) - X(1,2) ) / 2 / v;
n2y = ( X(1,1) - X(1,3) ) / 2 / v;
n3y = ( X(1,2) - X(1,1) ) / 2 / v;


D = [ n1x n2x n3x ; n1y n2y n3y];
V = [ v/3. v/3. v/3.];
H = da * K * D ;

endfunction
  
