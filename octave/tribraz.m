function [V,H,D] = tribraz(X,K)
################################################
## Evaluates the properties of a triangle element
## using the brazilian papers formulation,
################################################

del = [1/3., -1/6., -1/6.; -1/6., 1/3., -1/6.];
b = [-1, -1; 1, 0; 0, 1];
j = X * b ;
r = [ 0, -1; 1, 0 ];

D = transpose( b * inverse(j) );
H = zeros(3, 3);
for i = 1 : 3
  H(i,:) = transpose (b * inverse(j) * transpose(K) * r * j * del(:,i) );
endfor
V = [ det(j)/6. det(j)/6. det(j)/6.];

endfunction
  
