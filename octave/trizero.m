function [V,H] = trizero(X,K)
################################################
## Evaluates the properties of a triangle element
## using my report 0 formulation,
################################################

del = [-1/6., 1/3., -1/6.; -1/6., -1/6., 1/3.];
b = [-1, -1; 1, 0; 0, 1];
x1=X(1,1); x2=X(1,2); x3=X(1,3);
y1=X(2,1); y2=X(2,2); y3=X(2,3);
r = [ 0, -1; 1, 0 ];
j = [x2-x1 x3-x1 ; y2-y1 y3-y1];
H = zeros(3, 3);
for i = 1 : 3
  H(i,:) = transpose (b * inverse(j) * transpose(K) * r * j * del(:,i) );
endfor
V = [ det(j)/6. det(j)/6. det(j)/6.];

endfunction
   
