function [V,H,D] = quadbraz(X,K)
################################################
## Evaluates the properties of a quad element
## using the brazilian papers formulation,
################################################

del = [.5, 0, -.5 0; 0, .5 0 -.5];
fip = transpose( [ .25,.5  ; .5,.25  ; .75,.5  ; .5,.75  ] ) ;
vip = transpose( [ .25,.25 ; .75,.25 ; .75,.75 ; .25,.75 ] ) ;
r = [ 0, -1; 1, 0 ];

b = sderi([.5,.5]);
D = transpose( b * inverse(X * b) );

H = zeros(4, 4);
V = zeros(1,4);
for i = 1 : 4
  b = sderi( fip(:,i) );
  H(i,:) = transpose (b * inverse( X * b ) * transpose(K) * r * X *  b * del(:,i) );
  b = sderi ( vip(:,i) );
  V(i) = det (X * b) / 4;
endfor
endfunction

function B = sderi (Z)
###############################################
## shape function derivative matrix
###############################################

zeta = Z(1);
eta = Z(2);
B = [ -(1-eta) -(1-zeta) ; 1-eta -zeta ; eta zeta ; -eta 1-zeta ];

endfunction
