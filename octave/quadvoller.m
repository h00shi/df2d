function [V,H,D] = quadvoller(X,K)
################################################
## Evaluates the properties of a quad element
## using the brazilian papers formulation,
################################################

##usefull points
one=   [ X(1,1) , X(2,1) ];
two=   [ X(1,2) , X(2,2) ];
three= [ X(1,3) , X(2,3) ];
four=  [ X(1,4) , X(2,4) ];
mone =  (four + one  ) / 2.;
mtwo=   (one + two   ) / 2.;
mthree= (two + three ) / 2.;
mfour = (three + four) / 2.;
cent=   (one + two + three + four) / 4. ;
f = zeros(4,2);
f(1,:) = [.5 0]; # zeta , eta
f(2,:) = [0 .5]; # zeta , eta
f(3,:) = [-.5 0]; # zeta , eta
f(4,:) = [0 -.5]; # zeta , eta

##volume matrix
V(1) = vol(one, mtwo, cent) + vol(one , cent, mone);
V(2) = vol(mtwo, two, mthree) + vol(mtwo, mthree, cent);
V(3) = vol(cent, mthree, three) + vol(cent, three, mfour);
V(4) = vol(cent, mfour, mone) + vol(mone, mfour, four);

##area matrix
dx1 =  ( cent - mone) (1);
dx2 =  ( cent - mtwo) (1);
dx3 =  ( cent - mthree) (1);
dx4 =  (cent - mfour) (1);
dy1 =  ( cent - mone) (2);
dy2 =  ( cent - mtwo) (2);
dy3 =  ( cent - mthree) (2);
dy4 =  (cent - mfour) (2);
da = [-dy1 dx1; -dy2 dx2; -dy3 dx3; -dy4 dx4]; 

##derivative 
D = deri(X,[0,0]);

## H matrix
H = zeros(4,4);
for i = 1 : 4
	##deri(X,f(i,:)) ## print derivative matrix
	H(i,:) = da(i,:) * K * deri(X,f(i,:));
endfor

endfunction

function D = deri(X,Z)
###################################################
## finds the derivative matrix at a certain point
###################################################

zeta = Z(1);
eta =  Z(2);

j = .25 * [ 1+eta -1-eta -1+eta 1-eta ; 1+zeta 1-zeta -1+zeta -1-zeta ] * transpose( X );
jinv = inverse( j );

zx = jinv(1,1);
ex = jinv(1,2);
zy = jinv(2,1);
ey = jinv(2,2);

n1x= .25 * ( +zx*(1+eta) + (1+zeta)*ex ) ; 
n2x= .25 * ( -zx*(1+eta) + (1-zeta)*ex ) ;
n3x= .25 * ( -zx*(1-eta) - (1-zeta)*ex ) ;
n4x= .25 * ( +zx*(1-eta) - (1+zeta)*ex ) ;

n1y= .25 * ( +zy*(1+eta) + (1+zeta)*ey ) ; 
n2y= .25 * ( -zy*(1+eta) + (1-zeta)*ey ) ;
n3y= .25 * ( -zy*(1-eta) - (1-zeta)*ey ) ;
n4y= .25 * ( +zy*(1-eta) - (1+zeta)*ey ) ;

D = [n1x n2x n3x n4x ; n1y n2y n3y n4y];

endfunction

function ans = vol(r1,r2,r3)
###################################################
## finds the volume of a triangle
###################################################

ans = .5 * det ( [1 r1(1) r1(2); 1 r2(1) r2(2); 1 r3(1) r3(2)] );

endfunction
