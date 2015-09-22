//points
Point(1) = {-.8,0,0};
Point(2) = {0,0,0};
Point(3) = {.4,0,0};
Point(6) = {-.8,1,0};
Point(5) = {0,1,0};
Point(4) = {.4,1,0};

//lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {6,5};
Line(4) = {5,4};
Line(5) = {1,6}; //inject
Line(6) = {2,5};
Line(7) = {3,4}; //production

//surfaces
Line Loop(8) = {1,6,-3,-5};
Line Loop(9) = {2,7,-4,-6};
Plane Surface(10) = {8}; //left mat
Plane Surface(11) = {9}; //right mat

//physical stuff
Physical Line(100) = {5};
Physical Line(101) = {7};
Physical Surface(200) = {10};
Physical Surface(201) = {11};