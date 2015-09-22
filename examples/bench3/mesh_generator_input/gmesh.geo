// This is the gmsh file used to create the mesh for this problem.
// The most important thing is that the fracture and  matrix elements
// and injection and production nodes should have physical labels.
// df2d recognizes them through this physical label.

//fracture endpoints
Point(1) = {.18,.4,0,.02}; 
Point(2) = {.75,.7,0,.02};
Point(3) = {.3,.83,0,.02};
Point(4) = {.85,.33,0,.02};
Point(5) = {.55,.74,0,.02};
Point(6) = {.87,.53,0,.02};
Point(7) = {.5,.75,0,.02};
Point(8) = {.4,.16,0,.02};
Point(9) = {.25,.7,0,.02};
Point(10) = {.65,.9,0,.02};
Point(11) = {.35,.3,0,.02};
Point(12) = {.8,.15,0,.02};

//corner points
Point(13) = {0,0,0};
Point(14) = {1,0,0};
Point(15) = {1,1,0};
Point(16) = {0,1,0};

//fracture coincidence points
Point(17) = {.37452,.76226,0,.02};
Point(18) = {.48505,.66178,0,.02};
Point(19) = {.46621,.55040,0,.02};
Point(20) = {.55557,.59769,0,.02};
Point(21) = {.67284,.65939,0,.02};
Point(22) = {.41979,.27674,0,.02};

//lines
Line(1) = {13,14};
Line(2) = {14,15};
Line(3) = {15,16};
Line(4) = {16,13};
Line(5) = {9,17};
Line(6) = {3,17};
Line(7) = {10,17};
Line(8) = {17,18};
Line(9) = {7,18};
Line(10) = {18,19};
Line(11) = {1,19};
Line(12) = {19,20};
Line(13) = {20,21};
Line(14) = {18,20};
Line(15) = {5,21};
Line(16) = {21,2};
Line(17) = {21,6};
Line(18) = {20,4};
Line(19) = {19,22};
Line(20) = {11,22};
Line(21) = {22,12};
Line(22) = {22,8};

//surfaces
Line Loop(23) = {1,2,3,4};
Plane Surface(24) = {23};
Line{5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22} In Surface{24};

//physical tags
//Injection point
Physical Point(300) = {13};
//Production point
Physical Point(301) = {15};
//Matrix surface
Physical Surface(400) = {24};
//Fracture lines
Physical Line(401) = {5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};