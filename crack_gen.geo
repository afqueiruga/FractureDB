//SetFactory("OpenCASCADE");



a = 0.1;
clcrack = a/5.0;

// Inputs
x1 = -a;
y1 = -a;
x2 =  a;
y2 =  a;

//clbig = clcrack;
clbig = 0.5;

// Rectangle(1) = {-1, -1, 0, 1, 1, 0};
Point(1) = {-1,-1,0,clbig};
Point(2) = { 1,-1,0,clbig};
Point(3) = { 1, 0,0,clbig};
Point(4) = { 1, 1,0,clbig};
Point(5) = {-1, 1,0,clbig};
Point(6) = {-1, 0,0,clcrack};

Point(7) = { x1, y1,0,clcrack};
Point(8) = { x2, y2,0,clcrack};


//  5<-------4
//  v    ,>8 ^
//  6->7'    3
//  v        ^
//  1------->2

Line(1) = {1,2}; // beginning the boundary loop
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1}; // closes the boundary

// Line(7) = {6,7}; // the crack starts here
// Line(8) = {7,8};
Spline(7) = {6, 7, 8};
Line(9) = {8,3}; // Not part of the crack, closes the bisection

// For the mesher, the body is split in two,
// but the crack is not like that
Line Loop(1) = {6,1,2, -9,-7};
Line Loop(2) = {3,4,5,  7,9};

Plane Surface(1) = { 1 };
Plane Surface(2) = { 2 };

Physical Surface(1) = {1,2};
Physical Line(2) = { 7,8 };
Physical Line(3) = { 1,2,6,3,4,5 };

Mesh 2;

Plugin(Crack).Dimension = 1;
Plugin(Crack).PhysicalGroup = 2;
Plugin(Crack).OpenBoundaryPhysicalGroup = 0;
Plugin(Crack).Run ;

Physical Line(4) = { 10 };

Save "crack2.msh";
//+
