// Definition of a sphere in a sphere
/////////////////////////////////////


// Inner sphere
///////////////

Point(101) = {0.0               , 0.0               , 0.0               };
Point(102) = {0.0 + 10000000000000.0, 0.0               , 0.0               };
Point(103) = {0.0               , 0.0 + 10000000000000.0, 0.0               };
Point(104) = {0.0 - 10000000000000.0, 0.0               , 0.0               };
Point(105) = {0.0               , 0.0 - 10000000000000.0, 0.0               };
Point(106) = {0.0               , 0.0               , 0.0 - 10000000000000.0};
Point(107) = {0.0               , 0.0               , 0.0 + 10000000000000.0};

Circle(101) = {102,101,103};
Circle(102) = {103,101,104};
Circle(103) = {104,101,105};
Circle(104) = {105,101,102};
Circle(105) = {103,101,106};
Circle(106) = {106,101,105};
Circle(107) = {105,101,107};
Circle(108) = {107,101,103};
Circle(109) = {102,101,107};
Circle(110) = {107,101,104};
Circle(111) = {104,101,106};
Circle(112) = {106,101,102};

Line Loop(113) = {  102,  108, -110};
Line Loop(115) = {  110,  103,  107};
Line Loop(117) = { -108, -109,  101};
Line Loop(119) = { -111, -102,  105};
Line Loop(121) = { -105, -112, -101};
Line Loop(123) = { -103,  111,  106};
Line Loop(125) = { -107,  104,  109};
Line Loop(127) = { -104,  112, -106};

Surface(114) = {113};
Surface(116) = {115};
Surface(118) = {117};
Surface(120) = {119};
Surface(122) = {121};
Surface(124) = {123};
Surface(126) = {125};
Surface(128) = {127};

Surface Loop(129) = {128,126,116,114,120,124,122,118};


// Outer sphere
///////////////

Point(201) = {0.0                , 0.0                , 0.0                };
Point(202) = {0.0 + 7.8e+16, 0.0                , 0.0                };
Point(203) = {0.0                , 0.0 + 7.8e+16, 0.0                };
Point(204) = {0.0 - 7.8e+16, 0.0                , 0.0                };
Point(205) = {0.0                , 0.0 - 7.8e+16, 0.0                };
Point(206) = {0.0                , 0.0                , 0.0 - 7.8e+16};
Point(207) = {0.0                , 0.0                , 0.0 + 7.8e+16};

Circle(201) = {202,201,203};
Circle(202) = {203,201,204};
Circle(203) = {204,201,205};
Circle(204) = {205,201,202};
Circle(205) = {203,201,206};
Circle(206) = {206,201,205};
Circle(207) = {205,201,207};
Circle(208) = {207,201,203};
Circle(209) = {202,201,207};
Circle(210) = {207,201,204};
Circle(211) = {204,201,206};
Circle(212) = {206,201,202};

Line Loop(213) = {  202,  208, -210};
Line Loop(215) = {  210,  203,  207};
Line Loop(217) = { -208, -209,  201};
Line Loop(219) = { -211, -202,  205};
Line Loop(221) = { -205, -212, -201};
Line Loop(223) = { -203,  211,  206};
Line Loop(225) = { -207,  204,  209};
Line Loop(227) = { -204,  212, -206};

Surface(214) = {213};
Surface(216) = {215};
Surface(218) = {217};
Surface(220) = {219};
Surface(222) = {221};
Surface(224) = {223};
Surface(226) = {225};
Surface(228) = {227};

Surface Loop(229) = {228,226,216,214,220,224,222,218};


// Outer sphere minus inner sphere
//////////////////////////////////

Volume(130) = {229, 129};


Field[1]   = MathEval;
Field[1].F = "5e-15 * (x*x + y*y + z*z)";

Field[2] = Max;
Field[2].FieldsList = {1};
Background Field = 2;

Mesh.CharacteristicLengthFactor = 1.0;

Mesh.CharacteristicLengthMin = 500000000000.0;
Mesh.CharacteristicLengthMax = 3900000000000000.0;

Mesh.CharacteristicLengthExtendFromBoundary = 2;

