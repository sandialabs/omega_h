#include "Omega_h_array_ops.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_compare.hpp"
#include "Omega_h_vtk.hpp"
#include "Omega_h_xml_lite.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

#ifdef OMEGA_H_USE_GMSH
#include <gmsh.h>
#include "Omega_h_element.hpp"
#include "Omega_h_shape.hpp"
#endif  // OMEGA_H_USE_GMSH

using namespace Omega_h;

static const char* GMSH_SQUARE_GEO = R"GMSH(
Point(1) = {0, 0, 0, .03};
Point(2) = {1, 0, 0, .03};
Point(3) = {1, 1, 0, .03};
Point(4) = {0, 1, 0, .03};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};
Mesh 2;
)GMSH";

static const char* GMSH_SQUARE_MSH2 = R"GMSH(
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
29
1 0 0 0
2 1 0 0
3 1 1 0
4 0 1 0
5 0.2499999999994083 0 0
6 0.4999999999986718 0 0
7 0.7499999999993359 0 0
8 1 0.2499999999994083 0
9 1 0.4999999999986718 0
10 1 0.7499999999993359 0
11 0.7500000000004128 1 0
12 0.5000000000013305 1 0
13 0.2500000000006652 1 0
14 0 0.7500000000004128 0
15 0 0.5000000000013305 0
16 0 0.2500000000006652 0
17 0.5000000000000081 0.4999999999999986 0
18 0.7083333333332368 0.7083333333334162 0
19 0.2916666666665839 0.7083333333332364 0
20 0.708333333333436 0.2916666666667733 0
21 0.2916666666667735 0.291666666666564 0
22 0.4999999999998797 0.7797619047619079 0
23 0.2202380952381002 0.4999999999998676 0
24 0.7797619047619139 0.500000000000129 0
25 0.5000000000001404 0.2202380952380942 0
26 0.8385416666668851 0.8385416666665038 0
27 0.1614583333334962 0.8385416666668857 0
28 0.8385416666664522 0.1614583333331029 0
29 0.1614583333331026 0.1614583333335478 0
$EndNodes
$Elements
60
1 15 2 0 1 1
2 15 2 0 2 2
3 15 2 0 3 3
4 15 2 0 4 4
5 1 2 0 1 1 5
6 1 2 0 1 5 6
7 1 2 0 1 6 7
8 1 2 0 1 7 2
9 1 2 0 2 2 8
10 1 2 0 2 8 9
11 1 2 0 2 9 10
12 1 2 0 2 10 3
13 1 2 0 3 3 11
14 1 2 0 3 11 12
15 1 2 0 3 12 13
16 1 2 0 3 13 4
17 1 2 0 4 4 14
18 1 2 0 4 14 15
19 1 2 0 4 15 16
20 1 2 0 4 16 1
21 2 2 0 6 19 14 23
22 2 2 0 6 18 11 22
23 2 2 0 6 21 5 25
24 2 2 0 6 20 8 24
25 2 2 0 6 13 19 22
26 2 2 0 6 16 21 23
27 2 2 0 6 10 18 24
28 2 2 0 6 7 20 25
29 2 2 0 6 14 15 23
30 2 2 0 6 11 12 22
31 2 2 0 6 5 6 25
32 2 2 0 6 8 9 24
33 2 2 0 6 12 13 22
34 2 2 0 6 15 16 23
35 2 2 0 6 9 10 24
36 2 2 0 6 6 7 25
37 2 2 0 6 17 18 22
38 2 2 0 6 19 17 22
39 2 2 0 6 17 19 23
40 2 2 0 6 18 17 24
41 2 2 0 6 21 17 23
42 2 2 0 6 17 20 24
43 2 2 0 6 20 17 25
44 2 2 0 6 17 21 25
45 2 2 0 6 11 18 26
46 2 2 0 6 14 19 27
47 2 2 0 6 8 20 28
48 2 2 0 6 5 21 29
49 2 2 0 6 20 7 28
50 2 2 0 6 21 16 29
51 2 2 0 6 18 10 26
52 2 2 0 6 19 13 27
53 2 2 0 6 16 1 29
54 2 2 0 6 7 2 28
55 2 2 0 6 10 3 26
56 2 2 0 6 13 4 27
57 2 2 0 6 3 11 26
58 2 2 0 6 4 14 27
59 2 2 0 6 2 8 28
60 2 2 0 6 1 5 29
$EndElements
)GMSH";

static const char* GMSH_SQUARE_MSH40 = R"GMSH(
$MeshFormat
4 0 8
$EndMeshFormat
$Entities
4 4 1 0
1 0 0 0 0 0 0 0
2 1 0 0 1 0 0 0
3 1 1 0 1 1 0 0
4 0 1 0 0 1 0 0
1 0 0 0 1 0 0 0 2 1 -2
2 1 0 0 1 1 0 0 2 2 -3
3 0 1 0 1 1 0 0 2 3 -4
4 0 0 0 0 1 0 0 2 4 -1
6 0 0 0 1 1 0 0 4 4 1 2 3
$EndEntities
$Nodes
9 29
1 0 0 1
1 0 0 0
2 0 0 1
2 1 0 0
3 0 0 1
3 1 1 0
4 0 0 1
4 0 1 0
1 1 0 3
5 0.2499999999994083 0 0
6 0.4999999999986718 0 0
7 0.7499999999993359 0 0
2 1 0 3
8 1 0.2499999999994083 0
9 1 0.4999999999986718 0
10 1 0.7499999999993359 0
3 1 0 3
11 0.7500000000004128 1 0
12 0.5000000000013305 1 0
13 0.2500000000006652 1 0
4 1 0 3
14 0 0.7500000000004128 0
15 0 0.5000000000013305 0
16 0 0.2500000000006652 0
6 2 0 13
17 0.5000000000000081 0.4999999999999986 0
18 0.7083333333332368 0.7083333333334162 0
19 0.2916666666665839 0.7083333333332364 0
20 0.708333333333436 0.2916666666667733 0
21 0.2916666666667735 0.291666666666564 0
22 0.4999999999998797 0.7797619047619079 0
23 0.2202380952381002 0.4999999999998676 0
24 0.7797619047619139 0.500000000000129 0
25 0.5000000000001404 0.2202380952380942 0
26 0.8385416666668851 0.8385416666665038 0
27 0.1614583333334962 0.8385416666668857 0
28 0.8385416666664522 0.1614583333331029 0
29 0.1614583333331026 0.1614583333335478 0
$EndNodes
$Elements
9 60
1 0 15 1
1 1
2 0 15 1
2 2
3 0 15 1
3 3
4 0 15 1
4 4
1 1 1 4
5 1 5
6 5 6
7 6 7
8 7 2
2 1 1 4
9 2 8
10 8 9
11 9 10
12 10 3
3 1 1 4
13 3 11
14 11 12
15 12 13
16 13 4
4 1 1 4
17 4 14
18 14 15
19 15 16
20 16 1
6 2 2 40
21 19 14 23
22 18 11 22
23 21 5 25
24 20 8 24
25 13 19 22
26 16 21 23
27 10 18 24
28 7 20 25
29 14 15 23
30 11 12 22
31 5 6 25
32 8 9 24
33 12 13 22
34 15 16 23
35 9 10 24
36 6 7 25
37 17 18 22
38 19 17 22
39 17 19 23
40 18 17 24
41 21 17 23
42 17 20 24
43 20 17 25
44 17 21 25
45 11 18 26
46 14 19 27
47 8 20 28
48 5 21 29
49 20 7 28
50 21 16 29
51 18 10 26
52 19 13 27
53 16 1 29
54 7 2 28
55 10 3 26
56 13 4 27
57 3 11 26
58 4 14 27
59 2 8 28
60 1 5 29
$EndElements
)GMSH";

static const char* GMSH_SQUARE_MSH41 = R"GMSH(
$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
4 4 1 0
1 0 0 0 0
2 1 0 0 0
3 1 1 0 0
4 0 1 0 0
1 0 0 0 1 0 0 0 2 1 -2
2 1 0 0 1 1 0 0 2 2 -3
3 0 1 0 1 1 0 0 2 3 -4
4 0 0 0 0 1 0 0 2 4 -1
6 0 0 0 1 1 0 0 4 4 1 2 3
$EndEntities
$Nodes
9 29 1 29
0 1 0 1
1
0 0 0
0 2 0 1
2
1 0 0
0 3 0 1
3
1 1 0
0 4 0 1
4
0 1 0
1 1 0 3
5
6
7
0.2499999999994083 0 0
0.4999999999986718 0 0
0.7499999999993359 0 0
1 2 0 3
8
9
10
1 0.2499999999994083 0
1 0.4999999999986718 0
1 0.7499999999993359 0
1 3 0 3
11
12
13
0.7500000000004128 1 0
0.5000000000013305 1 0
0.2500000000006652 1 0
1 4 0 3
14
15
16
0 0.7500000000004128 0
0 0.5000000000013305 0
0 0.2500000000006652 0
2 6 0 13
17
18
19
20
21
22
23
24
25
26
27
28
29
0.5000000000000081 0.4999999999999986 0
0.7083333333332368 0.7083333333334162 0
0.2916666666665839 0.7083333333332364 0
0.708333333333436 0.2916666666667733 0
0.2916666666667735 0.291666666666564 0
0.4999999999998797 0.7797619047619079 0
0.2202380952381002 0.4999999999998676 0
0.7797619047619139 0.500000000000129 0
0.5000000000001404 0.2202380952380942 0
0.8385416666668851 0.8385416666665038 0
0.1614583333334962 0.8385416666668857 0
0.8385416666664522 0.1614583333331029 0
0.1614583333331026 0.1614583333335478 0
$EndNodes
$Elements
9 60 1 60
0 1 15 1
1 1
0 2 15 1
2 2
0 3 15 1
3 3
0 4 15 1
4 4
1 1 1 4
5 1 5
6 5 6
7 6 7
8 7 2
1 2 1 4
9 2 8
10 8 9
11 9 10
12 10 3
1 3 1 4
13 3 11
14 11 12
15 12 13
16 13 4
1 4 1 4
17 4 14
18 14 15
19 15 16
20 16 1
2 6 2 40
21 19 14 23
22 18 11 22
23 21 5 25
24 20 8 24
25 13 19 22
26 16 21 23
27 10 18 24
28 7 20 25
29 14 15 23
30 11 12 22
31 5 6 25
32 8 9 24
33 12 13 22
34 15 16 23
35 9 10 24
36 6 7 25
37 17 18 22
38 19 17 22
39 17 19 23
40 18 17 24
41 21 17 23
42 17 20 24
43 20 17 25
44 17 21 25
45 11 18 26
46 14 19 27
47 8 20 28
48 5 21 29
49 20 7 28
50 21 16 29
51 18 10 26
52 19 13 27
53 16 1 29
54 7 2 28
55 10 3 26
56 13 4 27
57 3 11 26
58 4 14 27
59 2 8 28
60 1 5 29
$EndElements
)GMSH";

static const char* GMSH_PHYSICAL_MSH2 = R"GMSH(
$MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
2
2 1 "Left"
2 2 "Right"
$EndPhysicalNames
$Nodes
43
1 0 0 0
2 0.5 0 0
3 1 0 0
4 1 1 0
5 0.5 1 0
6 0 1 0
7 0.1249999999997042 0 0
8 0.2499999999993359 0 0
9 0.3749999999996679 0 0
10 0.6249999999996034 0 0
11 0.7499999999993106 0 0
12 0.8749999999997486 0 0
13 1 0.2499999999994083 0
14 1 0.4999999999986718 0
15 1 0.7499999999993359 0
16 0.8750000000003966 1 0
17 0.7500000000006894 1 0
18 0.6250000000002514 1 0
19 0.3750000000002064 1 0
20 0.2500000000006652 1 0
21 0.1250000000003326 1 0
22 0 0.7500000000004128 0
23 0 0.5000000000013305 0
24 0 0.2500000000006652 0
25 0.5 0.2499999999994083 0
26 0.5 0.4999999999986718 0
27 0.5 0.7499999999993359 0
28 0.229166666666788 0.7317708333334298 0
29 0.2591107871718918 0.276193513119521 0
30 0.2603375960466909 0.5111917059928515 0
31 0.3266723356012474 0.8854078089567969 0
32 0.1551398337115055 0.8945297737151066 0
33 0.1577684645284056 0.1059280855200176 0
34 0.3403132086165502 0.1053535997731578 0
35 0.1260125147443213 0.3899617461021146 0
36 0.1245062878724942 0.6123971410580239 0
37 0.7499999999999261 0.2853043002912767 0
38 0.75 0.75 0
39 0.7499999999999968 0.5232240169151434 0
40 0.6562499999997282 0.107446550048445 0
41 0.8385416666664522 0.1071251417231883 0
42 0.8284438775514187 0.8869579081631257 0
43 0.6589073129253932 0.8978263180270769 0
$EndNodes
$Elements
60
1 2 2 1 10 26 27 30
2 2 2 1 10 29 26 30
3 2 2 1 10 25 26 29
4 2 2 1 10 27 28 30
5 2 2 1 10 22 28 32
6 2 2 1 10 25 29 34
7 2 2 1 10 28 27 31
8 2 2 1 10 29 24 33
9 2 2 1 10 19 20 31
10 2 2 1 10 7 8 33
11 2 2 1 10 20 21 32
12 2 2 1 10 8 9 34
13 2 2 1 10 23 24 35
14 2 2 1 10 28 22 36
15 2 2 1 10 30 28 36
16 2 2 1 10 23 35 36
17 2 2 1 10 35 30 36
18 2 2 1 10 29 30 35
19 2 2 1 10 22 23 36
20 2 2 1 10 24 29 35
21 2 2 1 10 31 20 32
22 2 2 1 10 33 8 34
23 2 2 1 10 28 31 32
24 2 2 1 10 29 33 34
25 2 2 1 10 1 7 33
26 2 2 1 10 24 1 33
27 2 2 1 10 9 2 34
28 2 2 1 10 2 25 34
29 2 2 1 10 5 19 31
30 2 2 1 10 27 5 31
31 2 2 1 10 21 6 32
32 2 2 1 10 6 22 32
33 2 2 2 11 14 37 13
34 2 2 2 11 25 37 26
35 2 2 2 11 14 39 37
36 2 2 2 11 37 39 26
37 2 2 2 11 15 39 14
38 2 2 2 11 26 39 27
39 2 2 2 11 38 39 15
40 2 2 2 11 27 39 38
41 2 2 2 11 37 41 13
42 2 2 2 11 38 43 27
43 2 2 2 11 25 40 37
44 2 2 2 11 15 42 38
45 2 2 2 11 11 40 10
46 2 2 2 11 17 42 16
47 2 2 2 11 12 41 11
48 2 2 2 11 18 43 17
49 2 2 2 11 11 41 40
50 2 2 2 11 17 43 42
51 2 2 2 11 40 41 37
52 2 2 2 11 42 43 38
53 2 2 2 11 10 40 2
54 2 2 2 11 2 40 25
55 2 2 2 11 3 41 12
56 2 2 2 11 13 41 3
57 2 2 2 11 4 42 15
58 2 2 2 11 16 42 4
59 2 2 2 11 5 43 18
60 2 2 2 11 27 43 5
$EndElements
)GMSH";

static const char* GMSH_PHYSICAL_MSH40 = R"GMSH(
$MeshFormat
4 0 8
$EndMeshFormat
$PhysicalNames
2
2 1 "Left"
2 2 "Right"
$EndPhysicalNames
$Entities
6 7 2 0
1 0 0 0 0 0 0 0
2 0.5 0 0 0.5 0 0 0
3 1 0 0 1 0 0 0
4 1 1 0 1 1 0 0
5 0.5 1 0 0.5 1 0 0
6 0 1 0 0 1 0 0
1 0 0 0 0.5 0 0 0 2 1 -2
2 0.5 0 0 1 0 0 0 2 2 -3
3 1 0 0 1 1 0 0 2 3 -4
4 0.5 1 0 1 1 0 0 2 4 -5
5 0 1 0 0.5 1 0 0 2 5 -6
6 0 0 0 0 1 0 0 2 6 -1
7 0.5 0 0 0.5 1 0 0 2 2 -5
10 0 0 0 0.5 1 0 1 1 4 1 7 5 6
11 0.5 0 0 1 1 0 1 2 4 -7 2 3 4
$EndEntities
$Nodes
15 43
1 0 0 1
1 0 0 0
2 0 0 1
2 0.5 0 0
3 0 0 1
3 1 0 0
4 0 0 1
4 1 1 0
5 0 0 1
5 0.5 1 0
6 0 0 1
6 0 1 0
1 1 0 3
7 0.1249999999997042 0 0
8 0.2499999999993359 0 0
9 0.3749999999996679 0 0
2 1 0 3
10 0.6249999999996034 0 0
11 0.7499999999993106 0 0
12 0.8749999999997486 0 0
3 1 0 3
13 1 0.2499999999994083 0
14 1 0.4999999999986718 0
15 1 0.7499999999993359 0
4 1 0 3
16 0.8750000000003966 1 0
17 0.7500000000006894 1 0
18 0.6250000000002514 1 0
5 1 0 3
19 0.3750000000002064 1 0
20 0.2500000000006652 1 0
21 0.1250000000003326 1 0
6 1 0 3
22 0 0.7500000000004128 0
23 0 0.5000000000013305 0
24 0 0.2500000000006652 0
7 1 0 3
25 0.5 0.2499999999994083 0
26 0.5 0.4999999999986718 0
27 0.5 0.7499999999993359 0
10 2 0 9
28 0.229166666666788 0.7317708333334298 0
29 0.2591107871718918 0.276193513119521 0
30 0.2603375960466909 0.5111917059928515 0
31 0.3266723356012474 0.8854078089567969 0
32 0.1551398337115055 0.8945297737151066 0
33 0.1577684645284056 0.1059280855200176 0
34 0.3403132086165502 0.1053535997731578 0
35 0.1260125147443213 0.3899617461021146 0
36 0.1245062878724942 0.6123971410580239 0
11 2 0 7
37 0.7499999999999261 0.2853043002912767 0
38 0.75 0.75 0
39 0.7499999999999968 0.5232240169151434 0
40 0.6562499999997282 0.107446550048445 0
41 0.8385416666664522 0.1071251417231883 0
42 0.8284438775514187 0.8869579081631257 0
43 0.6589073129253932 0.8978263180270769 0
$EndNodes
$Elements
2 60
10 2 2 32
1 26 27 30
2 29 26 30
3 25 26 29
4 27 28 30
5 22 28 32
6 25 29 34
7 28 27 31
8 29 24 33
9 19 20 31
10 7 8 33
11 20 21 32
12 8 9 34
13 23 24 35
14 28 22 36
15 30 28 36
16 23 35 36
17 35 30 36
18 29 30 35
19 22 23 36
20 24 29 35
21 31 20 32
22 33 8 34
23 28 31 32
24 29 33 34
25 1 7 33
26 24 1 33
27 9 2 34
28 2 25 34
29 5 19 31
30 27 5 31
31 21 6 32
32 6 22 32
11 2 2 28
33 14 37 13
34 25 37 26
35 14 39 37
36 37 39 26
37 15 39 14
38 26 39 27
39 38 39 15
40 27 39 38
41 37 41 13
42 38 43 27
43 25 40 37
44 15 42 38
45 11 40 10
46 17 42 16
47 12 41 11
48 18 43 17
49 11 41 40
50 17 43 42
51 40 41 37
52 42 43 38
53 10 40 2
54 2 40 25
55 3 41 12
56 13 41 3
57 4 42 15
58 16 42 4
59 5 43 18
60 27 43 5
$EndElements
)GMSH";

static const char* GMSH_PHYSICAL_MSH41 = R"GMSH(
$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
2
2 1 "Left"
2 2 "Right"
$EndPhysicalNames
$Entities
6 7 2 0
1 0 0 0 0
2 0.5 0 0 0
3 1 0 0 0
4 1 1 0 0
5 0.5 1 0 0
6 0 1 0 0
1 0 0 0 0.5 0 0 0 2 1 -2
2 0.5 0 0 1 0 0 0 2 2 -3
3 1 0 0 1 1 0 0 2 3 -4
4 0.5 1 0 1 1 0 0 2 4 -5
5 0 1 0 0.5 1 0 0 2 5 -6
6 0 0 0 0 1 0 0 2 6 -1
7 0.5 0 0 0.5 1 0 0 2 2 -5
10 0 0 0 0.5 1 0 1 1 4 1 7 5 6
11 0.5 0 0 1 1 0 1 2 4 -7 2 3 4
$EndEntities
$Nodes
15 43 1 43
0 1 0 1
1
0 0 0
0 2 0 1
2
0.5 0 0
0 3 0 1
3
1 0 0
0 4 0 1
4
1 1 0
0 5 0 1
5
0.5 1 0
0 6 0 1
6
0 1 0
1 1 0 3
7
8
9
0.1249999999997042 0 0
0.2499999999993359 0 0
0.3749999999996679 0 0
1 2 0 3
10
11
12
0.6249999999996034 0 0
0.7499999999993106 0 0
0.8749999999997486 0 0
1 3 0 3
13
14
15
1 0.2499999999994083 0
1 0.4999999999986718 0
1 0.7499999999993359 0
1 4 0 3
16
17
18
0.8750000000003966 1 0
0.7500000000006894 1 0
0.6250000000002514 1 0
1 5 0 3
19
20
21
0.3750000000002064 1 0
0.2500000000006652 1 0
0.1250000000003326 1 0
1 6 0 3
22
23
24
0 0.7500000000004128 0
0 0.5000000000013305 0
0 0.2500000000006652 0
1 7 0 3
25
26
27
0.5 0.2499999999994083 0
0.5 0.4999999999986718 0
0.5 0.7499999999993359 0
2 10 0 9
28
29
30
31
32
33
34
35
36
0.229166666666788 0.7317708333334298 0
0.2591107871718918 0.276193513119521 0
0.2603375960466909 0.5111917059928515 0
0.3266723356012474 0.8854078089567969 0
0.1551398337115055 0.8945297737151066 0
0.1577684645284056 0.1059280855200176 0
0.3403132086165502 0.1053535997731578 0
0.1260125147443213 0.3899617461021146 0
0.1245062878724942 0.6123971410580239 0
2 11 0 7
37
38
39
40
41
42
43
0.7499999999999261 0.2853043002912767 0
0.75 0.75 0
0.7499999999999968 0.5232240169151434 0
0.6562499999997282 0.107446550048445 0
0.8385416666664522 0.1071251417231883 0
0.8284438775514187 0.8869579081631257 0
0.6589073129253932 0.8978263180270769 0
$EndNodes
$Elements
2 60 1 60
2 10 2 32
1 26 27 30
2 29 26 30
3 25 26 29
4 27 28 30
5 22 28 32
6 25 29 34
7 28 27 31
8 29 24 33
9 19 20 31
10 7 8 33
11 20 21 32
12 8 9 34
13 23 24 35
14 28 22 36
15 30 28 36
16 23 35 36
17 35 30 36
18 29 30 35
19 22 23 36
20 24 29 35
21 31 20 32
22 33 8 34
23 28 31 32
24 29 33 34
25 1 7 33
26 24 1 33
27 9 2 34
28 2 25 34
29 5 19 31
30 27 5 31
31 21 6 32
32 6 22 32
2 11 2 28
33 14 37 13
34 25 37 26
35 14 39 37
36 37 39 26
37 15 39 14
38 26 39 27
39 38 39 15
40 27 39 38
41 37 41 13
42 38 43 27
43 25 40 37
44 15 42 38
45 11 40 10
46 17 42 16
47 12 41 11
48 18 43 17
49 11 41 40
50 17 43 42
51 40 41 37
52 42 43 38
53 10 40 2
54 2 40 25
55 3 41 12
56 13 41 3
57 4 42 15
58 16 42 4
59 5 43 18
60 27 43 5
$EndElements
)GMSH";

static void test_file_components(bool is_compressed, bool needs_swapping) {
  using namespace binary;
  std::stringstream stream;
  std::string s = "foo";
  LO n = 10;
  I8 a = 2;
  write_value(stream, a, needs_swapping);
  I32 b = 42 * 1000 * 1000;
  write_value(stream, b, needs_swapping);
  I64 c = I64(42) * 1000 * 1000 * 1000;
  write_value(stream, c, needs_swapping);
  Real d = 4.2;
  write_value(stream, d, needs_swapping);
  Read<I8> aa(n, 0, a);
  write_array(stream, aa, is_compressed, needs_swapping);
  Read<I32> ab(n, 0, b);
  write_array(stream, ab, is_compressed, needs_swapping);
  Read<I64> ac(n, 0, c);
  write_array(stream, ac, is_compressed, needs_swapping);
  Read<Real> ad(n, 0, d);
  write_array(stream, ad, is_compressed, needs_swapping);
  write(stream, s, needs_swapping);
  I8 a2;
  read_value(stream, a2, needs_swapping);
  OMEGA_H_CHECK(a == a2);
  I32 b2;
  read_value(stream, b2, needs_swapping);
  OMEGA_H_CHECK(b == b2);
  I64 c2;
  read_value(stream, c2, needs_swapping);
  OMEGA_H_CHECK(c == c2);
  Real d2;
  read_value(stream, d2, needs_swapping);
  OMEGA_H_CHECK(d == d2);
  Read<I8> aa2;
  read_array(stream, aa2, is_compressed, needs_swapping);
  OMEGA_H_CHECK(aa2 == aa);
  Read<I32> ab2;
  read_array(stream, ab2, is_compressed, needs_swapping);
  OMEGA_H_CHECK(ab2 == ab);
  Read<I64> ac2;
  read_array(stream, ac2, is_compressed, needs_swapping);
  OMEGA_H_CHECK(ac2 == ac);
  Read<Real> ad2;
  read_array(stream, ad2, is_compressed, needs_swapping);
  OMEGA_H_CHECK(ad2 == ad);
  std::string s2;
  read(stream, s2, needs_swapping);
  OMEGA_H_CHECK(s == s2);
}

static void test_file_components() {
  test_file_components(false, false);
  test_file_components(false, true);
#ifdef OMEGA_H_USE_ZLIB
  test_file_components(true, false);
  test_file_components(true, true);
#endif
}

static void build_empty_mesh(Mesh* mesh, Int dim) {
  build_from_elems_and_coords(mesh, OMEGA_H_SIMPLEX, dim, LOs({}), Reals({}));
}

static void test_file(Library* lib, Mesh* mesh0) {
  std::stringstream stream;
  binary::write(stream, mesh0);
  Mesh mesh1(lib);
  mesh1.set_comm(lib->self());
  binary::read(stream, &mesh1, binary::latest_version);
  mesh1.set_comm(lib->world());
  auto opts = MeshCompareOpts::init(mesh0, VarCompareOpts::zero_tolerance());
  compare_meshes(mesh0, &mesh1, opts, true, true);
  OMEGA_H_CHECK(*mesh0 == mesh1);
}

static void test_file(Library* lib) {
  {
    auto mesh0 = build_box(lib->world(), OMEGA_H_SIMPLEX, 1., 1., 1., 1, 1, 1);
    test_file(lib, &mesh0);
  }
  {
    Mesh mesh0(lib);
    build_empty_mesh(&mesh0, 3);
    test_file(lib, &mesh0);
  }
}

#ifdef OMEGA_H_USE_GMSH
Omega_h_Comparison light_compare_meshes(Mesh& a, Mesh& b) {
  OMEGA_H_CHECK(a.comm()->size() == b.comm()->size());
  OMEGA_H_CHECK(a.comm()->rank() == b.comm()->rank());
  const auto comm = a.comm();
  const auto should_print = comm->rank() == 0;
  if (a.family() != b.family()) {
    if (should_print) {
      std::clog << "mesh element families differ\n";
    }
    return OMEGA_H_DIFF;
  }
  if (a.dim() != b.dim()) {
    if (should_print) {
      std::clog << "mesh dimensions differ\n";
    }
    return OMEGA_H_DIFF;
  }

  Omega_h_Comparison result = OMEGA_H_SAME;
  for (Int dim = 0; dim <= a.dim(); ++dim) {
    const auto anents = a.nglobal_ents(dim);
    const auto bnents = b.nglobal_ents(dim);
    if (anents != bnents) {
      if (should_print) {
        std::clog << "global " << topological_singular_name(a.family(), dim)
                  << " counts differ (" << anents << " != " << bnents << ")\n";
      }
      result = OMEGA_H_DIFF;
    }
  }
  const auto lo_measures_a = measure_elements_real(&a);
  const auto lo_measures_b = measure_elements_real(&b);
  const auto go_measure_a =
      get_sum(a.comm(), a.owned_array(a.dim(), lo_measures_a, 1));
  const auto go_measure_b =
      get_sum(b.comm(), b.owned_array(b.dim(), lo_measures_b, 1));
  if (std::abs(go_measure_a / go_measure_b - 1) > 1e-6) {
    if (should_print) {
      std::clog << "total measures of mesh differ (" << go_measure_a
                << " != " << go_measure_b << ")\n";
    }
    result = OMEGA_H_DIFF;
  }
  return result;
}

static void test_gmsh_parallel(Library* lib) {
  ::gmsh::initialize();

  {
    std::ofstream oss("square.geo");
    oss << GMSH_SQUARE_GEO;
  }
  ::gmsh::open("square.geo");
  ::gmsh::write("square.msh");
  ::gmsh::clear();

  auto mesh = Omega_h::gmsh::read("square.msh", lib->world());
  Omega_h::gmsh::write_parallel("square_parallel", mesh);

  auto pmesh = Omega_h::gmsh::read_parallel("square_parallel", lib->world());

  ::gmsh::finalize();
  OMEGA_H_CHECK(light_compare_meshes(mesh, pmesh) == OMEGA_H_SAME);
}

#endif  // OMEGA_H_USE_GMSH

static void test_gmsh(Library* lib) {
  {
    const std::vector<const char*> meshes{
        GMSH_SQUARE_MSH2, GMSH_SQUARE_MSH40, GMSH_SQUARE_MSH41};
    for (const auto& msh : meshes) {
      std::istringstream iss(msh);
      const auto& mesh = Omega_h::gmsh::read(iss, lib->world());
      OMEGA_H_CHECK(mesh.dim() == 2);
      OMEGA_H_CHECK(mesh.nelems() == 40);
      OMEGA_H_CHECK(mesh.nedges() == 68);
      OMEGA_H_CHECK(mesh.nverts() == 29);
      OMEGA_H_CHECK(mesh.class_sets.empty());
    }
  }
  {
    const std::vector<const char*> meshes{
        GMSH_PHYSICAL_MSH2, GMSH_PHYSICAL_MSH40, GMSH_PHYSICAL_MSH41};
    for (const auto& msh : meshes) {
      std::istringstream iss(msh);
      const auto& mesh = Omega_h::gmsh::read(iss, lib->world());
      OMEGA_H_CHECK(mesh.dim() == 2);
      OMEGA_H_CHECK(mesh.nelems() == 60);
      OMEGA_H_CHECK(mesh.nedges() == 102);
      OMEGA_H_CHECK(mesh.nverts() == 43);

      OMEGA_H_CHECK(mesh.class_sets.size() == 2);
      Int id{};
      {
        const auto& classes = mesh.class_sets.find("Left");
        OMEGA_H_CHECK(classes != mesh.class_sets.end());
        const auto class_size = classes->second.size();
        OMEGA_H_CHECK(class_size == 32);
        for (const auto& p : classes->second) {
          OMEGA_H_CHECK(p.dim == 2);
          OMEGA_H_CHECK(p.id == id++);
        }
      }
      {
        const auto classes = mesh.class_sets.find("Right");
        OMEGA_H_CHECK(classes != mesh.class_sets.end());
        OMEGA_H_CHECK(classes->second.size() == 28);
        for (const auto& p : classes->second) {
          OMEGA_H_CHECK(p.dim == 2);
          OMEGA_H_CHECK(p.id == id++);
        }
      }
    }
  }
}

static void test_xml() {
  xml_lite::Tag tag;
  OMEGA_H_CHECK(!xml_lite::parse_tag("AQAAAAAAAADABg", &tag));
  OMEGA_H_CHECK(!xml_lite::parse_tag("   <Foo bar=\"qu", &tag));
  OMEGA_H_CHECK(!xml_lite::parse_tag("   <Foo bar=", &tag));
  OMEGA_H_CHECK(xml_lite::parse_tag("   <Foo bar=\"quux\"   >", &tag));
  OMEGA_H_CHECK(tag.elem_name == "Foo");
  OMEGA_H_CHECK(tag.attribs["bar"] == "quux");
  OMEGA_H_CHECK(tag.type == xml_lite::Tag::START);
  OMEGA_H_CHECK(
      xml_lite::parse_tag("   <Elem att=\"val\"  answer=\"42\" />", &tag));
  OMEGA_H_CHECK(tag.elem_name == "Elem");
  OMEGA_H_CHECK(tag.attribs["att"] == "val");
  OMEGA_H_CHECK(tag.attribs["answer"] == "42");
  OMEGA_H_CHECK(tag.type == xml_lite::Tag::SELF_CLOSING);
  OMEGA_H_CHECK(xml_lite::parse_tag("</Foo>", &tag));
  OMEGA_H_CHECK(tag.elem_name == "Foo");
  OMEGA_H_CHECK(tag.type == xml_lite::Tag::END);
}

static void test_read_vtu(Mesh* mesh0) {
  std::stringstream stream;
  vtk::write_vtu(
      stream, mesh0, mesh0->dim(), vtk::get_all_vtk_tags(mesh0, mesh0->dim()));
  Mesh mesh1(mesh0->library());
  vtk::read_vtu(stream, mesh0->comm(), &mesh1);
  auto opts = MeshCompareOpts::init(mesh0, VarCompareOpts::zero_tolerance());
  OMEGA_H_CHECK(
      OMEGA_H_SAME == compare_meshes(mesh0, &mesh1, opts, true, false));
}

static void test_read_vtu(Library* lib) {
  auto mesh0 = build_box(lib->world(), OMEGA_H_SIMPLEX, 1., 1., 1., 1, 1, 1);
  test_read_vtu(&mesh0);
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  OMEGA_H_CHECK(std::string(lib.version()) == OMEGA_H_SEMVER);
  if (lib.world()->size() == 1) {
    test_file_components();
    test_file(&lib);
    test_xml();
    test_read_vtu(&lib);
    test_gmsh(&lib);
  }
#ifdef OMEGA_H_USE_GMSH
  test_gmsh_parallel(&lib);
#endif  // OMEGA_H_USE_GMSH
}
