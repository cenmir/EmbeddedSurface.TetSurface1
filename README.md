# EmbeddedSurface.TetSurface1
Linear (P1) Tetrahedral embedded Surface repository. Includes extraction of discrete surfaces from an implicit surfacae and visualization.

## Contents
**TetSurface1** Tetrahedral embedded P1 Surface

    o = TetSurface1(H, phi, level)

H is a Mesh.Tet1 class.
phi is a function_handle describing the implicit function
phi is later redifined to be a scalar distance function
level is the iso-level

o contains the surface object with the following properties and methods:

### Properties

- Surface 		- A struct array containing the extracted surfaces
- Count 		- Number of extracted surfaces
- CutElements 	- Vector containing indices to the cut elements
- nCutEle 		- Total number of cut elements

### Methods

- TetSurface1() - The constructor
- vizP1Surf() 

## Demo
```matlab
%% Mesh
T = Mesh.Tet1(x0,x1,nxe,y0,y1,nye,z0,z1,nze,'fishbone');

%% Surface function
xc = mean([x0,x1]); yc = mean([y0,y1]); zc = mean([z0,z1]);
surfaceFunction = @(x,y,z) ((x-xc).^2+(z - zc).^2+(y - yc).^2).^.5-R;

%% Discrete surface
xnod = T.XC;ynod = T.YC;znod = T.ZC;
phi = surfaceFunction(xnod,ynod,znod);

%% Create Surface
level = 0;
SO = EmbeddedSurface.TetSurface1(T, surfaceFunction, level);
SO.vizP1Surf()
```

