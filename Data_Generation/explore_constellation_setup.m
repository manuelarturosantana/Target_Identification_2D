
% curve = Build_Curve("circular cavity", 500, 0.75);
curve = Build_Curve("elliptic cavity", 500, 0.75);

%% Snowman 
xshift = [0,0];
yshift = [2.6, 5];

[lp,mc] = setup_constellation(curve, xshift, yshift, true)

%% Triangle

xshift = [5, 2.6];
yshift = [0, 2.4];

setup_constellation(curve, xshift, yshift, true)

%% Line offset 

xshift = [2.7, 5];
yshift = [0, 0];

setup_constellation(curve, xshift, yshift, true)

%% Diagonal
xshift = [2.5, 4.8];
yshift = [2.5, 4.8];

setup_constellation(curve, xshift, yshift, true)


%% Two Far

xshift = [10];
yshift = [0];

setup_constellation(curve, xshift, yshift, true)