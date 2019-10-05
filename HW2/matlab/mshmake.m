lx = 8;
ly = 8;
nx = 5;
ny = 5;
nsd = 2;
x = zeros(nsd,nx*ny);
ng = 3;

%% Allocate x and conectivity matrix

for i = 1:ny
    x(1,(1 + nx*(i-1):nx*i)) = linspace(0,lx,nx);    
end

tmp = linspace(0,ly,ny);
for i = 1:ny
    x(2,(1 + nx*(i-1):nx*i)) = tmp(i);
end

IEN = delaunay(x(1,:),x(2,:));

save x.mat x
save IEN.mat IEN