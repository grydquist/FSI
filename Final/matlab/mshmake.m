clear IEN x;
lx = 10;
ly = 1;
nx = 5;
ny = 4;
nsd = 2;
x = zeros(nsd,nx*ny);
nt = length(x(1,:));
nsd = 2;

%% Allocate x and conectivity matrix

for i = 1:ny
    x(1,(1 + nx*(i-1):nx*i)) = linspace(0,lx,nx);    
end

tmp = linspace(0,ly,ny);
for i = 1:ny
    x(2,(1 + nx*(i-1):nx*i)) = tmp(i);
end

IEN = delaunay(x(1,:),x(2,:));
eNoN = 3;
nEl = length(IEN(:,1));

save x.mat x
save IEN.mat IEN


id = fopen('x.txt','w');
fprintf(id,'%6.0f %6.0f \n',[nt,nsd]);
fprintf(id,'%6.7f %6.7f \n',x);
fclose(id);

id = fopen('IEN.txt','w');

fprintf(id,'%6.0f %6.0f\n',[nEl,eNoN]);
fprintf(id,'%6.0f %6.0f %6.0f\n',IEN');
fclose(id);