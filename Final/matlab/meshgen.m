%clear all
% Number of grid points in x and y direction
nx=13; %21
ny=19; %21

% Number of grid points for the flexible beam (solid)
bnx=6; % This includes the square block 6
bny=3; % 3

% The overall size of the domain
l=19.5;
h=12;

% Size of the square block
ol=1;
oh=ol;

% Size of the flexible beam
bl=4;
bh=0.06;

% Size and number of grid points for the entrance section
el=4.5;
enx=ny;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normally, you don't need to change anything below this line
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of nodes on the square block should be an integer
onx = round((bnx-1)*ol/(ol+bl)) + 1;
bnx = (onx-1)*(1+bl/ol) + 1;

ny = round((ny-bny)/2)*2 + bny;

% Total number of nodes
nNo = (bnx-2)*enx*2 + (2*enx+nx-onx)*ny + 2;
% Node coordinates
X=zeros(nNo,2);

% Boundary nodes list
wall=zeros(2*(ny+onx-1),1);
inlet=zeros(ny-2,1);
top=zeros(nx,1);
bottom=zeros(nx,1);
bar=zeros((bnx-onx)*bny,1);

% Meshing the entrance
a=0;
b=0;
c=0;
for i=1:enx
   for j=2:ny-1
      a = a + 1;
      X(a,1) = (i-1)*el/(enx-1);
      X(a,2) = (j-1)/(ny-1)*(i-1)/(enx-1)*ol ...
             + (j-1)/(ny-1)*(enx-i)/(enx-1)*h...
             + (i-1)/(enx-1)*(h-ol)/2;
      if (i==1) 
        b = b + 1;
        inlet(b) = a;
      end
      if (i==enx)
        c = c+1;
        wall(c) = a;
      end
   end
end

% Meshing the region above and below the block
% Entrance and wake length
wl = 2*el + ol + bl;
b=0;
for i=1:bnx
   for j=1:enx
      a = a + 1;
      X(a,1) = (i-1)/(bnx-1)*(j-1)/(enx-1)*(bl+ol) ...
             + (i-1)/(bnx-1)*(enx-j)/(enx-1)*wl ...
             + (j-1)/(enx-1)*el;
      X(a,2) = (j-1)/(enx-1)*(h-oh)/2;
      a = a + 1;
      X(a,1) = X(a-1,1);
      X(a,2) = h - X(a-1,2);
      if (j==1) 
        b = b + 1;
        top(b) = a;
        bottom(b) = a-1;
      end
      if (j==enx && i <= onx) 
        c = c + 2;
        wall(c-1:c) = [a-1 a];
      end
   end
end

% Meshing the region above and below the beam
d=0;
for i=1:bnx-onx+1
   for j=1:(ny-bny)/2
      a = a + 1;
      X(a,1) = (i-1)/(bnx-onx)*bl + el + ol;
      X(a,2) = (j-1)/(ny-bny)*(ol-bh) + (h-ol)/2;
      a = a + 1;
      X(a,1) = X(a-1,1);
      X(a,2) = h - X(a-1,2);
      if (i == 1) 
        c = c + 2;
        wall(c-1:c) = [a-1 a];
      end
   end
end

% Meshing the beam
for i=1:bnx-onx+1
   for j=1:bny
      a = a + 1;
      X(a,1) = (i-1)/(bnx-onx)*bl + el + ol;
      X(a,2) = (j-1)/(bny - 1)*bh + (h-bh)/2;
      d = d + 1;
      bar(d) = a;
      if (i == 1) 
        c = c + 1;
        wall(c) = a;
      end
   end
end

% Meshing the top and bottom part of the wake
for i=2:enx
   for j=2:(ny-bny)/2 + 1
      a = a + 1;
      X(a,1) = (i-1)*el/(enx-1) + el + ol + bl;
      X(a,2) = (j-1)/(ny-bny)*(enx-i)/(enx-1)*(ol-bh) ...
             + (j-1)/(ny-1)*(i-1)/(enx-1)*h...
             + (enx-i)/(enx-1)*(h-oh)/2;
      a = a + 1;
      X(a,1) = X(a-1,1);
      X(a,2) = h - X(a-1,2);
   end
end

% Meshing the wake
for i=2:enx
   for j=2:bny - 1
      a = a + 1;
      X(a,1) = (i-1)*el/(enx-1) + el + ol + bl;
      X(a,2) = (j-1)/(bny-1)*(enx-i)/(enx-1)*bh ...
             + (j-1)/(ny-1)*(i-1)/(enx-1)*h ...
             + (enx-i)/(enx-1)*(h-bh)/2 ...
             + (i-1)/(enx-1)*h/2*(ny-bny)/(ny-1);
   end
end

% Meshing the exit
c=b;
for i=2:nx-bnx+1
   for j=1:ny
      a = a + 1;
      X(a,1) = (i-1)/(nx-bnx)*(l-wl) + wl;
      X(a,2) = (j-1)/(ny-1)*h;
      if (j==1) 
        b = b + 1;
        bottom(b) = a;
      elseif(j==ny)
        c = c + 1;  
        top(b) = a;
      end
   end
end

% Constructing the connectivity array
IEN = delaunay(X(:,1),X(:,2));

% Removing elements inside the block
tmp=zeros(nNo,1);
tmp(wall)=1;
rml=[];
for e=1:size(IEN,1)
  if (all(tmp(IEN(e,:))==1))
    rml = [rml; e];
  end
end
IEN(rml,:)=[];
nEl = size(IEN,1);

% Plotting the grid
figure(1)
plot(X(1:a,1),X(1:a,2),'.',X(bar,1),X(bar,2),'r.',...
     X(wall,1),X(wall,2),'k.',X(inlet,1),X(inlet,2),'g.',...
     X(top,1),X(top,2),'c.',X(bottom,1),X(bottom,2),'m.')
axis equal

figure(2)
triplot (IEN, X(:,1), X(:,2));

% Writing the results into the disk
fid=fopen('x.txt','W');
fprintf(fid,'%i\n',nNo,2);
for a=1:nNo
  fprintf(fid,'%.16f %.16f\n',X(a,1),X(a,2));
end
fclose(fid);
fid=fopen('IEN.txt','W');
fprintf(fid,'%i\n',nEl,3);
for e=1:nEl
  fprintf(fid,'%i %i %i\n',double(IEN(e,1)),double(IEN(e,2)),double(IEN(e,3)));
end
fclose(fid);
fid=fopen('wall','W');
fprintf(fid,'%i\n',length(wall));
for a=1:length(wall)
  fprintf(fid,'%i\n',double(wall(a)));
end
fclose(fid);
fid=fopen('top','W');
fprintf(fid,'%i\n',length(top));
for a=1:length(top)
  fprintf(fid,'%i\n',double(top(a)));
end
fclose(fid);
fid=fopen('bottom','W');
fprintf(fid,'%i\n',length(bottom));
for a=1:length(bottom)
  fprintf(fid,'%i\n',double(bottom(a)));
end
fclose(fid);
fid=fopen('beam','W');
fprintf(fid,'%i\n',length(bar));
for a=1:length(bar)
  fprintf(fid,'%i\n',double(bar(a)));
end
fclose(fid);
fid=fopen('inlet','W');
fprintf(fid,'%i\n',length(inlet));
for a=1:length(inlet)
  fprintf(fid,'%i\n',double(inlet(a)));
end
fclose(fid);
