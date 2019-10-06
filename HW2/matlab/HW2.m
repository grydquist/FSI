load x.mat
load IEN.mat;

%% Get variables all set
lx = max(x(1,:));
ly = max(x(2,:));
nt = length(x(1,:));
nx = sqrt(nt);
ny = nx;
nsd = 2;
ng = 3;
eNoN = 3;


%% Element cnsts for triangular
eNoN = 3;
nEl = length(IEN(:,1));

%% Get values for Gauss points, shps, etc.

[xg,Wg,Ng] = gausser;

% Now we need dx/dxi/J
% Note that this is the same for all elements here, but isn't in general
dxdxi = zeros(nEl,nsd,nsd);
xa = zeros(nsd,eNoN);
J = zeros(nEl,1);
tmp = zeros(nsd,nsd);
for i = 1:nEl
    xa = x(:,IEN(i,:));
%   dx/dxi
    dxdxi(i,1,1) = xa(1,1) - xa(1,3);
%   dx/deta
    dxdxi(i,1,2) = xa(1,2) - xa(1,3);
%   dy/dxi
    dxdxi(i,2,1) = xa(2,1) - xa(2,3);
%   dy/deta
    dxdxi(i,2,2) = xa(2,2) - xa(2,3);
    tmp(:,:) = dxdxi(i,:,:);
    J(i) = det(tmp);
end

% Check for boundary nodes
bnd = zeros(nt,1);
for i = 1:nt
    if x(1,i) ==0
        bnd(i) = 2;
    end
    if (x(2,i) == 0||x(2,i) == ly)
        bnd(i) = 1;
    end
end

%% Actual function/evaluation of element matrix
lam = 4;
fun = @(x,y) cos(lam*pi*(x/lx+y/ly));
kab = zeros(eNoN,eNoN);
KAB = zeros(nt,nt);
KABt = KAB;
f = zeros(eNoN,1);
F = zeros(nt,1);
Ft = F;
reasm = zeros(eNoN,eNoN,nsd);

xgt = zeros(nsd,1);

for i = 1:nEl
%   Get element coordinates
    xa = x(:,IEN(i,:));
    [kab,f] = localk(bnd(IEN(i,:)),xa,eNoN,nsd,Ng,J(i),Wg,ng,fun);
   
    IENt = IEN(i,:);
%   Put back into global

    [KABt,Ft] = globalk(kab,f,nt,IENt);
    KAB = KAB + KABt;
    F = F + Ft;
    
    kab(:) = 0;
    f(:) = 0;
end

for i = 1:nt
    if bnd(i) ==2
        KAB(i,i) = 1;
        F(i) = 1;
    end
    if bnd(i) ==1
        KAB(i,i) = 1;
        F(i) = 0;
    end
end
    

% Here's the solution!
d = KAB\F;

%% Some post-processing and output 
d = reshape(d,nx,ny);
X = reshape(x(1,:),nx,ny);
Y = reshape(x(2,:),nx,ny);
surf(X,Y,d);
figure
plot(X(:,(nx+1)/2),d(:,(nx+1)/2))
hold on
plot(X(:,(nx+1)/2),fun(X(:,(nx+1)/2),Y(:,(nx+1)/2)))

figure
plot(Y((nx+1)/2,:),d((nx+1)/2,:))
hold on
plot(X(:,(nx+1)/2),fun(X(:,(nx+1)/2),Y(:,(nx+1)/2)))

id = fopen('x.txt','w');
fprintf(id,'%6.0f %6.0f \n',[nt,nsd]);
fprintf(id,'%6.7f %6.7f \n',x);
fclose(id);

id = fopen('IEN.txt','w');

fprintf(id,'%6.0f %6.0f\n',[nEl,eNoN]);
fprintf(id,'%6.0f %6.0f %6.0f\n',IEN');
fclose(id);

