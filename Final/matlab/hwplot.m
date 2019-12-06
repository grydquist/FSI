

fileid = fopen('d.txt');
d = fscanf(fileid,'%f');
fclose(fileid);

fileid = fopen('gg.txt');
gg = fscanf(fileid,'%f');
fclose(fileid);
%%
fileid = fopen('ggt.txt');
ggt = fscanf(fileid,'%f');
fclose(fileid);
ggt = reshape(ggt,27,27);
ggt = ggt';



dm1 = ggt(1:3:end,:);
dm2 = ggt(2:3:end,:);
dc = ggt(3:3:end,:);

dm1d1= dm1(:,1:3:end);
dm1d2= dm1(:,2:3:end);
dm1dp= dm1(:,3:3:end);

dm2d1= dm2(:,1:3:end);
dm2d2= dm2(:,2:3:end);
dm2dp= dm2(:,3:3:end);

dcd1 = dc(:,1:3:end);
dcd2 = dc(:,2:3:end);
dcdp = dc(:,3:3:end);
%%
%d1 = d(1:2:end);
%d2 = d(2:2:end);
%d1 = reshape(d1,nx,ny);
%d2 = reshape(d2,nx,ny);
%d = reshape(d,nx,ny,200);
X = reshape(x(1,:),nx,ny);
Y = reshape(x(2,:),nx,ny);

% for i = 1:200
%     clf
%     contourf(X,Y,d(:,:,i),100,'LineColor','none')
%     pause(0.1)
% end

% figure
% plot(Y((nx+1)/2,:),d2((nx+1)/2,:),'r')
% hold on
% title('u^h(l_x/2,y), 17X17 nodes')
% xlabel('y')
% ylabel('u^h')


fileid = fopen('u.txt');
ut = fscanf(fileid,'%f');
fclose(fileid);
%reshape(u,2,nx,ny,4)
u = ut(1:2:end);
v = ut(2:2:end);

pts = 3;
u = reshape(u,nx,ny,pts);
v = reshape(v,nx,ny,pts);
for i=1:pts
    u(:,:,i) = u(:,:,i)';
    v(:,:,i) = v(:,:,i)';
    clf;
    quiver(u(:,:,i),v(:,:,i))
    hold on
    pause(0.1)
end