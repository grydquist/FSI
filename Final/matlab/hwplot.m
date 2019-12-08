

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

ggtmp = zeros(27,27);
ggt2 = ggt;
for i=1:27
    ggtmp(i,:) = ggt((1:27)+(i-1)*27);
end
ggt = ggtmp;
ggt(1:9,1:9) = ggtmp(1:3:end,1:3:end);
ggt(10:18,10:18) = ggtmp(2:3:end,2:3:end);
ggt(19:27,19:27) = ggtmp(3:3:end,3:3:end);
%ggt = reshape(ggt,27,27);
%ggt = ggt';



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
dNM_dUdot_11=alfa_m*(NANB_global+rho_f*tau_SUPS_udotGradNANB_global)+...
    alfa_f*gamma*dt*(rho_f*NAdotGradNBu_global+mu_f*Grad_Na_dot_Grad_Nb_global...
    +mu_f*GradNA_x1GradNB_x1_global+rho_f*tau_SUPS_udotGrad_NA_dot_udotGrad_NB_global...
    +rho_f*nu_LSIC_GradNa_dot_e1_GradNb_dot_e1_global);

dNM_dUdot_11= zeros(9,9)+ ...
    alfa_f*gamma*dt*mu_f*GradNA_x1GradNB_x1_global+...
    alfa_m*(NANB_global)+...
    alfa_f*gamma*dt*rho_f*NAdotGradNBu_global+...
    alfa_f*gamma*dt*mu_f*Grad_Na_dot_Grad_Nb_global+...
    zeros(9,9);


  a=alfa_f*gamma*dt*mu_f*GradNA_x1GradNB_x1_global;
  a2 = a(5,:);
  b=alfa_m*(NANB_global);
  b2 = b(5,:);
  c=alfa_f*gamma*dt*rho_f*NAdotGradNBu_global;
  c2 = c(5,:);
  d=alfa_f*gamma*dt*mu_f*Grad_Na_dot_Grad_Nb_global;
  d2 = d(5,:);

dNM_dUdot_11(5,:);
%% NS
%d1 = d(1:2:end);
%d2 = d(2:2:end);
%d1 = reshape(d1,nx,ny);
%d2 = reshape(d2,nx,ny);
%d = reshape(d,nx,ny,200);
X = reshape(x(1,:),nx,ny);
Y = reshape(x(2,:),nx,ny);
xx = X(1:end,1);
yy = Y(1,1:end)';

%[xx,yy] = meshgrid(X,Y);

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

pts = length(u)/nx/ny;
u = reshape(u,nx,ny,pts);
v = reshape(v,nx,ny,pts);
startx = linspace(X(1),X(end),10);
starty = linspace(Y(1),Y(end),10);

for i=1:pts
    u(:,:,i) = u(:,:,i)';
    v(:,:,i) = v(:,:,i)';
    clf;
    quiver(xx,yy(1:end-1),u(1:end-1,:,i),v(1:end-1,:,i))
    streamline(xx,yy,u(:,:,i),v(:,:,i),startx,starty)
    hold on
    pause(0.1)
end

%% Elasto
fileid = fopen('y.txt');
yt = fscanf(fileid,'%f');
fclose(fileid);
%reshape(u,2,nx,ny,4)
u = yt(1:2:end);
v = yt(2:2:end);

%nx = 9;
%ny = 9;
pts = length(u)/nx/ny;
u = reshape(u,nx,ny,pts);
v = reshape(v,nx,ny,pts);
for i=1:pts
    %u(:,:,i) = u(:,:,i)';
    %v(:,:,i) = v(:,:,i)';
    uu = reshape(u(:,:,i),nx*ny,1);
    vv = reshape(v(:,:,i),nx*ny,1);
    
    clf;
    %quiver(u(1:end-1,:,i),v(1:end-1,:,i))
    %scatter(uu,vv)
    triplot (IEN, uu, vv)
    pbaspect([10,2,1])
    set(gcf, 'Position',  [300, 300, 1200, 400])
    axis([0,10.1,-1,2])
    hold on
    pause(1)
end

%% Plot v field

%% Export solution ----------------------need work
% make sure you have right x and ien
fileid = fopen('u.txt');
ut = fscanf(fileid,'%f');
fclose(fileid);

u = ut(1:2:end);
v = ut(2:2:end);
eNo = 3;
%nNo = 1451;% 21,21,6,3
pts = length(u)/nNo;
u = reshape(u,nNo,pts);
v = reshape(v,nNo,pts);
for ii = 1:pts
    name = horzcat('vtks/result_',num2str(ii),'.vtk');
    fileID = fopen(name,'w');
    fprintf(fileID,'# vtk DataFile Version 3.0 \n');
    fprintf(fileID,'Results \n');
    fprintf(fileID,'ASCII \n');
    fprintf(fileID,'DATASET UNSTRUCTURED_GRID \n');
    fprintf(fileID,'POINTS \t%8.0f float \n',nNo);
    tmp = 0;
    for i = 1:nNo % write x
        fprintf(fileID,'%8.4f\t %8.4f\t',X(i,:), tmp);
        fprintf(fileID,'\n');
    end
    fprintf(fileID,'\n');
    fprintf(fileID,'CELLS \t%8.0f \t%8.0f \n',nEl, 4*nEl);
    for i = 1:nEl   % write IEN
        fprintf(fileID,'%8.0f\t %8.0f\t',eNon, IEN(i,:)-1);
        fprintf(fileID,'\n');
    end
    fprintf(fileID,'\n');
    fprintf(fileID,'CELL_TYPES \t%8.0f \n',nEl);
    for i = 1:nEl   % write IEN
        fprintf(fileID,'5 \n');
    end
    fprintf(fileID,'\n');
    fprintf(fileID,'POINT_DATA \t%8.0f \n',nNo);
    fprintf(fileID,'VECTORS U float \n');
    for i = 1:nNo % write u
        fprintf(fileID,'%8.4f\t %8.4f\t %8.4f\t',u(i,ii),v(i,ii), tmp);
        fprintf(fileID,'\n');
    end
    fclose(fileID);
end


disp('ran')
