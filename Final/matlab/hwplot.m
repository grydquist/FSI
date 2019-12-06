

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

pts = 25;
u = reshape(u,nx,ny,pts);
v = reshape(v,nx,ny,pts);
for i=1:pts
    %u(:,:,i) = u(:,:,i)';
    %v(:,:,i) = v(:,:,i)';
    clf;
    quiver(u(:,:,i),v(:,:,i))
    hold on
    pause(0.1)
end