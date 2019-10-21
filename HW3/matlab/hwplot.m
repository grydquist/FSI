

fileid = fopen('d.txt');
d = fscanf(fileid,'%f');
fclose(fileid);
d1 = d(1:2:end);
d2 = d(2:2:end);
d1 = reshape(d1,nx,ny);
d2 = reshape(d2,nx,ny);

X = reshape(x(1,:),nx,ny);
Y = reshape(x(2,:),nx,ny);

figure
plot(X(:,(nx+1)/2),d2(:,(nx+1)/2),'r')
hold on
plot(Y((nx+1)/2,:),d2((nx+1)/2,:),'--b')
hold on
legend('u^h(x,l_y/2)','u^h(l_x/2,y)')
title('3 X 3 nodes')
xlabel('x')
ylabel('u^h')

% figure
% plot(Y((nx+1)/2,:),d2((nx+1)/2,:),'r')
% hold on
% title('u^h(l_x/2,y), 17X17 nodes')
% xlabel('y')
% ylabel('u^h')