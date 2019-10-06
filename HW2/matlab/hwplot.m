

fileid = fopen('d.txt');
d = fscanf(fileid,'%f');
d = reshape(d,nx,ny);
fclose(fileid);

X = reshape(x(1,:),nx,ny);
Y = reshape(x(2,:),nx,ny);

figure
plot(X(:,(nx+1)/2),d(:,(nx+1)/2),'--r')
hold on
plot(X(:,(nx+1)/2),fun(X(:,(nx+1)/2),Y(:,(nx+1)/2)))
title('u^h(x,l_y/2), 17X17 nodes, \lambda = 4')
xlabel('x')
ylabel('u^h or f')
legend('u^h','f')

figure
plot(Y((nx+1)/2,:),d((nx+1)/2,:),'--r')
hold on
plot(X(:,(nx+1)/2),fun(X(:,(nx+1)/2),Y(:,(nx+1)/2)))
title('u^h(l_x/2,y), 17X17 nodes, \lambda = 4')
xlabel('y')
ylabel('u^h or f')
legend('u^h','f')