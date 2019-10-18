function [xg,Wg,Ng] = gausser
%% Get values for Gauss points, shps, etc.

% Gauss points
xg = [1/6,1/6,2/3;1/6,2/3,1/6]';
Wg = [1/6,1/6,1/6]';

N1 = @(x) x(1);
N2 = @(x) x(2);
N3 = @(x) 1 - x(1) - x(2);
Ng = zeros(3,length(xg(:,1)));

% Evaluate shps at Gauss fns
for j = 1:3
    Ng(1,j) = N1(xg(j,:));
    Ng(2,j) = N2(xg(j,:));
    Ng(3,j) = N3(xg(j,:));
end

end