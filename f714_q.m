n = 64;
nt = 64;

x = zeros((n+1)*(n+1),1);
y = zeros((n+1)*(n+1),1);
for i = 1:(n+1)
 for j = 1:(n+1)
  x((n+1)*(i-1)+j) = (i-1)/n;
  y((n+1)*(i-1)+j) = (j-1)/n;
 end
end
t = zeros(nt+1,1);
for i = 1:(nt+1)
 t(i) = (i-1)/nt;
end

% (ut - u)/dt = delta (ut + u)/2
% ut = u + dt delta (ut + u)/2
% (1 - dt delta/2) ut = (1 + dt delta/2) u, so
% v = (1 + dt delta/2) u;
% ut = (1 - dt delta/2) \ v;

delta = zeros((n+1)*(n+1),(n+1)*(n+1));
for i=2:n
 for j=2:n
  delta((n+1)*(i-1)+j, (n+1)*(i-2)+j)   =    n^2;
  delta((n+1)*(i-1)+j, (n+1)*i    +j)   =    n^2;
  delta((n+1)*(i-1)+j, (n+1)*(i-1)+j-1) =    n^2;
  delta((n+1)*(i-1)+j, (n+1)*(i-1)+j+1) =    n^2;
  delta((n+1)*(i-1)+j, (n+1)*(i-1)+j)   = -4*n^2;
 end
end

u = zeros(nt+1,(n+1)^2);
u(1,:) = (sin(3*pi*x).*sin(pi*y))';
for i = 2:(nt+1)
 v = (eye((n+1)*(n+1)) + delta/(2*nt))*u(i-1,:)';
 u(i,:) = (eye((n+1)*(n+1)) - delta/(2*nt)) \ v;
end

c = zeros(nt+1,(n+1)^2);
c(1,:) = u(1,:);
for i = 2:(nt+1)
 c(i,:) = exp(-10*pi*pi/nt)*c(i-1,:);
end

%u(nt+1,:)
max(abs(u-c),[],'all')