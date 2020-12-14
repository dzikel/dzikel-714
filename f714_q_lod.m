n = 16;
nt = 32;

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

delta = zeros(n+1,n+1); % just 2nd derivative
for i=2:n
 delta(i,i-1) =    n^2;
 delta(i,i+1) =    n^2;
 delta(i,i)   = -2*n^2;
end

u = zeros(nt+1,(n+1)^2);
u(1,:) = (x.*sin(y) + y.*cos(x))';
for i = 1:(n+1)
 u(1, i) = 0;
 u(1, (n+1)*n + i) = 0;
 u(1, (n+1)*(i-1) + 1) = 0;
 u(1, (n+1)*i) = 0;
end
for i = 2:(nt+1)
 uu = zeros((n+1)^2,1);
 for x = 2:n % boundaries are constant
  a = zeros(n+1,1);
  for j = 1:(n+1)
   a(j) = u(i-1,(n+1)*(x-1)+j);
  end
  v = (eye(n+1) + delta/(2*nt))*a;
  q = (eye(n+1) - delta/(2*nt)) \ v;
  for j = 1:(n+1)
   uu((n+1)*(x-1)+j) = q(j);
  end
 end
 %uu
 %pause();
 for y = 2:n
  a = zeros(n+1,1);
  for j = 1:(n+1)
   a(j) = uu((n+1)*(j-1)+y);
  end
  v = (eye(n+1) + delta/(2*nt))*a;
  q = (eye(n+1) - delta/(2*nt)) \ v;
  for j = 1:(n+1)
   u(i,(n+1)*(j-1)+y) = q(j);
  end
 end
end

c = zeros(nt+1,(n+1)^2);
c(1,:) = u(1,:);
for i = 2:(nt+1)
 c(i,:) = exp(-10*pi*pi/nt)*c(i-1,:);
end

% copy comparison code from f714_q.m if necessary
%u(nt+1,:)
%w = abs(u-c);
%for i = 1:(nt+1)
% for j = 1:((n+1)^2)
%  if w(i,j) > 0.24
%   disp (i*10000000+j);
%  end
% end
%end
%max(abs(u-c),[],'all')
u(:,50)
