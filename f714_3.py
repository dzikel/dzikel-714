from math import sin, cos, pi, exp

def xderiv(arr, N):
 xa = [cos(pi*i/N) for i in range(N+1)]
 out = []
 for _ in range(N+1): out += [[0]*(N+1)]
 for i in range(N+1):   # row
  for j in range(N+1):  # output col, index into out
   for k in range(N+1): # current index into row for computation
    if j==k:
     if j==0: c = (2*N**2 + 1)/6
     elif j==N: c = -(2*N**2 + 1)/6
     else: c = -xa[j]/(2*(1 - xa[j]**2))
    else:
     c = (-1)**(j+k)/(xa[j]-xa[k])
     c *= (2 if (j==0 or j==N) else 1)/(2 if (k==0 or k==N) else 1)
    out[i][j] += 2*c*arr[i][k] # xa[] from textbook is over [-1, 1] rather than [0, 1], so output is lower
 return out

def yderiv(arr, N):
 xa = [cos(pi*i/N) for i in range(N+1)]
 out = []
 for _ in range(N+1): out += [[0]*(N+1)]
 for j in range(N+1):   # col
  for i in range(N+1):  # output row, index into out
   for k in range(N+1): # current index into col for computation
    if i==k:
     if i==0: c = (2*N**2 + 1)/6
     elif i==N: c = -(2*N**2 + 1)/6
     else: c = -xa[i]/(2*(1 - xa[i]**2))
    else:
     c = (-1)**(i+k)/(xa[i]-xa[k])
     c *= (2 if (i==0 or i==N) else 1)/(2 if (k==0 or k==N) else 1)
    out[i][j] += 2*c*arr[k][j]
 return out

def chebyshev(N, dt, B):
 # get coefficients for derivative from textbook ch.6,
 # keep N at a reasonable value
 xa = [cos(pi*i/N) for i in range(N+1)]
 prev = []
 curr = []
 new = []
 inter = []
 xdv = []
 ydv = []
 for _ in range(N+1):
  prev += [[0]*(N+1)]; curr += [[0]*(N+1)]; new += [[0]*(N+1)]; inter += [[0]*(N+1)]; xdv += [[0]*(N+1)]; ydv += [[0]*(N+1)]
 for i in range(N+1):
  for j in range(N+1):
   x = (1+xa[i])/2
   y = (1+xa[j])/2
   #curr[i][j] = exp(-400*((x-.5)**2+(y-.5)**2))*dt
   curr[i][j] = sin(B*pi*x)*sin(B*pi*y)*dt
 for _ in range(int(.75//dt)):
  xdv = xderiv(curr, N)
  xdv = xderiv(xdv, N)
  ydv = yderiv(curr, N)
  ydv = yderiv(ydv, N)
  for i in range(N+1):
   for j in range(N+1):
    inter[i][j] = curr[i][j] + (dt**2)*(xdv[i][j]+ydv[i][j])
  xdv = xderiv(inter, N)
  xdv = xderiv(xdv, N)
  ydv = yderiv(inter, N)
  ydv = yderiv(ydv, N)
  # new - 2*curr + prev = (dt**2)*(xdv+ydv), so
  # new = (dt**2)*(xdv+ydv) + 2*curr - prev
  for i in range(N+1):
   for j in range(N+1):
    new[i][j] = 2*curr[i][j] - prev[i][j] + (dt**2)*(xdv[i][j]+ydv[i][j])
  prev = curr
  curr = new
  #print(curr)
 return curr