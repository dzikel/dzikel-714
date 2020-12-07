from math import sin, pi

def solver(n, nt, ntfin): # ntfin = total time * nt
 out = []
 for _ in range(ntfin+1): out += [[0]*n]
 # u(x,0) = 3/2 + sin(2*pi*x)
 # u_t + u u_x = 0
 for x in range(n): out[0][x] = 1.5 + sin(2*pi*x/n)
 for t in range(ntfin):
  for x in range(n):
   # upwind in space
   if out[t][x] >= 0:
    delta_x = out[t][(x+1)%n] - out[t][x]
   else:
    delta_x = out[t][x] - out[t][(x-1)%n]
   # (out[t+1][x] - out[t][x])*nt = -out[t][x]*delta_x*n, so
   # out[t+1][x] = out[t][x] - out[t][x]*delta_x*n/nt
   out[t+1][x] = out[t][x]*(1 - delta_x*n/nt)
 return out

def gsolver(n, nt, ntfin):
 out = []
 for _ in range(ntfin+1): out += [[0]*n]
 flux = [0]*n
 # integral of f_t + (ff/2)_x = 0, so
 # integral across x of f(x, t + delta t) - integral across x of f(x, t) +
 # flux integral across t of f(x + delta x / 2, t) - flux integral across t of f(x - delta x / 2, t) = 0
 # avg. at t+dt = avg. at t + (dt/dx) * avg. flux at x-dx/2 - (dt/dx) * avg. flux at x+dx/2
 for x in range(n): out[0][x] = 1.5 + sin(2*pi*x/n)
 for t in range(ntfin):
  for x in range(n):
   l = out[t][x]
   r = out[t][(x+1)%n]
   if l >= 0 and r >= 0: flux[x] = l**2 / 2
   elif l < 0 and r < 0: flux[x] = r**2 / 2
   elif l > 0 and r < 0:
    if l > -r: flux[x] = l**2 / 2
    elif l < -r: flux[x] = r**2 / 2
    else: flux[x] = 0
   else: # l < 0, r > 0
    # 0 is -l/(r-l) from l to r, so
    # value is always (rl-rl)/(r-l) = 0
    flux[x] = 0
  for x in range(n):
   # flux[i] stores flux at i + dx/2
   out[t+1][x] = out[t][x] + (n/nt)*flux[(x-1)%n] - (n/nt)*flux[x]
 return out