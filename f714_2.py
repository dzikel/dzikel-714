# three-point 2nd derivative:
# 1/h^2 -2/h^2 1/h^2
# five-point Laplacian:
#     0  1/h^2     0
# 1/h^2 -4/h^2 1/h^2
#     0  1/h^2     0

from math import exp

def solver(n, nt): # set n according to problem b
    current = []
    for _ in range(nt+1):
        current_add = []
        for __ in range(n+1):
            current_add += [[0]*(n+1)]
        current += [current_add]
        #print("finished _ =", _)
    for x in range(n+1):
        for y in range(n+1):
            current[1][x][y] = exp(-400*((x/n-.5)**2 + (y/n-.5)**2))/nt
    for t in range(1, nt):
        # compute laplacian
        lapl = []
        for _ in range(n+1): lapl += [[0]*(n+1)]
        for x in range(n+1):
            for y in range(n+1):
                lapl_add = -4*current[t][x][y] + (current[t][x-1][y] if x>0 else 0) + \
                    (current[t][x+1][y] if x<n else 0) + (current[t][x][y-1] if y>0 else 0) + \
                    (current[t][x][y+1] if y<n else 0)
                # n^2 (c[t+1] - 2c[t] + c[t-1]) = nt^2 lapl_add, so
                # c[t+1] = (nt/n)^2 lapl_add + 2c[t] - c[t-1]
                current[t+1][x][y] = lapl_add*(nt**2)/(n**2) + 2*current[t][x][y] - current[t-1][x][y]
        #print("finished t =", t)
    return current