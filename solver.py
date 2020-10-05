from math import pi, cos

sgn = lambda x: 2*(x >= 0) - 1

def solver(n, iters=10000):
    # jacobi, h=1/n
    current = []
    for _ in range(n+1): current += [[0]*(n+1)]
    next = []
    for _ in range(n+1): next += [[0]*(n+1)]
    for _ in range(iters):
        # handle dirichlet
        for y in range(n+1):
            next[0][y] = cos(2*pi*y/n)
            #next[0][y] = sgn(cos(2*pi*y/n))
            next[n][y] = 0
        # handle neumann, don't redo dirichlet
        for x in range(1, n):
            next[x][0] = current[x][1]
            next[x][n] = current[x][n-1]
        # handle laplacian
        for x in range(1, n):
            for y in range(1, n):
                next[x][y] = (current[x-1][y] + current[x+1][y] + current[x][y-1] + current[x][y+1])/4
        current = next
    return current

def solver_test(n, max_err = 0.1, best = [], adjust_best=True):
    if not best: best = solver(256)
    if adjust_best:
        best = best[::256//n]
        best = [v[::256//n] for v in best]
    # jacobi, h=1/n
    current = []
    for _ in range(n+1): current += [[0]*(n+1)]
    next = []
    for _ in range(n+1): next += [[0]*(n+1)]
    err = 1000000
    iters = 0
    while err > max_err:
        # handle dirichlet
        for y in range(n+1):
            next[0][y] = cos(2*pi*y/n)
            #next[0][y] = sgn(cos(2*pi*y/n))
            next[n][y] = 0
        # handle neumann, don't redo dirichlet
        for x in range(1, n):
            next[x][0] = current[x][1]
            next[x][n] = current[x][n-1]
        # handle laplacian
        for x in range(1, n):
            for y in range(1, n):
                next[x][y] = (current[x-1][y] + current[x+1][y] + current[x][y-1] + current[x][y+1])/4
        # compute current error:
        err = 0
        for x in range(n+1):
            for y in range(n+1):
                err = max(err, abs(next[x][y]-best[x][y]))
        current = next
        iters += 1
    return (iters, current)

def solver_nested(iters=2000):
    n = 4
    current = [[0]*5,[0]*5,[0]*5,[0]*5,[0]*5]
    while True:
        next = []
        for _ in range(n+1): next += [[0]*(n+1)]
        for _ in range(iters):
            # handle dirichlet
            for y in range(n+1):
                next[0][y] = cos(2*pi*y/n)
                #next[0][y] = sgn(cos(2*pi*y/n))
                next[n][y] = 0
            # handle neumann, don't redo dirichlet
            for x in range(1, n):
                next[x][0] = current[x][1]
                next[x][n] = current[x][n-1]
            # handle laplacian
            for x in range(1, n):
                for y in range(1, n):
                    next[x][y] = (current[x-1][y] + current[x+1][y] + current[x][y-1] + current[x][y+1])/4
            current = next
        # upscale
        n *= 2
        if n==512: break # don't upscale if at target size
        next = []
        for _ in range(n+1): next += [[0]*(n+1)]
        for x in range(n//2+1):
            for y in range(n//2+1):
                next[2*x][2*y] = current[x][y]
        for x in range(n//2+1):
            for y in range(n//2):
                next[2*x][2*y+1] = (current[x][y]+current[x][y+1])/2
        for x in range(n//2):
            for y in range(n//2+1):
                next[2*x+1][2*y] = (current[x][y]+current[x+1][y])/2
        for x in range(n//2):
            for y in range(n//2):
                next[2*x+1][2*y+1] = (current[x][y]+current[x][y+1]+current[x+1][y]+current[x+1][y+1])/4
        current = next
    return current
