import numpy as np
import matplotlib.pyplot as plt

#=======================================================================================
# Input parameters

# Geometrical properties
lx = 1.0  # Width (in meters)
ly = 1.0  # Height (in meters)
Nx = 20   # Number of mesh points along the x axis
Ny = 20   # Number of mesh points along the y axis
dx = lx / (Nx - 1)  # Mesh spacing in x dir
dy = ly / (Ny - 1)  # Mesh spacing in y dir
d = dx              # Since mesh spacing is the same in both dirs we use only one parameter delta (d)

# Boundary conditions
T = np.zeros(shape=(Nx, Ny))
T[:, 0] = 0     # Temperature on left
T[:, -1] = 0    # Temperature on right
T[0, :] = 100   # Temperature on top
T[-1, :] = 0    # Temperature on bottom

# Grid generation
x = np.linspace(0, lx, Nx)
y = np.linspace(0, ly, Ny)
X, Y = np.meshgrid(x, np.flipud(y))

#=======================================================================================
#Define functions and solvers

########################################### Function to perform one iteration of ADI method
def adi_iteration(T, tol=1e-5):
    Nx, Ny = T.shape

    
    for j in range(1, Ny - 1): # Step 1: Solve in x direction
        a = np.ones(Nx - 2)
        a[0] = 0
        b = -4 * np.ones(Nx - 2)
        c = np.ones(Nx - 2)
        c[-1] = 0
        d = np.zeros(Nx - 2)

        for i in range(1, Nx - 1):
            d[i - 1] = -T[i, j - 1] - T[i, j + 1]

        T[1:-1, j] = TDMA(a, b, c, d)

   
    for i in range(1, Nx - 1): # Step 2: Solve in y direction
        a = np.ones(Ny - 2)
        a[0] = 0
        b = -4 * np.ones(Ny - 2)
        c = np.ones(Ny - 2)
        c[-1] = 0
        d = np.zeros(Ny - 2)

        for j in range(1, Ny - 1):
            d[j - 1] = -T[i - 1, j] - T[i + 1, j]

        T[i, 1:-1] = TDMA(a, b, c, d)

    return T


########################################## Function to perform ADI method until convergence
def adi_solver(T, tol=1e-5):
    error = tol + 1
    count = 0

    while error > tol:
        T_g = np.copy(T)
        T = adi_iteration(T)
        error = np.max(np.abs(T - T_g))
        count += 1

    return T, count, error


########################################## TDMA function
def TDMA(a,b,c,d):
    n = len(a)
    ϕ = np.empty(n)
    P = np.empty(n)
    Q = np.empty(n)
    
    P[0] = -c[0]/b[0]
    Q[0] = d[0]/b[0]
    
    
    for i in range (1, n):
        P[i] = - c[i] / (b[i] + a[i]*P[i-1])
        Q[i] = (d[i] - a[i]*Q[i-1])/(b[i] + a[i]*P[i-1])
    
    
    ϕ[-1] = Q[-1]
    
    
    for i in range(n-2, -1, -1):
        ϕ[i] = P[i]*ϕ[i+1] + Q[i]

    return ϕ

#=======================================================================================
#Setup solver and do the calculations

T_solution, count, error = adi_solver(T)

#=======================================================================================
#Show results

#Simulation stats
print('Number of iterations:', count)
print('Value of final:', error)

# Adjust levels to every 10 degrees
levels = np.linspace(0, 100, 11)

# Plotting
plt.figure(3, dpi=180)
cp1 = plt.contourf(X, Y, T, 10, cmap='jet')
plt.colorbar()
cp1 = plt.contour(X, Y, T, 10, colors='k')
plt.clabel(cp1, inline=True, fontsize=10)
plt.show()