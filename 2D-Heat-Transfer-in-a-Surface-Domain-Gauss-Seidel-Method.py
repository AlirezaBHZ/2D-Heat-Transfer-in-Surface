import numpy as np
import matplotlib.pyplot as plt

#=======================================================================================
# Input parameters

# Geometrical properties
lx = 1.0 # Width (in meters)
ly = 1.0 # Height (in meters)
Nx = 20  # Number of mesh points along the x axis
Ny = 20  # Number of mesh points along the y axis
dx = lx / (Nx - 1) # Mesh spacing in x dir
dy = ly / (Ny - 1) # Mesh spacing in y dir
d = dx # Since mesh spacing is the same in both dirs we use only one parameter delta (d)

# Boundary conditions
T = np.zeros(shape=(Nx, Ny))
T[:, 0] = 0 # Temperature on left
T[:, -1] = 0 # Temperature on right
T[0, :] = 100 # Temperature on top
T[-1, :] = 0 # Temperature on bottom

#=======================================================================================
# Grid generation

x = np.linspace(0, lx, Nx)
y = np.linspace(0, ly, Ny)
X, Y = np.meshgrid(x, np.flipud(y))

#=======================================================================================
#Setup solver and do the calculations

# Guess array for comparison
Tg = T.copy()

# Initial error or entry in loop
error = 1

# Iteration counter
count = 1

# Comparison loop
while error > 1.E-5:
   
    T_prev = T.copy() # Copying T to store previous values

  
    for i in range(1, Nx - 1):  # Gauss-Seidel Iteration
        for j in range(1, Ny - 1):
            T[i, j] = (T_prev[i, j - 1] + T[i, j + 1] + T_prev[i - 1, j] + T[i + 1, j]) / 4

   
    error = np.sqrt(np.sum(np.abs(T - Tg)**2)) # Evaluating error


    Tg = T.copy() # Updating guess for next iteration


    count += 1    # Incrementing counter
    
#=======================================================================================
#Show results

#Simulation stats
print('Number of iterations:', count)
print('Value of final:', error)

# Result Plotting (Temperature contour plot)
plt.figure(3, dpi=180)
cp1 = plt.contourf(X, Y, T, 10, cmap='jet')
plt.colorbar()
cp1 = plt.contour(X, Y, T, 10, colors='k')
plt.clabel(cp1, inline=True, fontsize=10)
plt.show()