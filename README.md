# ♨ 2D Heat Transfer in a Surface Domain (Gauss-Seidel and ADI method) ♨
🟢 Python script to solve the 2D heat equation and gain temperature distribution contours, using Gauss-Seidel and ADI (Alternating-direction implicit) method .  
🟢 This solution is based on finite difference method.  

  
## 🧬 Input variables to define the problem:  
`lx`: Width (in meters)  
`ly`: Height (in meters)  
`Nx`: Number of nodes along the x axis  
`Ny`: Number of nodes along the y axis  

`T[:, 0]`: Temperature on left  
`T[:, -1]`: Temperature on right  
`T[0, :]`: Temperature on top  
`T[-1, :]`: Temperature on bottom  
    
## 🤖 Usage  
Just replace the existing values of input parameters with your desired values and run the code.
