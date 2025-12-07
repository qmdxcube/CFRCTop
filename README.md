CFRCTop: An Efficient MATLAB Implementation for Topology Optimization of Continuous Fiber-Reinforced Composite Structures

We present CFRCTop, a MATLAB implementation for topology optimization of continuous fiber-reinforced composite structures. The implementation includes density and fiber-orientation filtering, finite element analysis, sensitivity analysis, design variable updating, verification of optimality of fiber orientations, and visualization of results. This code is built upon the well-known topology optimization code top88. The template stiffness matrices (TSMs)-based method is employed for efficient finite element analysis and sensitivity analysis. The density and fiber-orientation variables are updated separately. Visualization of spatially varying fiber orientations is provided. Extensions to solving various problems are also discussed. Computational performance and scalability are studied to showcase the high efficiency of this implementation. CFRCTop is intended for students and newcomers in the field of topology optimization.

The paper describing the theoretical foundation, implementation, numerical examples and extensions can be found at https://www.mdpi.com/2076-3417/15/17/9242.

Zhao J. CFRCTop: An Efficient MATLAB Implementation for Topology Optimization of Continuous Fiber-Reinforced Composite Structures. Applied Sciences, 2025, 15(17): 9242.

1. Optimized design for the MBB problem. (a) Topology and fiber orientations. (b) Collinearity of fiber orientations with the dominant principal stress directions.
<img width="2764" height="2545" alt="applsci-15-09242-g006" src="https://github.com/user-attachments/assets/33c4b0bc-4187-4f05-9faf-4f9e207bde3d" />

2. Optimized design for the MBB problem obtained by sequentially optimizing of the topology and fiber orientations. (a) Topology and fiber orientations. (b) Collinearity of fiber orientations with the dominant principal stress directions.
<img width="2659" height="2454" alt="applsci-15-09242-g008" src="https://github.com/user-attachments/assets/89e755e1-24f3-422a-84f9-6a23ecad4f8c" />

3. Optimized design obtained from different initialization of the fiber orientations. (a,c,e) The optimized designs for X, Y and random initializations, respectively. (b,d,f) Collinearity of optimized fiber orientations with the dominant principal stress directions for X, Y and random initializations, respectively.
<img width="3987" height="3243" alt="applsci-15-09242-g010" src="https://github.com/user-attachments/assets/f0d38953-76f1-4967-a8f6-0f39646c353f" />

4. Optimized design for the short cantilever problem with a passive region. (a) Topology and fiber orientations. (b) Collinearity of fiber orientations with the dominant principal stress directions.
<img width="3502" height="1576" alt="applsci-15-09242-g020" src="https://github.com/user-attachments/assets/6adc3f3b-2a8e-463a-b0eb-a4739744b97e" />

5. Optimized design for the cantilever problem with two load cases. For both load cases, the load directions are modified to be horizontal. (a) Topology and fiber orientations. (b) Collinearity of fiber orientations with the dominant principal stress directions under the first load case. (c) Collinearity of fiber orientations with the dominant principal stress directions under the second load case.
<img width="4172" height="1454" alt="applsci-15-09242-g018" src="https://github.com/user-attachments/assets/bd68156a-4f02-4151-8606-40266a9c8f28" />

6. Average computation time per iteration for the MBB problem with different resolutions and filter radii.
   ![applsci-15-09242-g021-550](https://github.com/user-attachments/assets/51529db2-772b-4254-ab97-e00d5cee98ad)

8. Average computation time per iteration obtained with the updated version of the codes for the MBB problem with different resolutions and filter radii.
<img width="2411" height="2009" alt="applsci-15-09242-g022" src="https://github.com/user-attachments/assets/9ddc0eb8-750c-419c-b3dd-4c5a31ae0ed4" />
