An immersed boundary simulation of the Red Blood Cell (RBC) and endothelial surface layer (ESL) interaction using the coarse-grained model.


List of MATLAB programs:
1. To run the simulation and change parameters, use the main routine:

     IBdriver.m

2. Fluid solver:

    solveNonDimenNSeqn.m

3. The standard discrete 4-point delta function for interpolating and spreading Lagrangian forces

    evalPhi.m

4. Evaluate the modified discrete 4-point delta function with no-slip BCs on the top and bottom of the channel
  
    evalDeltaPhysBCs.m

5. Spread the Lagrangian force density on the Eulerian grid using the modified discrete 4-point delta function 

   spreadLagForcePhysBCs.m

6. Evaluate the Lagrangian force density associated with the RBC
  
   getLagForceRLA.m

7. Compute the porous slip velocity based on Darcy's law
  
   getPorousSlipV.m

8. Discrete differential operators for the fluid solver:
  
   D02DDirPer.m, D02DNeumann.m, L2DNeumannPer.m, Lu2DDirPer.m

9. Helper routine:
  
   blktridiag.m
