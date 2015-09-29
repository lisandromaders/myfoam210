//volScalarField PVsource = pow((pv*(1-pv)),2)*fanzy;

//solve
//(
//      fvm::ddt(rho, pv)
//    + fvm::div(phi, pv)
//    - fvm::laplacian(Dt, pv)
//   ==
//     fgmTable.sourcePV2D(zt,pv)
//    
//);
fvScalarMatrix pvEqn
(
      fvm::ddt(rho, pv)
    + fvm::div(phi, pv)
    - fvm::laplacian(Dt, pv)
   ==
     fgmTable.sourcePV2D(zt,pv)
    
);

pvEqn.relax();
pvEqn.solve();

//pv.max(0.0);

            
