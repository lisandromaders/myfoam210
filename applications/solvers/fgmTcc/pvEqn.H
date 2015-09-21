//volScalarField PVsource = pow((pv*(1-pv)),2)*fanzy;

solve
(
      fvm::ddt(rho, pv)
    + fvm::div(phi, pv)
    - fvm::laplacian(alfa, pv)
    ==
      fgmTable.sourcePV2D(zt,pv)
    
);

pv.max(0.0);

            
