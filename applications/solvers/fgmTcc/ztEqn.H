fvScalarMatrix ztEqn
(
    fvm::ddt(rho, zt)
    + fvm::div(phi, zt)
    - fvm::laplacian(Dt, zt)
);

ztEqn.solve();

