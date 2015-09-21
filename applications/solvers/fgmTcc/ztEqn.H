//Resolve a equacao de transporte da fracao de mistura

//ver depois da utilizacao dos smartpointers (ver quando se aplica)
// e se eh necessario para este caso.. No caso, seria assim:
// tmp<fvScalarMatrix> ztEqn ao inves do que esta abaixo

fvScalarMatrix ztEqn
(
    fvm::ddt(rho, zt)
    + fvm::div(phi, zt)
    - fvm::laplacian(DT, zt)
);

ztEqn.solve();
