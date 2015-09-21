{
    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, h)
      + fvm::div(phi, h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
        dpdt
      - (fvc::ddt(rho, K) + fvc::div(phi, K))
    );

    hEqn.relax();
    hEqn.solve();

    cv1 = h; // do this ONLY if cv1 equals enthalpy! 
    thermo.correct();

    Info<< "T min / max   = " << min(T).value() << " / "
        << max(T).value() << endl;
}
