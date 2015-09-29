if (fgmTable.debugOutput())
{
    Info<< "### Write thermo property fields ###" << endl;
    Info<< "Writing Cp" << endl;
    volScalarField CpDebug
    (
        IOobject
        (
            "Cp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fgmTable.getField2D(zt,pv,fgmThermoTransportIndices[2])
    );
    CpDebug.write();

    Info<< "Writing kappa" << endl;
    volScalarField kappaDebug
    (
        IOobject
        (
            "kappa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fgmTable.getField2D(zt,pv,fgmThermoTransportIndices[3])
    );
    kappaDebug.write();

    Info<< "Writing mu" << endl;
    volScalarField muDebug
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.mu()
    );
    muDebug.write();

    Info<< "Writing rho" << endl;
    volScalarField rhoDebug
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.rho()
    );
    rhoDebug.write();

   Info<< "Writing Dt" << endl;
   volScalarField DtDebug
   (
       IOobject
       (
           "Dt",
           runTime.timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::NO_WRITE
       ),
       thermo.Dt()
   );
   DtDebug.write();
}
