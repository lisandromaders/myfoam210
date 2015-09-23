Info<< "### Write fields from the FGM table ###" << endl;

volScalarField Tfgm
(
    IOobject
    (
        "T_fgm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fgmTable.getField2D(zt,pv,fgmThermoTransportIndices[1])
);
Tfgm.write();

volScalarField sourcePV
(
    IOobject
    (
        "PVsource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fgmTable.sourcePV2D(zt,pv)
);
sourcePV.write();
