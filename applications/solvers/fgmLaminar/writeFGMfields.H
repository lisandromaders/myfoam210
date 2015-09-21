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
    fgmTable.getField2D(h,pv,fgmThermoTransportIndices[4])
);
Tfgm.write();

volScalarField hfgm
(
    IOobject
    (
        "h_fgm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fgmTable.getField2D(h,pv,fgmThermoTransportIndices[5])
);
hfgm.write();

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
    fgmTable.sourcePV2D(h,pv)
);
sourcePV.write();
