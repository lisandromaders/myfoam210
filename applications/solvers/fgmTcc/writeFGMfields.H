Info<< "### Write fields from the FGM table ###" << endl;

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

