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

volScalarField Dtfgm
(
    IOobject
    (
        "Dt_fgm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fgmTable.getField2D(zt,pv,fgmThermoTransportIndices[3])/(fgmTable.getField2D(zt,pv,fgmThermoTransportIndices[2])+0.01)
);
Dtfgm.write();

volScalarField rhofgm
(
    IOobject
    (
        "rho_fgm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fgmTable.getField2D(zt,pv,fgmThermoTransportIndices[0])
);
rhofgm.write();

volScalarField cpfgm
(
    IOobject
    (
        "cp_fgm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fgmTable.getField2D(zt,pv,fgmThermoTransportIndices[2])
);
cpfgm.write();

volScalarField kappafgm
(
    IOobject
    (
        "kappa_fgm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fgmTable.getField2D(zt,pv,fgmThermoTransportIndices[3])
);
kappafgm.write();

volScalarField mufgm
(
    IOobject
    (
        "mu_fgm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fgmTable.getField2D(zt,pv,fgmThermoTransportIndices[4])
);
mufgm.write();
