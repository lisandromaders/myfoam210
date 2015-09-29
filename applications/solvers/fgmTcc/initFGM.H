Info<< "Initialize fanzyLookUp\n" << endl;
fanzyLookUp fgmTable(mesh);

// List with FGM thermo and transport quantities indices (only for DEBUG output!)
labelList fgmThermoTransportIndices;
{
    IOdictionary fgmThermoProperties
    (
        IOobject
        (
            "thermophysicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    fgmThermoTransportIndices = labelList
                                (
                                    fgmThermoProperties.lookup("fgmThermoTransportIndices")
                                );
}

Info<< "Reading FGM control variable zt\n" << endl;
volScalarField zt
(
    IOobject
    (
        "zt",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading FGM progress variable\n" << endl;
volScalarField pv
(
    IOobject
    (
        "PV",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// COSTRUCAO DO OBJETO thermo QUE SERA USADO PARA CHAMAR TODAS AS FUNCOES DE CALCULO DAS PROPRIEDADES QUE PERTENCEM A CLASSE hPsiFanzy


Info<< "Initialize fgmThermo\n" << endl;
hPsiFanzy thermo
                (
                    mesh,
                    fgmTable,
                    zt,
                    pv
                );

const volScalarField& T = thermo.T(); 
const volScalarField& rho = thermo.rho(); 
//const volScalarField& kappa = thermo.kappa(); 
//const volScalarField& Cp = fgmTable.getField2D(zt,pv,fgmThermoTransportIndices[2]);
const volScalarField& Dt = thermo.Dt(); 

{
    Info<< "### Write initial temperature field ###" << endl;
    volScalarField Tinit
    (
        IOobject
        (
            "Tinit",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T
    );
    Tinit.write();

    Info<< "### Write initial density field ###" << endl;
    volScalarField rhoinit
    (
        IOobject
        (
            "rhoinit",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho
    );
    rhoinit.write();

//   Info<< "### Write initial kappa field ###" << endl;
//   volScalarField kappainit
//   (
//       IOobject
//       (
//           "kappainit",
//           runTime.timeName(),
//           mesh,
//           IOobject::NO_READ,
//           IOobject::NO_WRITE
//       ),
//       kappa
//   );
//   kappainit.write();

//   Info<< "### Write initial Cp field ###" << endl;
//   volScalarField Cpinit
//   (
//       IOobject
//       (
//           "Cpinit",
//           runTime.timeName(),
//           mesh,
//           IOobject::NO_READ,
//           IOobject::NO_WRITE
//       ),
//       Cp
//   );
//   Cpinit.write();

    Info<< "### Write initial Dt field ###" << endl;
    volScalarField Dtinit
    (
        IOobject
        (
            "Dtinit",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Dt
    );
    Dtinit.write();

  Info<< "### Initial T min / max = " << min(T).value() << " / "
        << max(T).value() << endl;
}
