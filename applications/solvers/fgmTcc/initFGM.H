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



  Info<< "### Initial T min / max = " << min(T).value() << " / "
        << max(T).value() << endl;
}
