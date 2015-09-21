Info<< "Initialize tccLookUp\n" << endl;
tccLookUp fgmTable(mesh);

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

Info<< "Reading FGM control variable 1\n" << endl;
volScalarField cv1
(
    IOobject
    (
        "CV1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE // switch to AUTO_WRITE if transport equation is solved
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
tccFanzy thermo
                (
                    mesh,
                    fgmTable,
                    cv1,
                    pv
                );
const volScalarField& T = thermo.T(); // a funcao T() eh uma funcao membro da classe basicThermo. Ela eh base para basicPsiThermo, que eh base para a classe hPsiFanzy. Essa funcao retorna o membro provado T_, que eh construido com o objeto IOobject, que deve ser lido e AUTO_WRITE. Porem,ainda nao entendi pq ele le a temperatura se em teoria a equacao da energia nao eh resolvida. Logo abaixo ele associara o objeto cv1 a h_, que eh o resultado da funcao membro h() (esta eh funcao membro da classe hPsiFanzy). O membro privado h_ eh construido no arquivo hPsiFanzy.C, apos os : da inicializacao do construtor. Tambem fica a duvida de porque a variavel de controle eh igualada a entalpia neste caso...

cv1 = thermo.h(); // do this ONLY if cv1 equals enthalpy! 

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
