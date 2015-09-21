       Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
       	IOobject
       	(
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
       	)
    );


/*    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();
*/
    volScalarField& p = thermo.p();
    volScalarField& h = thermo.h();
    const volScalarField& psi = thermo.psi();

    dimensionedScalar fanzy
    (
       	transportProperties.lookup("fanzy")
    );

//-------------- IMPLEMENTACAO FGM-TCC -------------

// Vou comecar com o termo difusivo constante.. Depois faco a implementacao do termo correto olhando as variaveis na tabela do FGM...
// Neste caso, a minha equacao possui um termo de difusao que eu chamei de DT, e ele procura no arquivo transportProperties pelo nome DT
// e associa o valor encontrado a esta variavel do tipo dimensionedScalar... Eh um termo constante e as unidades sao as unidades basicas
// no SI de qualquer termo de difusao.

    dimensionedScalar DT
    (
        transportProperties.lookup("DT")
    );

// Definindo o objeto da fracao de mistura zt como sendo um volScalarField

    Info<< "Reading field zt\n" << endl;
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

    volScalarField rho_output
    (
        IOobject
        (
            "rho_output",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );


    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    volScalarField alfa
    (
     	IOobject
        (
            "alfa",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
	thermo.alpha()
    );


    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    #include "compressibleCreatePhi.H"

    dimensionedScalar rhoMax(pimple.dict().lookup("rhoMax"));
    dimensionedScalar rhoMin(pimple.dict().lookup("rhoMin"));

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    Info<< "Creating field dpdt\n" << endl;
    volScalarField dpdt("dpdt", fvc::ddt(p));

    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K("K", 0.5*magSqr(U));
