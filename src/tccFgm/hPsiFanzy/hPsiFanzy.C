/*---------------------------------------------------------------------------*\
  
 ####################                                   ###################
 #                  #                                   #                 #   
 #       /||=====   #                                   #   /||  |||||||  #
 #      //||        #                                   #  //||       ||  #
 #     // ||        #     The FanzyFgm tool             #    ||       ||  #
 #    //==||===     #                                   #    ||    |||||  #
 #   //   ||        #     Copyright (C) 2013 by a.f.    #    ||       ||  #
 #  //    ||anzy    #                                   #    ||       ||  #
 #                  #                                   #    ||   ||||||  #
 ####################                                   ===================

-------------------------------------------------------------------------------
License
    This file is part of the TheFanzyFgm library.

       
\*---------------------------------------------------------------------------*/

#include "hPsiFanzy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hPsiFanzy, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::hPsiFanzy::calculate()
{
    scalarField& TCells = this->T_.internalField();
    scalarField& psiCells = this->psi_.internalField();
    scalarField& muCells = this->mu_.internalField();
    scalarField& alphaCells = this->alpha_.internalField();

    const scalarField& foamCV1Cells = foamCV1_.internalField();
    const scalarField& foamCV2Cells = foamCV2_.internalField();
    
    forAll(TCells, celli)
    {
        scalar Ttab = fgmTable_.getValue2D
                      (
                          foamCV1Cells[celli],
                          foamCV2Cells[celli],
                          fgmThermoTransportIndices_[4]
                      );
        TCells[celli] = Ttab;
        
        // psi = 1.0/(R*T);
       psiCells[celli] = 1.0/(TCells[celli]*fgmTable_.getValue2D
                                            (
                                                foamCV1Cells[celli],
                                                foamCV2Cells[celli],
                                                fgmThermoTransportIndices_[1]
                                            ));

        muCells[celli] = fgmTable_.getValue2D
                         (
                             foamCV1Cells[celli],
                             foamCV2Cells[celli],
                             fgmThermoTransportIndices_[2]
                         );
        alphaCells[celli] = fgmTable_.getValue2D
                            (
                                foamCV1Cells[celli],
                                foamCV2Cells[celli],
                                fgmThermoTransportIndices_[3]
                            );
    }

    forAll(T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& ppsi = this->psi_.boundaryField()[patchi];

        fvPatchScalarField& pmu = this->mu_.boundaryField()[patchi];
        fvPatchScalarField& palpha = this->alpha_.boundaryField()[patchi];

        const fvPatchScalarField& pfoamCV1 = foamCV1_.boundaryField()[patchi];
        const fvPatchScalarField& pfoamCV2 = foamCV2_.boundaryField()[patchi];
        
        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                scalar Ttab = fgmTable_.getValue2D
                                (
                                    pfoamCV1[facei],
                                    pfoamCV2[facei],
                                    fgmThermoTransportIndices_[4]
                                );

               ppsi[facei] = 1.0/(pT[facei]*fgmTable_.getValue2D
                                               (
                                                   pfoamCV1[facei],
                                                   pfoamCV2[facei],
                                                   fgmThermoTransportIndices_[1]
                                               ));

                pmu[facei] = fgmTable_.getValue2D
                                (
                                    pfoamCV1[facei],
                                    pfoamCV2[facei],
                                    fgmThermoTransportIndices_[2]
                                );
                palpha[facei] = fgmTable_.getValue2D
                                (
                                    pfoamCV1[facei],
                                    pfoamCV2[facei],
                                    fgmThermoTransportIndices_[3]
                                );
            }
        }
        else
        {
            forAll(pT, facei)
            {
                scalar Ttab = fgmTable_.getValue2D
                                (
                                    pfoamCV1[facei],
                                    pfoamCV2[facei],
                                    fgmThermoTransportIndices_[4]
                                );

                pT[facei] = Ttab;

               ppsi[facei] = 1.0/(pT[facei]*fgmTable_.getValue2D
                                               (
                                                   pfoamCV1[facei],
                                                   pfoamCV2[facei],
                                                   fgmThermoTransportIndices_[1]
                                               ));

                pmu[facei] = fgmTable_.getValue2D
                                (
                                    pfoamCV1[facei],
                                    pfoamCV2[facei],
                                    fgmThermoTransportIndices_[2]
                                );
                palpha[facei] = fgmTable_.getValue2D
                                (
                                    pfoamCV1[facei],
                                    pfoamCV2[facei],
                                    fgmThermoTransportIndices_[3]
                                );
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hPsiFanzy::hPsiFanzy
(
    const fvMesh& mesh,
    const fanzyLookUp& fanzyLookUp,
    const volScalarField& foamCV1,
    const volScalarField& foamCV2
)
:
    basicPsiThermo(mesh),

    fgmTable_(fanzyLookUp),
    foamCV1_(foamCV1),
    foamCV2_(foamCV2),
    
    fgmThermoTransportIndices_(lookup("fgmThermoTransportIndices"))
{
    scalarField& TCells = this->T_.internalField();

    const scalarField& foamCV1Cells = foamCV1_.internalField();
    const scalarField& foamCV2Cells = foamCV2_.internalField();

    forAll(TCells, celli)
    {
        // Initialize enthalpy and temperature with values from the FGM table
        TCells[celli] = fgmTable_.getValue2D
                        (
                            foamCV1Cells[celli],
                            foamCV2Cells[celli],
                            fgmThermoTransportIndices_[4]
                        );
    }

    forAll(T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        
        const fvPatchScalarField& pfoamCV1 = foamCV1_.boundaryField()[patchi];
        const fvPatchScalarField& pfoamCV2 = foamCV2_.boundaryField()[patchi];
            
        // Initialize enthalpy and temperature with values from the FGM table
        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                scalar Ttab = fgmTable_.getValue2D
                                (
                                    pfoamCV1[facei],
                                    pfoamCV2[facei],
                                    fgmThermoTransportIndices_[4]
                                );
            }
        }
        else
        {
            forAll(pT, facei)
            {
                pT[facei] = fgmTable_.getValue2D
                                (
                                    pfoamCV1[facei],
                                    pfoamCV2[facei],
                                    fgmThermoTransportIndices_[4]
                                );
            }
        }
    }

    // Correct enthalpy and temperature boundary conditions after initilisation
    T_.correctBoundaryConditions();
    
    Info << "dataFieldName[Cp] = "
         << fgmTable_.dataNames()[fgmThermoTransportIndices_[0]] << nl
         << "dataFieldName[R] = "
         << fgmTable_.dataNames()[fgmThermoTransportIndices_[1]] << nl
         << "dataFieldName[mu] = "
         << fgmTable_.dataNames()[fgmThermoTransportIndices_[2]] << nl
         << "dataFieldName[alpha] = "
         << fgmTable_.dataNames()[fgmThermoTransportIndices_[3]] << nl
         << "dataFieldName[T] = "
         << fgmTable_.dataNames()[fgmThermoTransportIndices_[4]] << nl
         << endl;
    
    calculate();

    // Switch on saving old time
    this->psi_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hPsiFanzy::~hPsiFanzy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hPsiFanzy::correct()
{
    if (debug)
    {
        Info<< "entering hPsiFanzy::correct()" << endl;
    }

    // force the saving of the old-time values
    this->psi_.oldTime();

    calculate();

    if (debug)
    {
        Info<< "exiting hPsiFanzy::correct()" << endl;
    }
}


Foam::tmp<Foam::scalarField> Foam::hPsiFanzy::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCp(new scalarField(T.size()));
    scalarField& cp = tCp();

    const fvPatchScalarField& pfoamCV1 = foamCV1_.boundaryField()[patchi];
    const fvPatchScalarField& pfoamCV2 = foamCV2_.boundaryField()[patchi];

    forAll(T, facei)
    {
        cp[facei] = fgmTable_.getValue2D
                    (
                        pfoamCV1[facei],
                        pfoamCV2[facei],
                        fgmThermoTransportIndices_[0]
                    );
    }

    return tCp;
}


Foam::tmp<Foam::volScalarField> Foam::hPsiFanzy::Cp() const
{
    const fvMesh& mesh = this->T_.mesh();
    const scalarField& foamCV1Cells = foamCV1_.internalField();
    const scalarField& foamCV2Cells = foamCV2_.internalField();

    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimEnergy/dimMass/dimTemperature
        )
    );

    volScalarField& cp = tCp();

    forAll(this->T_, celli)
    {
        cp[celli] = fgmTable_.getValue2D
                    (
                        foamCV1Cells[celli],
                        foamCV2Cells[celli],
                        fgmThermoTransportIndices_[0]
                    );
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pCp = cp.boundaryField()[patchi];

        const fvPatchScalarField& pfoamCV1 = foamCV1_.boundaryField()[patchi];
        const fvPatchScalarField& pfoamCV2 = foamCV2_.boundaryField()[patchi];

        forAll(pT, facei)
        {
            pCp[facei] = fgmTable_.getValue2D
                         (
                             pfoamCV1[facei],
                             pfoamCV2[facei],
                             fgmThermoTransportIndices_[0]
                         );
        }
    }

    return tCp;
}


// ************************************************************************* //
