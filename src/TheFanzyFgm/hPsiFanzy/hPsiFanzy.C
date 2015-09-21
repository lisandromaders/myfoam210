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
    const scalarField& hCells = h_.internalField();

    scalarField& TCells = this->T_.internalField();
    scalarField& psiCells = this->psi_.internalField();
    scalarField& muCells = this->mu_.internalField();
    scalarField& alphaCells = this->alpha_.internalField();

    const scalarField& foamCV1Cells = foamCV1_.internalField();
    const scalarField& foamCV2Cells = foamCV2_.internalField();
    
    forAll(TCells, celli)
    {
        scalar Cptab = fgmTable_.getValue2D
                      (
                          foamCV1Cells[celli],
                          foamCV2Cells[celli],
                          fgmThermoTransportIndices_[0]
                      );
        scalar Ttab = fgmTable_.getValue2D
                      (
                          foamCV1Cells[celli],
                          foamCV2Cells[celli],
                          fgmThermoTransportIndices_[4]
                      );
        scalar htab = fgmTable_.getValue2D
                      (
                          foamCV1Cells[celli],
                          foamCV2Cells[celli],
                          fgmThermoTransportIndices_[5]
                      );
        // T = T_tab - (h_tab - h)/Cp_tab;
        TCells[celli] = Ttab - (htab - hCells[celli])/Cptab;
        
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

        fvPatchScalarField& ph = h_.boundaryField()[patchi];

        fvPatchScalarField& pmu = this->mu_.boundaryField()[patchi];
        fvPatchScalarField& palpha = this->alpha_.boundaryField()[patchi];

        const fvPatchScalarField& pfoamCV1 = foamCV1_.boundaryField()[patchi];
        const fvPatchScalarField& pfoamCV2 = foamCV2_.boundaryField()[patchi];
        
        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                scalar Cptab = fgmTable_.getValue2D
                                (
                                    pfoamCV1[facei],
                                    pfoamCV2[facei],
                                    fgmThermoTransportIndices_[0]
                                );
                scalar Ttab = fgmTable_.getValue2D
                                (
                                    pfoamCV1[facei],
                                    pfoamCV2[facei],
                                    fgmThermoTransportIndices_[4]
                                );
                scalar htab = fgmTable_.getValue2D
                                (
                                    pfoamCV1[facei],
                                    pfoamCV2[facei],
                                    fgmThermoTransportIndices_[5]
                                );
                // h = h_tab + (T - T_tab)*Cp_tab;
                ph[facei] = htab + (pT[facei] - Ttab)*Cptab;
//                 ph[facei] = htab;

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
                scalar Cptab = fgmTable_.getValue2D
                                (
                                    pfoamCV1[facei],
                                    pfoamCV2[facei],
                                    fgmThermoTransportIndices_[0]
                                );
                scalar Ttab = fgmTable_.getValue2D
                                (
                                    pfoamCV1[facei],
                                    pfoamCV2[facei],
                                    fgmThermoTransportIndices_[4]
                                );
                scalar htab = fgmTable_.getValue2D
                                (
                                    pfoamCV1[facei],
                                    pfoamCV2[facei],
                                    fgmThermoTransportIndices_[5]
                                );
                pT[facei] = Ttab - (htab - ph[facei])/Cptab;

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

    h_
    (
        IOobject
        (
            "h",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->hBoundaryTypes()
    ),

    fgmTable_(fanzyLookUp),
    foamCV1_(foamCV1),
    foamCV2_(foamCV2),
    
    fgmThermoTransportIndices_(lookup("fgmThermoTransportIndices"))
{
    scalarField& hCells = h_.internalField();
    scalarField& TCells = this->T_.internalField();

    const scalarField& foamCV1Cells = foamCV1_.internalField();
    const scalarField& foamCV2Cells = foamCV2_.internalField();

    forAll(hCells, celli)
    {
        // Initialize enthalpy and temperature with values from the FGM table
        TCells[celli] = fgmTable_.getValue2D
                        (
                            foamCV1Cells[celli],
                            foamCV2Cells[celli],
                            fgmThermoTransportIndices_[4]
                        );
        hCells[celli] = fgmTable_.getValue2D
                        (
                            foamCV1Cells[celli],
                            foamCV2Cells[celli],
                            fgmThermoTransportIndices_[5]
                        );
    }

    forAll(h_.boundaryField(), patchi)
    {
        fvPatchScalarField& ph = h_.boundaryField()[patchi];
        fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        
        const fvPatchScalarField& pfoamCV1 = foamCV1_.boundaryField()[patchi];
        const fvPatchScalarField& pfoamCV2 = foamCV2_.boundaryField()[patchi];
            
        // Initialize enthalpy and temperature with values from the FGM table
        if (pT.fixesValue())
        {
            forAll(ph, facei)
            {
                scalar Cptab = fgmTable_.getValue2D
                                (
                                    pfoamCV1[facei],
                                    pfoamCV2[facei],
                                    fgmThermoTransportIndices_[0]
                                );
                scalar Ttab = fgmTable_.getValue2D
                                (
                                    pfoamCV1[facei],
                                    pfoamCV2[facei],
                                    fgmThermoTransportIndices_[4]
                                );
                scalar htab = fgmTable_.getValue2D
                                (
                                    pfoamCV1[facei],
                                    pfoamCV2[facei],
                                    fgmThermoTransportIndices_[5]
                                );
                // h = h_tab + (T - T_tab)*Cp_tab;
                ph[facei] = htab + (pT[facei] - Ttab)*Cptab;
//                 ph[facei] = htab;
            }
        }
        else
        {
            forAll(ph, facei)
            {
                pT[facei] = fgmTable_.getValue2D
                                (
                                    pfoamCV1[facei],
                                    pfoamCV2[facei],
                                    fgmThermoTransportIndices_[4]
                                );
                ph[facei] = fgmTable_.getValue2D
                                (
                                    pfoamCV1[facei],
                                    pfoamCV2[facei],
                                    fgmThermoTransportIndices_[5]
                                );
            }
        }
    }

    // Correct enthalpy and temperature boundary conditions after initilisation
    T_.correctBoundaryConditions();
    hBoundaryCorrection(h_);
    
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
         << "dataFieldName[h] = "
         << fgmTable_.dataNames()[fgmThermoTransportIndices_[5]] << nl
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


Foam::tmp<Foam::scalarField> Foam::hPsiFanzy::h
(
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th();

    const scalarField& foamCV1Cells = foamCV1_.internalField();
    const scalarField& foamCV2Cells = foamCV2_.internalField();

    forAll(T, celli)
    {
        scalar Cptab = fgmTable_.getValue2D
                      (
                          foamCV1Cells[cells[celli]],
                          foamCV2Cells[cells[celli]],
                          fgmThermoTransportIndices_[0]
                      );
        scalar Ttab = fgmTable_.getValue2D
                      (
                          foamCV1Cells[cells[celli]],
                          foamCV2Cells[cells[celli]],
                          fgmThermoTransportIndices_[4]
                      );
        scalar htab = fgmTable_.getValue2D
                      (
                          foamCV1Cells[cells[celli]],
                          foamCV2Cells[cells[celli]],
                          fgmThermoTransportIndices_[5]
                      );
        // h = h_tab + (T - T_tab)*Cp_tab;
        h[celli] = htab + (T[celli] - Ttab)*Cptab;
    }

    return th;
}


Foam::tmp<Foam::scalarField> Foam::hPsiFanzy::h
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th();

    const fvPatchScalarField& pfoamCV1 = foamCV1_.boundaryField()[patchi];
    const fvPatchScalarField& pfoamCV2 = foamCV2_.boundaryField()[patchi];
        
    forAll(T, facei)
    {
        scalar Cptab = fgmTable_.getValue2D
                        (
                            pfoamCV1[facei],
                            pfoamCV2[facei],
                            fgmThermoTransportIndices_[0]
                        );
        scalar Ttab = fgmTable_.getValue2D
                        (
                            pfoamCV1[facei],
                            pfoamCV2[facei],
                            fgmThermoTransportIndices_[4]
                        );
        scalar htab = fgmTable_.getValue2D
                        (
                            pfoamCV1[facei],
                            pfoamCV2[facei],
                            fgmThermoTransportIndices_[5]
                        );
        // h = h_tab + (T - T_tab)*Cp_tab;
        h[facei] = htab + (T[facei] - Ttab)*Cptab;
    }

    return th;
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
