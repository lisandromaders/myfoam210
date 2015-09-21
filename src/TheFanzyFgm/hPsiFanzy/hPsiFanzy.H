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

Class
    Foam::hPsiFanzy

Description
    * Enthalpy (h) for a mixture based on compressibility
    * Temperature calculation based on 'transported' h and tabulated T, h and Cp
    * Transport and thermo properties obtained from a FGM table

SourceFiles
    hPsiFanzy.C

\*---------------------------------------------------------------------------*/

#ifndef hPsiFanzy_H
#define hPsiFanzy_H

#include "basicPsiThermo.H"
#include "fanzyLookUp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class hPsiFanzy Declaration
\*---------------------------------------------------------------------------*/

class hPsiFanzy
:
    public basicPsiThermo	//inherited class (does not have in OF 2.3.0)
{
    // Private data DECLARATION --> they will be defined in the .C file. 
    // Ex: the object h_ is declared of being of type volScalarField, but it will be constructed from IOobject class in the .C file.

        //- Sensible enthalpy field [J/kg]
        volScalarField h_;

        //- Table look up routines
        const fanzyLookUp& fgmTable_;
        
        //- Reference to CV1 & CV2
        const volScalarField& foamCV1_;
        const volScalarField& foamCV2_;
        
        //- FGM table index maps (CpIndex, RIndex, muIndex, alphaIndex)
        const labelList fgmThermoTransportIndices_;
        
        
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();

        //- Construct as copy (not implemented)
        hPsiFanzy(const hPsiFanzy&);


public:

    //- Runtime type information
    TypeName("hPsiFanzy");


    // Constructors

        //- Construct from mesh
        hPsiFanzy
        (
            const fvMesh&,
            const fanzyLookUp&,
            const volScalarField&,
            const volScalarField&
        );


    //- Destructor
    ~hPsiFanzy();


    // Member functions

        //- Update properties
        void correct();


        // Access to thermodynamic state variables

            //- Sensible enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            volScalarField& h()
            {
                return h_;
            }

            //- Sensible enthalpy [J/kg]
            const volScalarField& h() const
            {
                return h_;
            }


        // Fields derived from thermodynamic state variables

            //- Enthalpy for cell-set [J/kg]
            tmp<scalarField> h
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Enthalpy for patch [J/kg]
            tmp<scalarField> h
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            tmp<volScalarField> Cp() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

// #include "hPsiFanzyI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
