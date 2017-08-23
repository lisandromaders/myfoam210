/*      
Class
    Foam::fanzyLookUp

Description
    Class for the FGM table looking-up

SourceFiles
    fanzyLookUpI.H
    fanzyLookUp.C
    fanzyLookUpAuxiliaryFunctions.C
*/

#ifndef fanzyLookUp_H
#define fanzyLookUp_H

#include "IOdictionary.H"
#include "Switch.H"
#include "scalar.H"
#include "scalarField.H"
#include "volFields.H"
#include "PtrList.H"
#include "zeroGradientFvPatchFields.H"
#include "OFstream.H"
#include "IFstream.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fanzyLookUp Declaration
\*---------------------------------------------------------------------------*/

class fanzyLookUp
:
    public IOdictionary
{
    // Private Member Functions

        //- Disallow copy constructor
        fanzyLookUp(const fanzyLookUp&);

        //- Disallow default bitwise assignment
        void operator=(const fanzyLookUp&);


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Switch debug output on/off
        Switch debugOutput_;
        
        //- Filename of the FGM table
        const fileName fgmFile_;

        //- Index of sourcePV
        const label sourcePVindex_;
        
      // Variables read from the FGM table
        
          //- Dimension of Control Variable (CV) 1
          label nFgmCV1_;

          //- Dimension of Control Variable (CV) 2
          label nFgmCV2_;

          //- Dimension of the data set (i.e. number of tabulated variables)
          label nFgmData_;

          //- Names of variables in the FGM table
          List<word> fgmVariableNames_;

          //- Control Variable (CV) 1
          List<scalar> fgmCV1_;
          
          //- Control Variable (CV) 2
          List<scalar> fgmCV2_;
          
          //- Max. value of Control Variable (CV) 1
          scalar maxCV1_;

          //- Max. value of Control Variable (CV) 2
          scalar maxCV2_;

          //- FGM data table
          List<List<scalar> >  fgmData_;
          
    // Protected Member Functions

        //- Read and discard to newline
        static void readToNewline(IFstream&);

        //- Read FGM table from file
        void readFgmFile(const fileName);

public:

    //- Runtime type information
    TypeName("fanzyLookUp");


    // Constructors

        //- Construct from components
        fanzyLookUp(const fvMesh&);


    //- Destructor
    ~fanzyLookUp();


    // Access general FGM table data

        //- Return const access to the mesh database
        inline const fvMesh& mesh() const;

        //- Switch debug output on/off
        inline Switch debugOutput() const;

        //- Names of the data variables in the FGM table
        inline List<word> fgmVariableNames() const;

        //- Names of the Control Variables
        inline List<word> cvNames() const;
        
        //- Names of the data variables in the FGM table
        inline List<word> dataNames() const;

    // Access FGM variables and fields
        
        //- Return source term for PV from the FGM table
        tmp<volScalarField> sourcePV2D
        (
            const volScalarField& foamCV1,
            const volScalarField& foamCV2
        ) const;
        
        //- Get a field with values from the FGM table
        tmp<volScalarField> getField2D
        (
            const volScalarField&,
            const volScalarField&,
            const label
        ) const;

        //- Get a value from the FGM table (with switch for PV scaling)
        scalar getValue2D
        (
            const scalar,
            const scalar,
            const label
        ) const;

        //- Interpolate FGM variable for CV1 & CV2
        scalar interpolate2D
        (
            const scalar,
            const scalar,
            const label
        ) const;

    // Auxiliary Functions
        
        //- Write FGM table to file
        void writeFgmFile(const fileName);
        
        //- Write FGM table to files in OpenFOAM format
        void writeFgmFileOpenfoam();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fanzyLookUpI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
