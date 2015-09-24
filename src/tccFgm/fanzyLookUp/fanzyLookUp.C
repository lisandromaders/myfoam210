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

#include "fanzyLookUp.H"
#include "fvMesh.H"
#include "Time.H"
#include "ListOps.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(fanzyLookUp, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fanzyLookUp::fanzyLookUp
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "fgmProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    mesh_(mesh),
    
    debugOutput_(lookup("debugOutput")),
    
    fgmFile_(fileName(lookup("fgmFile")).expand()),

    scalePV_(lookup("scalePV")),
    scalePVFile_(fileName(lookup("scalePVFile")).expand()),
    sourcePVindex_(readLabel(lookup("sourcePVindex")))
{
    Info<< "fanzyLookUp: FGM file = " << fgmFile_ << endl;
    // Read and write the FGM table on construction (class instantiated)
    readFgmFile(fgmFile_);
    fileName fgmTableOut = path()/"FGM.out";
//     fileName fgmTableOut = "FGMTable.out";
    writeFgmFile(fgmTableOut);
    
    // Read and write the ScalePV table if switch 'scalePV' is 'on'
    if (scalePV_)
    {
        readScalePVFile(scalePVFile_);
        fileName scalePVTableOut = path()/"ScalePV.out";
//         fileName scalePVTableOut = "ScalePV.out";
        writeScalePVFile(scalePVTableOut);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fanzyLookUp::~fanzyLookUp()
{}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

void Foam::fanzyLookUp::readToNewline(IFstream& is)
{
    char ch = '\n';
    do
    {
        (is).get(ch);
    }
    while ((is) && ch != '\n');
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
FLAMELET GENERATED MANIFOLD

[NUMBER_FLAMELETS]
 <nFlamelets>
[NUMBER_PV]
 <nPV>
[NUMBER_VARIABLES]
  <nVariables>

[DATA]
<lableVar01> <lableVar02> ... <lableVarN>
  <var01>   <var02>   ...   <varN>
\*---------------------------------------------------------------------------*/
void Foam::fanzyLookUp::readFgmFile(const fileName inputName)
{
    IFstream is(inputName);
    char c;
    label nFlamelets, nPV, nVariables;
    
    if (is.good())
    {
        Info<< "fanzyLookUp: Read FGM file..." << endl;
    }
    else
    {
        FatalErrorIn
        (
            "fanzyLookUp::readFgmFile(const fileName& inputName)"
        )
        << "cannot read " << is.name() << abort(FatalError);
    }

    //- Read header
    readToNewline(is);                       // FLAMELET GENERATED MANIFOLD
    readToNewline(is);                       // 
    readToNewline(is);                       // [NUMBER_FLAMELETS]
    is.read(nFlamelets); readToNewline(is);  // <nFlamelets>
    readToNewline(is);                       // [NUMBER_PV]
    is.read(nPV); readToNewline(is);         // <nPV>
    readToNewline(is);                       // [NUMBER_VARIABLES]
    is.read(nVariables); readToNewline(is);  // <nVariables>
    readToNewline(is);                       // 
    readToNewline(is);                       // [DATA]

    //- Read variable names
    fgmVariableNames_ = List<word>(nVariables);
    forAll(fgmVariableNames_, i)
    {
        is.read(fgmVariableNames_[i]);
        is.get(c);
       if (c != ' ')
       {
           FatalErrorIn
           (
               "fanzyLookUp::readFgmFile(const fileName& inputName)"
           )
           << "cannot read FGM variable names" << abort(FatalError);
       };
    }
    readToNewline(is);
    
    /* Initialize the arrays
     * =====================
     * Put Control Variable 1 and 2 to global 1D arrays (fgmCV1_ and fgmCV2_)
     * Put FGM data to a global 3D array (fgmData_): 
     * [        Dim 1  x        Dim 2  x                           Dim 3  ] = 
     * [ (Size of CV1) x (Size of CV2) x (Number of variables in the FGM) ]
     */
    nFgmCV1_ = nFlamelets;
    nFgmCV2_ = nPV;
    nFgmData_ = nVariables-2;

    scalar rawCV1;
    fgmCV1_ = List<scalar>(nFgmCV1_);
    fgmCV2_ = List<scalar>(nFgmCV2_);
    fgmData_  = List<List<List<scalar> > >
                (
                    nFgmCV1_,
                    List<List<scalar> >
                    (
                        nFgmCV2_,
                        List<scalar>
                        (
                            nFgmData_
                        )
                    )
                );

    //- Write out some information
    Info << "               Dimension of Control Variable 1: " 
         << fgmData_.size() << nl
         << "               Dimension of Control Variable 2: " 
         << fgmData_[0].size() << nl
         << "               Number of tabulated variables:   " 
         << fgmData_[0][0].size() << nl
         << endl;
            
    Info << "               fgmVariableNames[" << 0 << "] (cv1) = " 
         << fgmVariableNames_[0] << nl
         << "               fgmVariableNames[" << 1 << "] (cv2) = " 
         << fgmVariableNames_[1] << nl;
    for (label i=0; i<nFgmData_; i++)
    {
        Info << "               fgmVariableNames[" << i+2 << "] (dataField[" 
             << i << "]) = " << fgmVariableNames_[i+2] << nl;
    }
    Info << endl;
    
    //- Read the FGM data
	    forAll(fgmData_, i)
	    {
		forAll(fgmData_[0], j)
		{
		    is.read(rawCV1);
		    is.read(fgmCV2_[j]);
		    
		    forAll(fgmData_[0][0], k)
		    {
			is.read(fgmData_[i][j][k]);
		    }
		    readToNewline(is);
		}
		fgmCV1_[i] = rawCV1;
	    }

/*
    O ERRO DE LEITURA ESTA ANTES DESTA PARTE, POIS NAS PROXIMAS LINHAS,QUANDO ELE ESCREVE
    NA TELA O VALOR MAXIMO DAS VARIAVEIS DE CONTROLE, ELE ESCREVE COMO SENDO O VALOR 
    MAXIMO DE PV (CV2) = 2.XXE-10.. LOGO, NA LEITURA OCORRE ALGUMA COISA ESTRANHA..
*/
    
    //- Find max. value of CV 1 & 2; Needed to avoid segmentation faults if the 
    //  CV exceeds the tabulated values
    maxCV1_ = max(fgmCV1_);
    maxCV2_ = max(fgmCV2_);
    Info << "               Max. value of Control Variable 1: " << maxCV1_ << nl
         << "               Max. value of Control Variable 2: " << maxCV2_ << nl
         << endl;
}


/*---------------------------------------------------------------------------*\
MINIMUM- AND MAXIMUM VALUES FOR REACTION PROGRESS VARIABLE

[NUMBER_FLAMELETS]
 251
[DATA]
  <MixtureFraction>   <PVmin>   <PVmax>
\*---------------------------------------------------------------------------*/
void Foam::fanzyLookUp::readScalePVFile(const fileName inputName)
{
    IFstream is(inputName);
    label nFlamelets;
    
    if (is.good())
    {
        Info<< "fanzyLookUp: Read ScalePV file..." << endl;
    }
    else
    {
        FatalErrorIn
        (
            "fanzyLookUp::readScalePVFile(const fileName& inputName)"
        )
        << "cannot read " << is.name() << abort(FatalError);
    }

    //- Read header
    readToNewline(is);                       // MINIMUM- AND MAXIMUM VALUES FOR REACTION PROGRESS VARIABLE
    readToNewline(is);                       // 
    readToNewline(is);                       // [NUMBER_FLAMELETS]
    is.read(nFlamelets); readToNewline(is);  // <nFlamelets>
    readToNewline(is);                       // [DATA]
    
    /* Initialize global 2D array
     * ==========================
     * [ Dim 1 x        Dim 2  ] = 
     * [     3 x (Size of CV1) ]
     */
    scalePVTable_ = List<List<scalar> >(3,List<scalar>(nFlamelets));

    //- Check dimensions matching of CV1 in FGM-file and ScalePV-file
    if (nFlamelets == nFgmCV1_)
    {
        Info << "                   Number of flamelets: " << nFlamelets << nl
             << endl;
    }
    else
    {
        FatalErrorIn
        (
            "fanzyLookUp::readScalePVFile(const fileName& inputName)"
        )
        << "Dimension mismatch of " << is.name() << " and " << fgmFile_ << " ."
        << abort(FatalError);
    }

    //- Read scaling data
    forAll(scalePVTable_[0], i)
    {
        forAll(scalePVTable_, j)
        {
            is.read(scalePVTable_[j][i]);
        }
        readToNewline(is);
    }
}


Foam::scalar Foam::fanzyLookUp::scalePV
(
    const scalar CV1,
    const scalar CV2
) const
{
    //- Find lower index for CV1 and set the higher index
    label lo  = findLower(scalePVTable_[0], CV1), hi = lo+1;
    scalar CV2min, CV2max;
    if (scalePVTable_[0][hi] == CV1)
    {
        //- If we have a tabulated value no interpolation...
        CV2min = scalePVTable_[1][hi];
        CV2max = scalePVTable_[2][hi];
    }
    else
    {
        //- ... else do linear interpolation.
        CV2min = scalePVTable_[1][lo]
             + (scalePVTable_[1][hi]-scalePVTable_[1][lo])
             / (scalePVTable_[0][hi]-scalePVTable_[0][lo])
             * (CV1                  -scalePVTable_[0][lo]);
        CV2max = scalePVTable_[2][lo]
             + (scalePVTable_[2][hi]-scalePVTable_[2][lo])
             / (scalePVTable_[0][hi]-scalePVTable_[0][lo])
             * (CV1                  -scalePVTable_[0][lo]);
    }
    //- scaledPV = (unscaledPV - PVmin(Z)) / (PVmax(Z) - PVmin(Z))
    return min((CV2 - CV2min) / max(CV2max - CV2min, SMALL), scalar(1));
}


Foam::tmp<Foam::volScalarField> Foam::fanzyLookUp::getScaledPV
(
    const volScalarField& foamCV1,
    const volScalarField& foamCV2
) const
{
    tmp<volScalarField> tScaledPV
    (
        new volScalarField
        (
            IOobject
            (
                "scaledPV",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    volScalarField& scaledPV = tScaledPV();

    forAll(foamCV1, cellI)
    {
        scaledPV[cellI] = scalePV(foamCV1[cellI],foamCV2[cellI]);
    }

    forAll(foamCV1.boundaryField(), patchi)
    {
        const fvPatchScalarField& pFoamCV1 = foamCV1.boundaryField()[patchi];
        fvPatchScalarField& pScaledPV = scaledPV.boundaryField()[patchi];

        forAll(pFoamCV1, faceI)
        {
            pScaledPV[faceI] = scalePV(foamCV1[faceI],foamCV2[faceI]);
        }
    }

    return tScaledPV;
}


Foam::tmp<Foam::volScalarField> Foam::fanzyLookUp::sourcePV2D
(
    const volScalarField& foamCV1,
    const volScalarField& foamCV2
) const
{
    tmp<volScalarField> tsourcePV
    (
        new volScalarField
        (
            IOobject
            (
                "sourcePV",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0.0)
        )
    );

    volScalarField& sourcePV = tsourcePV();

    forAll(sourcePV, celli)
    {
        sourcePV[celli] = getValue2D(foamCV1[celli],foamCV2[celli],sourcePVindex_);
    }

    return tsourcePV;
}


Foam::tmp<Foam::volScalarField> Foam::fanzyLookUp::getField2D
(
    const volScalarField& foamCV1,
    const volScalarField& foamCV2,
    const label varI
) const
{
    tmp<volScalarField> tFgmField
    (
        new volScalarField
        (
            IOobject
            (
                "fgm_" + fgmVariableNames_[varI+2],
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    volScalarField& fgmField = tFgmField();

    forAll(foamCV1, celli)
    {
        fgmField[celli] = getValue2D(foamCV1[celli],foamCV2[celli],varI);
    }

    forAll(foamCV1.boundaryField(), patchi)
    {
        const fvPatchScalarField& pFoamCV1 = foamCV1.boundaryField()[patchi];
        const fvPatchScalarField& pFoamCV2 = foamCV2.boundaryField()[patchi];
        fvPatchScalarField& pFgmField = fgmField.boundaryField()[patchi];
        
        forAll(pFoamCV1, facei)
        {
            pFgmField[facei] = getValue2D(pFoamCV1[facei],pFoamCV2[facei],varI);
        }
    }

    return tFgmField;
}


Foam::scalar Foam::fanzyLookUp::getValue2D
(
    const scalar foamCV1,
    const scalar foamCV2,
    const label varI
) const
{
    if(scalePV_)
    {
        scalar foamCV2s = scalePV(foamCV1, foamCV2);
        return interpolate2D
               (
                   min(foamCV1,maxCV1_),
                   min(foamCV2s,maxCV2_),
                   varI
               );
    }
    else
    {
        return interpolate2D
               (
                   min(foamCV1,maxCV1_),
                   min(foamCV2,maxCV2_),
                   varI
               );
    }
}


Foam::scalar Foam::fanzyLookUp::interpolate2D
(
    const scalar foamCV1,
    const scalar foamCV2,
    const label varI
) const
{
    scalar interpolatedValue;
    //- Get index of the lower value
    label loCV1  = findLower(fgmCV1_, foamCV1);
    label loCV2  = findLower(fgmCV2_, foamCV2);
    
    //- Set index of the higher value
    label hiCV1  = loCV1+1;
    label hiCV2  = loCV2+1;

    /* Bi-linear interpolation
     * =======================
     *  implemented as in: 
     *  Numerical Recipes in Fortran 77, Second Edition (1992),
     *  '3.6 Interpolation in Two or More Dimensions' (p.117)
     */
    scalar t, u;
    
    //- Check if 'findLower' returned -1, i.e. the Control Variable is out of 
    //  range at the lower end.
    if (loCV1 == -1)
    {
        loCV1 = hiCV1;
        t = 1.0;
    }
    else
    {
        t = (foamCV1 - fgmCV1_[loCV1])/(fgmCV1_[hiCV1] - fgmCV1_[loCV1]);
    }
    if (loCV2 == -1)
    {
        loCV2 = hiCV2;
        u = 1.0;
    }
    else
    {
        u = (foamCV2 - fgmCV2_[loCV2])/(fgmCV2_[hiCV2] - fgmCV2_[loCV2]);
    }

    //- Calculate interpolated value for given control variable values
    interpolatedValue = (1-t) * (1-u) * fgmData_[loCV1][loCV2][varI]
                      +    t  * (1-u) * fgmData_[hiCV1][loCV2][varI]
                      +    t  *    u  * fgmData_[hiCV1][hiCV2][varI]
                      + (1-t) *    u  * fgmData_[loCV1][hiCV2][varI];

    return interpolatedValue;
}


// ************************************************************************* //
