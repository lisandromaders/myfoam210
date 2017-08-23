/*
Library that reads and interpolates the manifold
*/

#include "fanzyLookUp.H"
#include "fvMesh.H"
#include "Time.H"
#include "ListOps.H"

#define MN(a,c) ((a)<(c)?(a):(c))
#define MX(a,c) ((a)>(c)?(a):(c))
#define MNMX(a,b,c) ((a)<(b)?(b):((a)>(c)?(c):(a)))
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
    sourcePVindex_(readLabel(lookup("sourcePVindex")))
{
    Info<< "fanzyLookUp: FGM file = " << fgmFile_ << endl;
    // Read and write the FGM table on construction (class instantiated)
    readFgmFile(fgmFile_);
    fileName fgmTableOut = path()/"FGM.out";
//     fileName fgmTableOut = "FGMTable.out";
    writeFgmFile(fgmTableOut);
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
    nFgmData_ = nVariables;
    fgmCV1_ = List<scalar>
		(
		    nFgmCV1_*nFgmCV2_
		);
    fgmCV2_ = List<scalar>
		(
		    nFgmCV1_*nFgmCV2_
		);

    fgmData_  = List<List<scalar> > 
                (
                    nFgmCV1_*nFgmCV2_,
                    List<scalar> 
                    (
                        nFgmData_
                    )
                );

    //- Write out some information
    Info << "              Number of lines of Manifold : " 
         << fgmData_.size() << nl
         << "              Number of colums - variables (including CV1 and CV2): " 
         << fgmData_[0].size() << nl
         << endl;
            
    Info << "               fgmVariableNames[" << 0 << "] (cv1) = " 
         << fgmVariableNames_[0] << nl
         << "               fgmVariableNames[" << 1 << "] (cv2) = " 
         << fgmVariableNames_[1] << nl;
    for (label i=0; i<nFgmData_-2; i++)
    {
        Info << "               fgmVariableNames[" << i+2 << "] (dataField[" 
             << i+2 << "]) = " << fgmVariableNames_[i+2] << nl;
    }
    Info << endl;

/*
	O codigo original foi escrito para buscar os dados em uma configuracao
	especifica de manifold. O algoritmo de busca foi modificado para lidar 
	com o manifold que o Zimmer e o Cristian me passaram... Em caso de 
	duvidas, contatar um de nos para pegar o manifold usado
*/
 
    for (label i=0; i<nFlamelets*nPV; i++)
    {
	for (label j=0; j<nFgmData_; j++)
        {
	     is.read(fgmData_[i][j]);
	}
    }

/*
	Vou armazenar os valores de cv1 e cv2 em arrays 1D originais para se usar nas
	outras funcoes.. acho que vai facilitar. Provavelmente a forma que eu fiz nao eh
	a mais otimizada em termos de tempo computacional, mas nao vai alterar muita coisa
	uma vez que a leitura eh feita somente nesse momento..
*/

    forAll(fgmData_, i)
    {
	fgmCV1_[i] = fgmData_[i][0];
	fgmCV2_[i] = fgmData_[i][1];
    }

    //- Find max. value of CV 1 & 2; Needed to avoid segmentation faults if the 
    //  CV exceeds the tabulated values

    maxCV1_ = max(fgmCV1_);
    maxCV2_ = max(fgmCV2_);
    Info << "               Max. value of Control Variable 1: " << maxCV1_ << nl
         << "               Max. value of Control Variable 2: " << maxCV2_ << nl
         << endl;
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
//                "fgm_" + fgmVariableNames_[varI+2],
                "fgm_" + fgmVariableNames_[varI],
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
    return interpolate2D
           (
               min(foamCV1,maxCV1_),
               min(foamCV2,maxCV2_),
               varI
           );
}


Foam::scalar Foam::fanzyLookUp::interpolate2D
(
    const scalar foamCV1,
    const scalar foamCV2,
    const label varI
) const
{
    scalar interpolatedValue;
    
    label i1m, i2m;
    scalar cv2min_low, cv2max_low, cv2min_high, cv2max_high, w1p, w1m, w2p, w2m;
    scalar eta1, eta2;
    scalar xmin, xmax;

    /* ---------  CV1 ----------*/
    xmin = fgmData_[0][0];
    xmax = fgmData_[(nFgmCV1_*nFgmCV2_) - 1][0];
    eta1 = (foamCV1 - xmin)/(xmax - xmin);
    eta1 = MX(0.0,eta1);
    eta1 = MN(1.0,eta1);
    /* -- lower index -- */
    i1m = eta1*(nFgmCV1_ - 1);
    i1m = MX(0, i1m);
    i1m = MN(nFgmCV1_ - 2,i1m);
    /* -- weigh factors -- */
    w1p = (eta1*(nFgmCV1_ - 1.0)) - i1m;
    w1m = 1.0 - w1p;

    /* --------- CV2 ----------*/
    cv2min_low =  fgmData_[((i1m * nFgmCV2_) +   1.0)    -1.0][1];
    cv2max_low =  fgmData_[((i1m +  1.0    ) * nFgmCV2_) -1.0][1];
    cv2min_high = fgmData_[((i1m +  1.0    ) * nFgmCV2_)     ][1];
    cv2max_high = fgmData_[((i1m +  2.0    ) * nFgmCV2_) -1.0][1];

    xmin = w1p*cv2min_low + w1m*cv2min_high; 
    xmax = w1p*cv2max_low + w1m*cv2max_high; 
    eta2 = (foamCV2 - xmin)/(xmax - xmin); 
    eta2 = MX(0.0,eta2);
    eta2 = MN(1.0,eta2);
    /* -- lower index -- */
    i2m = eta2*(nFgmCV2_ - 1);
    i2m = MX(0,i2m);
    i2m = MN(nFgmCV2_-2,i2m);
    /* -- weigh factors -- */
    w2p = (eta2*(nFgmCV2_ - 1.0)) - i2m;
    w2m = 1.0 - w2p;

    interpolatedValue =   w1m * w2m * fgmData_[ i2m    +  i1m    * nFgmCV2_][varI]
			+ w1m * w2p * fgmData_[(i2m+1) +  i1m    * nFgmCV2_][varI]
			+ w1p * w2m * fgmData_[ i2m    + (i1m+1) * nFgmCV2_][varI]
			+ w1p * w2p * fgmData_[(i2m+1) + (i1m+1) * nFgmCV2_][varI];

    return interpolatedValue;


}


// ************************************************************************* //
