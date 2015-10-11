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

    label i,k1,k2;
    scalar cv2ns,cv2min_low,cv2max_low,cv2min_high,cv2max_high,dh_low,dh_high,b1,b2,d1,d2,e1,e2,cv1,cv2ns2,f10,f11,f30,f31,f32,f33;
    	
    cv1=MNMX(foamCV1,fgmCV1_[0],fgmCV1_[(nFgmCV1_*nFgmCV2_)-1]);
    i=MNMX((nFgmCV1_-1.)*cv1,0,nFgmCV1_-2.); 			/* lower index */
    b1=i/(nFgmCV1_-1.);
    b2=(i+1.)/(nFgmCV1_-1.);
    f10=MNMX(((b2-cv1)/(b2-b1)),0.,1.);
    f11=1.-f10;
    
    cv2min_low=fgmCV2_[((i*nFgmCV2_)+1)-1];
    cv2max_low=fgmCV2_[((i+1)*nFgmCV2_)-1];
    cv2ns=MNMX(foamCV2,cv2min_low,cv2max_low);
    dh_low=((cv2max_low-cv2min_low)/(nFgmCV2_-1));
    k1= MNMX(cv2ns/dh_low,0.0,nFgmCV2_-2);
    d1=k1*dh_low + cv2min_low;
    d2=(k1+1)*dh_low+cv2min_low;
    f30=MNMX(((d2-cv2ns)/(d2-d1)),0.,1.);
    f31=1.-f30;
    
    cv2min_high=fgmCV2_[(((i+1)*nFgmCV2_)+1)-1];
    cv2max_high=fgmCV2_[((i+2)*nFgmCV2_)-1];
    cv2ns2=MNMX(foamCV2,cv2min_high,cv2max_high);
    dh_high=((cv2max_high-cv2min_high)/(nFgmCV2_-1));
    k2= MNMX(cv2ns2/dh_high,0.0,nFgmCV2_-2);
    e1=k2*dh_high + cv2min_high;
    e2=(k2+1)*dh_high+cv2min_high;
    f32=MNMX(((e2-cv2ns2)/(e2-e1)),0.,1.); 
    f33=1.-f32; 
      
      
      
      /*IN THE CASE OF REACHING THE UPPER LIMIT OF THE MANIFOLD BOTH FACTOR MUST BE SET TO THE UPPER LIMIT:*/
      /**/
      if ((k1==nFgmCV2_-2) && (f31==1) && (k2!=nFgmCV2_-2) && (f33!=1))
      {
          k2=nFgmCV2_-2; f32=0; f33=1;
      }
      /**/
      else if ((k2==nFgmCV2_-2) && (f33==1) && (k1!=nFgmCV2_-2) && (f31!=1))
      {
          k1=nFgmCV2_-2; f30=0; f31=1;
      }
      /**/
      /*IN THE CASE OF REACHING THE LOWER LIMIT OF THE MANIFOLD BOTH FACTOR MUST BE SET TO THE LOWER LIMIT:*/
      /**/
      if ((k1==0) && (f30==1) && (k2!=0) && (f32!=1))
      {
          k2=0; f32=1;f33=0;
      }
      else if ((k2==0) && (f32==1) && (k1!=0) && (f30!=1))
      {
         k1=0; f30=1; f31=0;
      }
      
      
      interpolatedValue  =  f10*f30*fgmData_[ k1   + i   *nFgmCV2_][varI]
    		       +f10*f31*fgmData_[(k1+1)+ i   *nFgmCV2_][varI]
    		       +f11*f32*fgmData_[ k2   +(i+1)*nFgmCV2_][varI]
    		       +f11*f33*fgmData_[(k2+1)+(i+1)*nFgmCV2_][varI];
    
    return interpolatedValue; 

}


// ************************************************************************* //
