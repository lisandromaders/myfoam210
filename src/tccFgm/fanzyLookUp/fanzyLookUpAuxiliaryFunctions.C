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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::fanzyLookUp::writeFgmFile(const fileName outputName)
{
    //- Initialize output stream
    OFstream os(outputName);
    
    //- Set the width of field for integer and floats
    label intWidth = 4, floatWidth = 22;
    
    //- Set the write precision
    os.precision(14);
    
    //- Set the write format (scientific: 0.123e+00)
    os.setf(std::ios::scientific); 

    Info<< "fanzyLookUp: Write out FGM table to: " << outputName << nl
        << endl;
    
    //- Write file header
    os  <<
        "/*--------------------------------*- C++ "
        "-*----------------------------------*\\\n"
        "|                  Do not use as input-file for simulations!"
        "                  |\n"
        "|            Only intended to verify the reading of the FGM table!"
        "            |\n"
        "\\*-----------------------------------------"
        "----------------------------------*/\n";
    //- Write FGM header
    os << "FLAMELET GENERATED MANIFOLD" << nl
       << nl
       << "[NUMBER_FLAMELETS]" << nl;
    os.width(intWidth);
    os << nFgmCV1_ << nl
       << "[NUMBER_PV]" << nl;
    os.width(intWidth);
    os << nFgmCV2_ << nl
       << "[NUMBER_VARIABLES]" << nl;
    os.width(intWidth);
    os << nFgmData_+2 << nl
       << nl
       << "[DATA]" << nl;
    forAll(fgmVariableNames_, i)
    {
        os << fgmVariableNames_[i] << ' ';
    }
    os << nl;
    
    //- Write FGM data
    label dataLabel = 0;
    forAll(fgmData_, i)
    {
        forAll(fgmData_[0], j)
        {
            os.width(floatWidth);
            os << fgmCV1_[i] << ' ';
            os.width(floatWidth);
//            os << fgmCV2_[j] << ' ';
            os << rawCV2[i][j] << ' ';
            
            forAll(fgmData_[0][0], k)
            {
                os.width(floatWidth);
                os << (fgmData_[i][j][k]) << ' ';
            }
            os << nl;
            dataLabel++;
        }
    }
    os.flush();
}


void Foam::fanzyLookUp::writeScalePVFile(const fileName outputName)
{
    //- Initialize output stream
    OFstream os(outputName);
    
    //- Set the width of field for integer and floats
    label intWidth = 4, floatWidth = 22;
    
    //- Set the write precision
    os.precision(14);
    
    //- Set the write format (scientific: 0.123e+00)
    os.setf(std::ios::scientific); 

    Info<< "fanzyLookUp: Write out ScalePV table to: " 
        << outputName << nl << endl;
    
    //- Write file header
    os  <<
        "/*--------------------------------*- C++ "
        "-*----------------------------------*\\\n"
        "|                  Do not use as input-file for simulations!"
        "                  |\n"
        "|         Only intended to verify the reading of the ScalePV table!"
        "               |\n"
        "\\*-----------------------------------------"
        "----------------------------------*/\n";
    //- Write ScalePV header
    os << "MINIMUM- AND MAXIMUM VALUES FOR REACTION PROGRESS VARIABLE" << nl
       << nl
       << "[NUMBER_FLAMELETS]" << nl;
    os.width(intWidth);
    os << scalePVTable_[0].size() << nl
       << "[DATA]" << nl;
    
    //- Write ScalePV data
    forAll(scalePVTable_[0], i)
    {
        os.width(floatWidth);
        os << scalePVTable_[0][i] << ' ';
        os.width(floatWidth);
        os << scalePVTable_[1][i] << ' ';
        os.width(floatWidth);
        os << scalePVTable_[2][i]
           << nl;
    }
    os.flush();
}


void Foam::fanzyLookUp::writeFgmFileOpenfoam()
{
    Info<< "fanzyLookUp: Write out FGM table in OpenFOAM format." << endl;
    
    //- FGM data directory
    fileName fgmDataDir = "FGMData";

    //- Create directory if does not exist.
    mkDir(fgmDataDir);

    //- Initialize output streams
    OFstream osCV1(fgmDataDir/fgmVariableNames_[0]);
    OFstream osCV2(fgmDataDir/fgmVariableNames_[1]);
    OFstream osData(fgmDataDir/"fgmData");

    //- Write FGM data
    osCV1 << fgmCV1_;
    osCV2 << fgmCV2_;
    osData << fgmData_;
}


// ************************************************************************* //
