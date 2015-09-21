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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::fanzyLookUp::mesh() const
{
    return mesh_;
}

inline Foam::Switch Foam::fanzyLookUp::debugOutput() const
{
    return debugOutput_;
}

inline Foam::List<Foam::word> Foam::fanzyLookUp::fgmVariableNames() const
{
    return fgmVariableNames_;
}

inline Foam::List<Foam::word> Foam::fanzyLookUp::cvNames() const
{
    return SubList<word>(fgmVariableNames_,2,0);
}

inline Foam::List<Foam::word> Foam::fanzyLookUp::dataNames() const
{
    return SubList<word>(fgmVariableNames_,nFgmData_,2);
}

// ************************************************************************* //
