/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    particlesIdentificationFoam

Description
    This utility alylize alpha field of domain, identifies particles/droplets
    and writes them into a single text file
    Derived from interFoam
    Changed in createFields.H: gamma:NO_WRITE

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "incompressibleTwoPhaseMixture.H"
#include "interpolationTable.H"
#include "pimpleControl.H"
#include "preBubble.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createTimeControls.H"

pimpleControl pimple(mesh);

#include "initContinuityErrs.H"
#include "createFields.H"
#include "readTimeControls.H"
#include "createMyDict.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "\nStarting time loop\n" << endl;

/*
    const scalar& write_interval(
	  readScalar(runTime.controlDict().lookup("writeInterval"))
	  );

    runTime.setDeltaT
    (
        write_interval
    );
*/
Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

forAll(timeDirs, timei)
{
    // update the alpha field
    Info<< "Time = " << timeDirs[timei].value() << nl << endl;

    volScalarField alpha1
	(
	 IOobject
	 (
	  "alpha.fuel",
	  std::to_string(timeDirs[timei].value()),
	  mesh,
	  IOobject::READ_IF_PRESENT,
	  IOobject::NO_WRITE
	  ),
	 mesh
	 );
    

    // update velosity field
    volVectorField U
    (
        IOobject
        (
            "U",
            std::to_string(timeDirs[timei].value()),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh
    );

    //Variables qeuired for calculations

    scalar alphaLimMax = readScalar(myDict.lookup("alphaLimMax"));
    scalar alphaLimMin = readScalar(myDict.lookup("alphaLimMin"));



    // List which contains the cells for each particle

    // List of particles
    DynamicList<DynamicList<label>> Nn;
    //List of particles positions
    DynamicList<vector> partPos;
    //List of particles diameter
    DynamicList<scalar> partDiam;
    //List of particles Volume
    DynamicList<scalar> partVol;
    //List of particles velosity
    DynamicList<vector> partVel;
    //List of boundary points for particles
    DynamicList<DynamicList<vector>> partLimits;


    //initialization ID field with -1 value
    forAll (mesh.C(), celli) //loop through cell centres
    {
        ID[celli] = -1.0;
    }

    int ii = 1; // index which is the same for each cell into the same particles, starts with 1
    forAll (mesh.C(), celli) //loop through cell centres
    {
        
        if(alpha1[celli] > alphaLimMax && ID[celli] == -1)
        {
            DynamicList<label> neigh;
            neigh.append(celli);
            addNeighbours(neigh,celli,mesh,alpha1,alphaLimMin, ID, ii);
            Nn.append(neigh);
            ++ii;
        }
    }


    Info << "Number of particles = " << ii - 1  << endl;

    //Position of particle
    calcPos(Nn,mesh,partPos);
    //Info << "Position of particles = " << partPos  << endl;
    //Volume of the particle
    calcVolume(Nn,mesh,alpha1,partVol);
    //Info << "Volume of particles = " << partVol  << endl;
    //Diameters of particle
    calcDiam(partVol,partDiam);
    //Info << "Diameter of particles = " << partDiam  << endl;
    //Velocity of particle
    calcVel(Nn,mesh,U,partVel);
    //Info << "Velocity of particle = " << partVel  << endl;
    //Limits of particle
    //calcLim(Nn,mesh,partLimits);

    //Make a list of particle which shuld be added
    DynamicList<label> activeList = activeParticles(partPos);
    //Info << "Number of particle  = " << activeList.size()  << endl;
    //Write the results
    ID.write();
    #include "write.H"

}
}












































// ************************************************************************* //
