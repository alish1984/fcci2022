/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    osppt_daKa

Description
    Solver for combustion with chemical reactions.

\*---------------------------------------------------------------------------*/

#include "openSMOKE_headers.H" //add this

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
/*      argList::addNote
    (
        "Solver for combustion with chemical reactions"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "rhoEqn.H"

        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "YEqn.H"
            #include "EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                if (pimple.consistent())
                {
                    #include "pcEqn.H"
                }
                else
                {
                    #include "pEqn.H"
                }
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        rho = thermo.rho();

        runTime.write();

        runTime.printExecutionTime(Info);
    }
*/
    argList::addNote //a note for your solver description
    (
        "An OpenSMOKE++ post processing tool for evaluation of species prod/consumption rates"
    );

    timeSelector::addOptions(); //options for command line
    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"

    Info<< nl << "Starting time loop" << endl;
    instantList timeDirs = timeSelector::select0(runTime, args); //select (or get) times

    #include "createMesh.H"

    forAll(timeDirs, timeI) //loop over all time directories
    {
        const double tStartFull = OpenSMOKE::OpenSMOKEGetCpuTime(); //define a timer!

        runTime.setTime(timeDirs[timeI], timeI); //set the current time

        Info<< "\n------------------------------------\n" << endl;
        Info<< "\n         Calculations for time      \n" << endl;
        Info<< "\n         "<<runTime.timeName()<<"   \n" << endl;
        Info<< "\n------------------------------------\n" << endl;

        #include "createFields.H" //main create fields

        //------- OpenSMOKE -------//
        #include "createOpenSMOKEFieldsGlobal.H" //create main openSMOKE fields
        #include "createOpenSMOKEFields_daKa.H" //create specific openSMOKE fields
        # include "getAndWriteDaKa.H" //evaluate (and write) species prod/consumption rates
        //------- OpenSMOKE -------//

        const double tEndFull = OpenSMOKE::OpenSMOKEGetCpuTime(); // stop the timer
        Info << "CPU time full processes : " << tEndFull - tStartFull << " s " << endl;
    }
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
