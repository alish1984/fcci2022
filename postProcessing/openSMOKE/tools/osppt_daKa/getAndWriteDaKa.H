scalarField& tauChemCells = tauChem.ref();
scalarField& DaEtaCells = DaEta.ref();
scalarField& DaLCells = DaL.ref();
scalarField& KaCells = Ka.ref();

tmp<volScalarField> tmpnu(thermo.nu());
const volScalarField& nu = tmpnu();

tmp<volScalarField> tmpalpha(thermo.alpha());
const volScalarField& alpha = tmpalpha();

tmp<volScalarField> tmprho(thermo.rho());
const volScalarField& rho = tmprho();

chi=2*alpha/rho*magSqr(fvc::grad(Z) );
chi.correctBoundaryConditions();

forAll(T,celli) // loop over all cells and fill the values
{
    // ----- NOTE that OpenSMOKE::OpenSMOKEVectorDouble objects are "1" based objects ----- //
    OpenSMOKE::OpenSMOKEVectorDouble omegaCell(ns);
    // OpenSMOKE::OpenSMOKEVectorDouble xCell(ns);
    // OpenSMOKE::OpenSMOKEVectorDouble cCell(ns);
    // OpenSMOKE::OpenSMOKEVectorDouble RCell(ns);
    // OpenSMOKE::OpenSMOKEVectorDouble qCell(nr);

    //double hMean, MWMean;

    for(unsigned int i=0;i<ns;i++)
    {
             omegaCell[i+1] = Y[i].internalField()[celli];
    }

    //thermodynamicsMapXML->MoleFractions_From_MassFractions(xCell, MWMean, omegaCell);

    double pCell=p[celli];
    const double TCell = T.internalField()[celli];
    // double MWCell=MWMean;
    // double cToTCell = pCell/(PhysicalConstants::R_J_kmol * TCell);
    // double rhoCell = cToTCell*MWCell;
    // Product(cToTCell, xCell, &cCell);

    thermodynamicsMapXML->SetTemperature(TCell);
    thermodynamicsMapXML->SetPressure(pCell);
    kineticsMapXML->SetTemperature(TCell);
    kineticsMapXML->SetPressure(pCell);

    // kineticsMapXML->ReactionEnthalpiesAndEntropies();
    // kineticsMapXML->KineticConstants();
    // kineticsMapXML->ReactionRates(cCell);
    // kineticsMapXML->GetReactionRates(&qCell);
    //
    // kineticsMapXML->FormationRates(&RCell);


    // get Chemical time from OpenSMOKE
    const scalar tauChemCell = chemicalTimes->CalculateCharacteristicChemicalTime(TCell,pCell,omegaCell);
    tauChemCells[celli] = tauChemCell; // fill
    //

    const scalar nuCell = nu[celli];
    const scalar kCell = k[celli];
    const scalar epsCell = epsilon[celli];
    const scalar chiCell = chi[celli];

    const scalar tauEtaCell = Foam::sqrt(nuCell/(epsCell + SMALL));
    const scalar tauLCell = kCell/(epsCell + SMALL);

    scalar DaEtaCell = tauEtaCell/tauChemCell; //evaluate Da
    DaEtaCell = min(DaEtaCell, maximum_Da); // limit minimum
    DaEtaCell = max(DaEtaCell, minimum_Da); // limit maximum
    DaEtaCells[celli] = DaEtaCell; // fill

    scalar DaLCell = tauLCell/tauChemCell; //evaluate Da
    DaLCell = min(DaLCell, maximum_Da); // limit minimum
    DaLCell = max(DaLCell, minimum_Da); // limit maximum
    DaLCells[celli] = DaLCell; // fill

    scalar KaCell = tauChemCell/tauEtaCell; //evaluate Ka
    KaCell = min(KaCell, maximum_Ka); // limit minimum
    KaCell = max(KaCell, minimum_Ka); // limit maximum
    KaCells[celli] = KaCell; // fill
}

tauChem.correctBoundaryConditions();
DaL.correctBoundaryConditions();
DaEta.correctBoundaryConditions();
Ka.correctBoundaryConditions();
chi.correctBoundaryConditions();

tauChem.write();
DaL.write();
DaEta.write();
Ka.write();
chi.write();
