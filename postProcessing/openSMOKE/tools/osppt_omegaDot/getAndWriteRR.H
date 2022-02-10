forAll(T,celli) // loop over all cells and fill the values
{
    // ----- NOTE that OpenSMOKE::OpenSMOKEVectorDouble objects are "1" based objects ----- //
    OpenSMOKE::OpenSMOKEVectorDouble omegaCell(ns);
    OpenSMOKE::OpenSMOKEVectorDouble xCell(ns);
    OpenSMOKE::OpenSMOKEVectorDouble cCell(ns);
    OpenSMOKE::OpenSMOKEVectorDouble RCell(ns);
    OpenSMOKE::OpenSMOKEVectorDouble qCell(nr);

    double hMean, MWMean;

    for(unsigned int i=0;i<ns;i++)
    {
             omegaCell[i+1] = Y[i].internalField()[celli];
    }

    thermodynamicsMapXML->MoleFractions_From_MassFractions(xCell, MWMean, omegaCell);

    double pCell=p[celli];
    const double TCell = T.internalField()[celli];
    double MWCell=MWMean;
    double cToTCell = pCell/(PhysicalConstants::R_J_kmol * TCell);
    double rhoCell = cToTCell*MWCell;
    Product(cToTCell, xCell, &cCell);

    thermodynamicsMapXML->SetTemperature(TCell);
    thermodynamicsMapXML->SetPressure(pCell);
    kineticsMapXML->SetTemperature(TCell);
    kineticsMapXML->SetPressure(pCell);

    kineticsMapXML->ReactionEnthalpiesAndEntropies();
    kineticsMapXML->KineticConstants();
    kineticsMapXML->ReactionRates(cCell);
    kineticsMapXML->GetReactionRates(&qCell);

    kineticsMapXML->FormationRates(&RCell);

    for(unsigned int AUXK=0;AUXK<ns;AUXK++)
        omegaDotSpecies[AUXK].ref()[celli] = RCell[AUXK+1]; //kmol/m3/s
        //omegaDotSpecies[AUXK].ref()[celli] = RCell[AUXK+1]* thermodynamicsMapXML_.MW()[AUXK+1];//if kg/m3/s

    // Ropa
    OpenSMOKE::ROPA_Data ropa;
    kineticsMapXML->RateOfProductionAnalysis(ropa);
    for (unsigned int AUXK=0;AUXK<ns;AUXK++)
    {
        for(unsigned int k=0;k<ropa.production_coefficients[AUXK].size();k++)
        {
            omegaDotPosSpecies[AUXK].ref()[celli] += ropa.production_coefficients[AUXK][k];
        }

        for(unsigned int k=0;k<ropa.destruction_coefficients[AUXK].size();k++)
        {
            omegaDotNegSpecies[AUXK].ref()[celli] += ropa.destruction_coefficients[AUXK][k];
        }
    }

    for (unsigned int AUXK=0;AUXK<ns;AUXK++)
    {
        omegaDotNegSpecies[AUXK].ref()[celli]= -1.*omegaDotNegSpecies[AUXK].ref()[celli];
        xSpecies[AUXK].ref()[celli]=xCell[AUXK+1];
    }
}
for (unsigned int AUXK=0;AUXK<ns;AUXK++)
{
    omegaDotNegSpecies[AUXK].correctBoundaryConditions();
    omegaDotPosSpecies[AUXK].correctBoundaryConditions();
    omegaDotSpecies[AUXK].correctBoundaryConditions();
    xSpecies[AUXK].correctBoundaryConditions();

    omegaDotNegSpecies[AUXK].write();
    omegaDotPosSpecies[AUXK].write();
    omegaDotSpecies[AUXK].write();
    xSpecies[AUXK].write();
}
