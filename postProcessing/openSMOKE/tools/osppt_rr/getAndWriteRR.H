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
    kineticsMapXML->GetReactionRates(&qCell); //reaction rates already computed!

    kineticsMapXML->FormationRates(&RCell);

    //----added for osppt_rr----//
    for(unsigned int AUXJ=0;AUXJ<nr;AUXJ++)
        rReactions[AUXJ].ref()[celli] = qCell[AUXJ+1]; //kmol/m3/s
    //----//

}

//----added for osppt_rr----//
for (unsigned int AUXJ=0;AUXJ<nr;AUXJ++)
{
    rReactions[AUXJ].correctBoundaryConditions();

    rReactions[AUXJ].write();
}
//----//
