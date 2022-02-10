bool calcChemScales(OSppOptionsDictionary.lookupOrDefault("calcChemScales",true));


if(calcChemScales)
{
    Info << "----------  calcChemScales  ------------    "<< calcChemScales<<endl;

    threshold_temperature_chemical_time = readScalar(coeffsDictChemTimes.lookup("ChemicalTimesTemperatureThreshold"));
    threshold_chemical_time = readScalar(coeffsDictChemTimes.lookup("ChemicalTimesThreshold"));
    word chemicalTimesTypeString(coeffsDictChemTimes.lookup("ChemicalTimesType"));

    if (chemicalTimesTypeString == "formationRates")
    	characteristic_chemical_times_type = CharacteristicChemicalTimes::CHEMICAL_TIMES_FORMATION_RATES;
    else if (chemicalTimesTypeString == "reactionRates")
    	characteristic_chemical_times_type = CharacteristicChemicalTimes::CHEMICAL_TIMES_REACTION_RATES;
    else if (chemicalTimesTypeString == "eigenvalues")
    	characteristic_chemical_times_type = CharacteristicChemicalTimes::CHEMICAL_TIMES_EIGENVALUES;
    else
    {
    	Info << "Wrong chemical times type: formationRates | reactionRates | eigenvalues" << endl;
    	abort();
    }

    //minimum_Da  = readScalar(coeffsDictChemTimes.lookup("minimum_Da"));
}
