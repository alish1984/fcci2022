Info <<" Create OSpp fields"<<endl;
//declare placeholders (pointers)

volScalarField DaEta
(
	IOobject
	(
		"post_Da_eta",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar(dimless, Zero),
	extrapolatedCalculatedFvPatchScalarField::typeName
);

volScalarField DaL
(
	IOobject
	(
		"post_Da_L",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar(dimless, Zero),
	extrapolatedCalculatedFvPatchScalarField::typeName
);

volScalarField chi
(
    IOobject
    (
        "chi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("chi", dimensionSet(0,0,-1,0,0,0,0),0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField Ka
(
	IOobject
	(
		"post_Ka",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar(dimless, Zero),
	extrapolatedCalculatedFvPatchScalarField::typeName
);

volScalarField tauChem
(
    IOobject
    (
        "post_tauChem",
        runTime.timeName(),
        runTime,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
	dimensionedScalar(dimTime, Zero),
	extrapolatedCalculatedFvPatchScalarField::typeName
);

// characteristic chemical times objects
dictionary coeffsDictChemTimes(OSppOptionsDictionary.subDict("chemicalTimes"));
scalar threshold_temperature_chemical_time = 300.;
scalar threshold_chemical_time = 0.1;
const double minimum_Da = 1e-10;
const double maximum_Da = 10.;
const double minimum_Ka = 1e-10;
const double maximum_Ka = 100.;
CharacteristicChemicalTimes::CharacteristicChemicalTimesType characteristic_chemical_times_type = CharacteristicChemicalTimes::CHEMICAL_TIMES_FORMATION_RATES;

CharacteristicChemicalTimes* chemicalTimes = new CharacteristicChemicalTimes(*thermodynamicsMapXML, *kineticsMapXML);

#include "readOptions.H" //read options

chemicalTimes->SetThresholdChemicalTime(threshold_chemical_time);
chemicalTimes->SetType(characteristic_chemical_times_type);
chemicalTimes->SetThresholdTemperatureChemicalTime(threshold_temperature_chemical_time);

Info <<" End of Create OSpp fields"<<endl;
