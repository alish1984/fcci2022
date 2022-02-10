Info <<" Create OSpp fields"<<endl;
//declare placeholders (pointers)
PtrList<volScalarField> xSpecies;
PtrList<volScalarField> omegaDotSpecies;
PtrList<volScalarField> omegaDotPosSpecies;
PtrList<volScalarField> omegaDotNegSpecies;

{
	// Check the size of the mechanism and order of species
	if (thermodynamicsMapXML->NumberOfSpecies() != Y.size())
	FatalError << "Kinetic schemes do not match..." << endl;

	for (unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
	if (Y[i].name() != thermodynamicsMapXML->NamesOfSpecies()[i])
	{
		FatalError << "Kinetic schemes do not match..." << endl;
	}

    //resize the pointers in the first dimension to number of species
	omegaDotPosSpecies.resize(ns);
    omegaDotSpecies.resize(ns);
    omegaDotNegSpecies.resize(ns);
    xSpecies.resize(ns);

    //resize the pointers in the second dimension to number of computational cells
	Info<<"setting omegaDotSpecies"<<endl;
	forAll(Y,i)
	{
		omegaDotPosSpecies.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"post_OmegaDotPos_"+ Y[i].name(),
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
                mesh,
				dimensionedScalar("zero", dimensionSet(0, -3, -1, 0, 1), 0.), //kmol/m3/s
                extrapolatedCalculatedFvPatchScalarField::typeName
			)
		);
        omegaDotSpecies.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"post_OmegaDot_"+ Y[i].name(),
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
                mesh,
				dimensionedScalar("zero", dimensionSet(0, -3, -1, 0, 1), 0.), //kmol/m3/s
                extrapolatedCalculatedFvPatchScalarField::typeName
			)
		);
        omegaDotNegSpecies.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"post_OmegaDotNeg_"+ Y[i].name(),
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
                mesh,
				//dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0.), //kg/m3/s
                dimensionedScalar("zero", dimensionSet(0, -3, -1, 0, 1), 0.), //kmol/m3/s
                extrapolatedCalculatedFvPatchScalarField::typeName
			)
		);
        xSpecies.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"post_X_"+ Y[i].name(),
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
                mesh,
				dimensionedScalar(dimless, Zero),
                extrapolatedCalculatedFvPatchScalarField::typeName
			)
		);
	}
}
Info <<" End of Create OSpp fields"<<endl;
