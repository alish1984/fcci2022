Info <<" Create OSpp fields"<<endl;
//declare placeholders (pointers)
PtrList<volScalarField> rReactions; //----added for osppt_rr----//
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
	rReactions.resize(nr); //----added for osppt_rr----//

    //resize the pointers in the second dimension to number of computational cells
	//----added for osppt_rr----//
	Info<<"setting rates of reactions "<<endl;
	for (unsigned int i=0;i<nr;i++)
	{
		rReactions.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"post_RReaction_"+ std::to_string(i),
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
	}
}
Info <<" End of Create OSpp fields"<<endl;
