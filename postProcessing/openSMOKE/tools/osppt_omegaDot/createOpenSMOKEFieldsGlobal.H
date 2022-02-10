const word OSppOptionsDictionaryName("OSppOptions");

Info<< "Reading OSppOptions dictionary\n" << endl;
IOdictionary OSppOptionsDictionary
(
	IOobject
	(
		OSppOptionsDictionaryName,
		runTime.constant(),
		runTime.db(),
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

// Read the kinetic scheme in XML format
OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>* thermodynamicsMapXML;
OpenSMOKE::KineticsMap_CHEMKIN<double>* kineticsMapXML;
//OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>* transportPropertiesMapXML;
{
	word kinetics(OSppOptionsDictionary.lookup("kineticMechanism"));
    Info <<"kinetics  =  " <<  kinetics<<endl;
	boost::filesystem::path path_kinetics = kinetics;
    Info <<"-----------1"<<endl;
    std::cout <<"path_kinetics   =  "<<path_kinetics <<endl;

	rapidxml::xml_document<> doc;
	std::vector<char> xml_string;
	OpenSMOKE::OpenInputFileXML(doc, xml_string,path_kinetics / "kinetics.xml");
    Info <<"-----------2"<<endl;

	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
	thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>(doc);
	//transportPropertiesMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>(doc);
	kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN<double>(*thermodynamicsMapXML, doc);

	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();


	std::cout << " * Time to read XML file: " << tEnd-tStart << std::endl;
}
const unsigned int ns  = thermodynamicsMapXML->NumberOfSpecies();
const unsigned int nr = kineticsMapXML->NumberOfReactions();
