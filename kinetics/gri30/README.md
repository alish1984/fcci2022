# to convert to opensmoke format 
- ofv2106 
- openSMOKEppCHEMKINPreProcessor --input input.dic 

# to convert to cantera (cti) format 
- conda activate ct-fcci2022
- ck2cti --input=chemkin/grimech30.dat --thermo=chemkin/thermo30.dat --transport=chemkin/transport.dat --output=cti/gri30.cti --permissive
