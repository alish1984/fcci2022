# A repository for fcci2022 workshop on pre/post-processing with cantera and openSMOKE++

![image info](./.cantera-logo_s.png)
![image info](./.logo_helvetica_s.png)

## 0. Prerequisites 
- [x] Basic knowledge of Linux shell commands
- [x] Basic knowledge of python language
- [x] Basic knowledge of C++ language
- [x] Basic knowledge of OpenFOAM, case setup and coding
- [x] Advanced (MSc.) knowledge of theory of combustion 

## 1. Installation dependencies
| :warning: We use Linux! Specifically all below commands have been tested on ubuntu 18.04!✨|
|--|

- ### _For cantera (pre/post-processing)_
    * You need python (>=3.), anaconda (>=v3.), cantera (>=2.5)

- ### _For openSMOKE++ source (post-processing)_
    * You need boost, eigen and rapid xml libraries
    * eigen and rapid xml are available in (externalLibs)
    * openSMOKE++ source is available in (externalLibs)
    
	* If you use openSMOKE for your publications, we kindly ask you to cite the following two papers:
    	> Cuoci, A., Frassoldati, A., Faravelli, T., Ranzi, E., 
    	> OpenSMOKE++: An object-oriented framework for the numerical modeling of reactive systems with detailed kinetic mechanisms 
    	> (2015) Computer Physics Communications, 192, pp. 237-264, DOI: 10.1016/j.cpc.2015.02.01
	
	-----------

## 2. Installation guides
- ### 2.1  Define repo paths  
    ```sh
    echo "export fcci2022=${HOME}/projects/fcci2022" >> ${HOME}/.bashrc
    echo "export fcci2022_postostools=${HOME}/projects/fcci2022/postProcessing/openSMOKE/tools" >> ${HOME}/.bashrc
    echo "export fcci2022_postostuts=${HOME}/projects/fcci2022/postProcessing/openSMOKE/tutorials" >> ${HOME}/.bashrc
    echo "export fcci2022_postoscases=${HOME}/projects/fcci2022/postProcessing/openSMOKE/cases" >> ${HOME}/.bashrc
    ```
- ### 2.2 Install anaconda, jupyter lab, 
    - Anaconda 3 
    * https://docs.anaconda.com/anaconda/install/
        ```sh
        sudo apt-get install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6
        ```
    * Download (save) anaconda installer, you can use the latest version or choose from Archive, I have already tested the 2020.02 version which I got from Archive in 
    * https://repo.anaconda.com/archive/ 
	* https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
    * Then if you have downloaded it in ~/Downloads 
        ```sh
        bash ~/Downloads/Anaconda3-2020.02-Linux-x86_64.sh
        ```
    * During installation I usually set the location (<YOUR_HOME_PATH>/NumericalLibraries)
    * The installer prompts “Do you wish the installer to initialize Anaconda3 by running conda init?” We recommend “yes”.
        ```sh
        conda config --set auto_activate_base False
        ```
    * see installationScripts/conda_cantera_python.md
        ```sh
        cd $fcci2022/installationScripts
        conda env create -f conda_cantera_python.yaml
        python -m ipykernel install --user --name ct-fcci2021-2
        ```
    * to open Jpyterlab 
        ```sh
        cd $fcci2022/preProcessing/cantera
        jupyter lab
        ```

- ### 2.3 Install gcc 8.4 
	* see installationScripts/compile-README-gcc8.4.0.txt
- ### 2.4 Install boost 
	* see installationScripts/compileBoost.sh
- ### 2.5 Install OpenFOAMv2106 
    * see installationScripts/preReqCompileOpenFOAM.md
	* see installationScripts/compileOpenFOAM.sh



