# 2.1 Install anaconda, jupyter lab, 
- Anaconda 3 
- https://docs.anaconda.com/anaconda/install/linux/
	* sudo apt-get install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6
    * Download (save) anaconda installer, you can use the latest version or choose from Archive, I have already tested the 2020.02 version which I got from Archive in 
    * https://repo.anaconda.com/archive/ 
	* https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
    * Then if you have downloaded it in ~/Downloads 
	* bash ~/Downloads/Anaconda3-2020.02-Linux-x86_64.sh
    * During installation I usually set the location (<YOUR_HOME_PATH>/NumericalLibraries)
    * The installer prompts “Do you wish the installer to initialize Anaconda3 by running conda init?” We recommend “yes”.
	* conda config --set auto_activate_base False

## to create the environment with all dependencies 
conda env create -f conda_cantera_python.yaml

## to check the list of environments 
conda env list



conda activate ct-fcci2021-2
## to check the list of installed libraries 
conda list

## create ipython kernel 
python -m ipykernel install --user --name ct-fcci2021-2

## launch jupyter lab 
jupyter lab

## to remove kernels 
jupyter kernelspec list 
jupyter kernelspec uninstall ct-fcci2021-2


## to remove the environment (not recommended!) 
conda remove --name ct-fcci2021-2 --all





https://docs.anaconda.com/anaconda/install/uninstall/
