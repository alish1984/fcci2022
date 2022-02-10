# install gcc 8.4.0 
sudo apt-get install gcc-8 g++-8

# gcc --version 
# it may give you sth like the below output (first line)
# gcc (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0
# which shows you have gcc 7 
# now we want to add an alternative just compiled, i.e. gcc 8.4.0 
# depending on your default version, you may need to change the version (e.g. gcc-8  gcc++-8) in the first command 

sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-8 1 --slave /usr/bin/g++ g++ /usr/bin/g++-8
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 2 --slave /usr/bin/g++ g++ /usr/bin/g++-7


#now that you have introduced your options, choose gcc-5 by running below command and choosing the right option

sudo update-alternatives --config gcc

#here you can choose your default version, after pressing Enter, the default gcc of your ubuntu will be modified
# for example select 2
#check you did correctly

gcc --version 


