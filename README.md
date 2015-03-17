Galaxy-M
========

Metabolomics Tools for Galaxy

The folders stored here corresond to folders found in a standard Galaxy installation. Included are the tool files (both original code and .xml wrappers) for Metabolomics analysis and the Galaxy config files from a working installation of Galaxy. Not all of the config files have been altered, but all are included here for completeness. 

The Galaxy version that these files have been tested on is from the Stable branch on Mercurial: changeset 14929:5fc83b69fe24, date Thu Dec 11.

Availability and Requirements
=============================

###Programming languages (versions given were used for development, other versions may also be compatible): 
*Python (version 2.7), 
*Matlab (version 2012a), 
*PLS-Toolbox for PCA tools (version 7.0.3) 
*R programming language (version 3.0.1, x86 64bit)

###Other requirements: 
*Galaxy [www.getgalaxy.org - Stable branch on Mercurial: changeset 14929:5fc83b69fe24, date Thu Dec 11], 
*MI-Pack [https://github.com/Viant-Metabolomics/MI-Pack], 
*WineHQ (version 1.6.2, [winehq.org]), 
*XCMS [https://metlin.scripps.edu/xcms/] 
*MSFileReader package (Thermo Scientific [http://sjsupport.thermofinnigan.com/public/detail.asp?id=703]).

###License: 
*GNU General Public License version 3.0 (GPLv3)

###Any restrictions to use by non-academics: 
*no restrictions on use

Installation instructions for Ubuntu 13.10 64bit
================================================

Setting up GalaxyM on Ubuntu

Requirements/System

Ubuntu 13.10 64bit edition (we’re running it on VMware so it can be uploaded to AWS)

Step 0. Update the package manager

You may need to run this with superuser privileges, I’ll assume that you do so my commands will begin ‘sudo’. This requires you to provide your password and Ubuntu will check that you have superuser privileges

If you are running a 64bit (amd) version of Ubuntu (13.10) then it may be problematic to install WINE because it expects an i386 architecture. A workaround is:

$> sudo dpkg —add-architecture i386
$> sudo apt-get update

In any event, you should start with 

$> sudo apt-get update

Step 1. Install WINE

You may wish to see the official instructions here: https://www.winehq.org/download/ubuntu

If using Ubuntu 64bit, there are problems with WINE and Python2.7 sometimes. WINE is notoriously difficult. Some nice people have set up an Ubuntu WINE repository and so it’s advised to set that up first:

$> sudo add-apt-repository ppa:ubuntu-wine/ppa

Then update

$> sudo apt-get update

It is then advised to install the latest version of WINE (made available by the new repository). I’ll be installing WINE1.7 but you may wish to check which versions are available to you using:

$> sudo apt-cache search wine

We then perform the standard ‘apt-get install’ on our chosen WINE version. On running the command, the package manager will ask you if you are sure you want to install the package- type Y to confirm.

$> sudo apt-get install wine1.7

Once it’s finished downloading all the packages, it might present you with an End User Agreement for Microsoft software that you can scroll through with the arrow keys and select the ‘yes’/‘no’ options by pressing Tab and then Enter or Space to select. Having agreed, more installation will proceed. On my most recent WINE1.7 install, none of this happened, but on previous attempts it did.

Step 2. Install Windows packages in WINE

$> mkdir ~/wine-soft
$> cd ~/wine-soft

Download packages: Python 2.7, 32 bit

$> wget https://www.python.org/ftp/python/2.7.8/python-2.7.8.msi

Download packages: numpy 1.8.1

$> wget http://sourceforge.net/projects/numpy/files/NumPy/1.8.1/numpy-1.8.1-win32-superpack-python2.7.exe

Download packages: scipy 0.14.0

$> wget http://sourceforge.net/projects/scipy/files/scipy/0.14.0/scipy-0.14.0-win32-superpack-python2.7.exe

Download packages: comtypes 0.6.1

$> wget http://sourceforge.net/projects/comtypes/files/comtypes/0.6.1/comtypes-0.6.1.win32.exe

Download software: MSFileReader 3.0 sp2

For MSFileReader, you need to go to ThermoFisher’s site and register ( https://thermo.flexnetoperations.com/control/thmo/RegisterMemberToAccount ) then wait for a login to be emailed to you… then you can download the software by following the link in the email - which will take you to a product list, MSFileReader is under Utility Software - I’ve chosen to use version 3.0 sp2

Extract the downloaded .zip file into your ~/wine-soft folder to join the other packages. It will produce a couple of folders, one for standard x86 architecture and one for 64bit. WINE is working off 32bit so the file you are looking for is MSFileReader.exe in the standard (non 64bit) folder. 

Install packages: python2.7
WINE can’t install .MSI files directly. You need to use the msiexec command as follows:

$> wine msiexec /i python-2.7.8.msi

Install packages: numpy

$> wine numpy-1.8.1-win32-superpack-python2.7.exe

Install packages: scipy

$> wine scipy-0.14.0-win32-superpack-python2.7.exe

Install packages: comtypes

$> wine comtypes-0.6.1.win32.exe

Install software: MSFileReader 3.0 sp2
Below I’m printing out the full path for clarity. If your downloaded/extracted files look like mine, this might help. If your file system doesn’t look like this, replace my path directions with wherever you’ve stored the 32bit ‘standalone’ version of MSFileReader.exe

$> wine ~/wine-soft/MSFileReader\ 3.0\ SP2/MSFileReader_x86_Standalone/MSFileReader.exe


Install visual C using winetricks

At this stage, Wine may produce errors because it can’t register the XRawFile2.dll required for MSFileReader. You may already have installed winetricks during the previous step when you installed wine1.7 but in case you haven’t , you can install it using apt-get. Then you need to install Visual C. These two steps can be done as follows (the vc install will ask you to accept various agreements etc, just say yes!):

$> sudo apt-get install winetricks
$> winetricks vcrun2008




WINE-32BIT-ATTEMPT
It may be that MSFileReader requires 32bit (DLLs etc) and that WINE is automatically running in 64bit - being modern etc. In order to get WINE to run in 32bit, remove the .wine folder in your home directory (or move it to e.g. ~/.wine_backup) which will remove all installed software and data from WINE. Then set the ‘WINEARCH’ system variable by updating your ~/.profile with the line “export WINEARCH=win32” (place at end of file, don’t include the quotation marks). Then try installing the software as before…



Step 3: Install MI Pack python package
Download the MI Pack python package (there is a version installed on the GalaxyM VM already, the original package files are at /home/galaxym/Galaxy_MI_Pack/MI_Pack_Python_package/) and install as you would any standard Python package. That is to say:

sudo python [path/to/package/]setup.py build
sudo python [path/to/package/]setup.py install

Then if you have multiple processing cores on your system and wish to make use of parallel processing to speed up e.g. Empirical Formula search, install Parallel Python. This can be downloaded from parallelpython.com or there is a version within the MI Pack folder on the VM. Install this package in the same way as MI Pack (sudo python setup.py install). It is NOT necessary to have parallelpython installed for MI Pack to work, but it is advised if you have multiple cores. 



Step 4: Install Galaxy

You may wish to check out the instructions at Get Galaxy: https://wiki.galaxyproject.org/Admin/GetGalaxy

First, install ‘mercurial’ if it’s not already on your VM. 

$> sudo apt-get install mercurial

Upon download, the Galaxy distribution (using the following settings) will create a folder called ‘galaxy-dist’ in the directory that you call the command from. So if you are copying these instructions verbatim, make sure you’re in a directory you’re happy to work from. I’ll be in my user’s home directory /home/galaxym. 

$> hg clone https://bitbucket.org/galaxy/galaxy-dist/
$> cd galaxy-dist
$> hg update stable

Step 5: Add GalaxyM to Galaxy

Obtain the latest version of GalaxyM from …

1. Merge or replace the 'galaxy-dist/tools/' directory with 'GalaxyM/tools/'
2. Merge or replace the 'galaxy-dist/test-data/' directory with 'GalaxyM/test-data/'
3. Merge or replace the 'galaxy-dist/tool-data/' directory with 'GalaxyM/tool-data/'

4. edit or replace ‘galaxy-dist/config/tool_conf.xml’ with [the contents of] ‘GalaxyM/config/tool_conf.xml’ 

5. edit/replace/create ‘galaxy-dist/config/datatypes_conf.xml’ with the SQLite ‘datatype extension’ tags from ‘GalaxyM/config/datatypes_conf.xml’ (or simply copy the GalaxyM file to your galaxy-dist/config/ directory).

6. edit ‘galaxy-dist/lib/galaxy/datatypes/binary.py’ to include the 5 SQLite ‘classes’ and subsequent ‘Binary.register_sniffable…’ commands from  ‘GalaxyM/lib/galaxy/datatypes/binary.py’ or simply replace the standard file with ours. 

7. If you were already running Galaxy, restart now. 

NB: to run Galaxy, you need to call the run.sh bash script within the main galaxy-dist folder e.g.

$> sh ~/galaxy-dist/run.sh

If run in this way, Galaxy can be stopped by typing Ctrl-C.

Step 6. Install Matlab

a) install Java RunTime Environment and WebStart. 
This enables you to use the Matlab installers as downloaded from the Mathworks website. Alternatively you can just download the product files directly. To install Java JRE, first check which version is available:

$> sudo apt-cache search openjdk-*

Within the results you should see packages with names like ‘openjdk-6-jre’. The number indicates the Java version. I’ll be installing the latest one available in my package list, version 7.

$> sudo apt-get install openjdk-7-jre

In order to use the Mathworks’ ‘Download Agent’, you’ll also need java webstart. On Ubuntu, you can download this from the apt-get repository as follows:

$> sudo apt-get install icedtea-netx

b) Download Matlab and Statistics Toolbox (version 2012a) using Download Agent
Within your Mathworks website account, navigate to your licenses and start the download products process. The site intends you to use their Java ‘webstart’ based Download Agent but there are alternatives. The website will try to see that you have Java installed before download. If it can’t find your Java JRE but you installed it in the previous step, just click ‘i have java’ and continue. A file called Download Agent will start to download. 

Your internet browser may ask if you want to open with IcedTea - feel free to do so. Alternatively download the file and then open with IcedTea. In theory you ought to be able to do this via the commandline using ‘javaws /path/to/my/file’ but I found this produced an error. However, simply opening up a file browser and double clicking on the file automatically opened it with the Iced Tea tool. 

Simply follow all defaults to install Matlab. The only exception will be that Matlab will try to create a folder called /usr/local/MATLAB/R2012a but you may not have started the process with superuser privileges and so the installer won’t be allowed to do this. If there is an error, simply use the following commands (at the commandline) and then tell it to try again:

$> sudo mkdir -p /usr/local/MATLAB/R2012a
$> sudo chmod 777 /usr/local/MATLAB/R2012a

Of course, you should make sure you’ve created the same path that the installer wanted to make or that you change the path when running the installer. Assuming that this allows you to progress, and you can activate your software etc. it is a good idea to go back and change the permissions on that folder, undoing the ‘chmod 777’ that we just did to allow anyone (including the installer) to write to that location. 

$> sudo chmod 755 /usr/local/MATLAB/R2012a

Now, for ease, you should put a symbolic link between the matlab binary that you’ll want to run and the /usr/local/bin directory that holds most of these programs. This is optional but it means you can simply type ‘matlab’ at the command line and Linux will know what you mean.

$> sudo ln -s /usr/local/MATLAB/R2012a/bin/matlab /usr/local/bin/matlab

NB: for some versions of 64 bit Ubuntu (13.10 for example), Matlab declares an error when it runs - ‘blah blah blah… /lib64/libc.so.6: not found’. Technically it will run and seems to perform the GalaxyM tasks but this warning is seen as an error by Galaxy and so despite storing the correct output in the history, it appears in red (i.e. failed job) and that’s not very helpful. A workaround is to find the libc.so.6 file and place a symbolic link in the location that Matlab is looking. 

First find the libc.so.6 file:

$> locate libc.so.6

In my case this produces the result: /lib/x86_64-linux-gnu/libc.so.6
Next step is to create the symbolic link between that result (where the file actually is) and the directory that Matlab is looking for it (/lib64/libc.so.6)

$> sudo ln -s /lib/x86_64-linux-gnu/libc.so.6 /lib64/libc.so



Step 7. Install PLS Toolbox

Log into the Eigenvector.com website and navigate to software downloads. Choose your version (we’re using 7.0.3) and download the .zip file. Assuming you’ve downloaded that to the ‘/home/your_user/Downloads’ folder (and that it’s version 7.0.3 and that you put your Matlab install in the default location), you can unzip it and place it in the Matlab toolbox folder as follows:

$> unzip ~/Downloads/PLS_Toolbox_703.zip -d /usr/local/MATLAB/R2012a/toolbox/

We now need to run the install program which is a matlab script. Open up Matlab from the commandline:

$> matlab

Once inside the Matlab environment, you’ll have to specify the location of your PLS_Toolbox folder (i.e. /usr/local/MATLAB/R2012a/toolbox/PLS_Toolbox_703/ ) so Matlab knows where the PLS scripts are. Using the GUI, this can be done from drop down menus: File-> Set Path -> Add with subfolders. Or it can be done from the Matlab command line as follows:

>> addpath(genpath(‘/usr/local/MATLAB/R2012a/toolbox/PLS_Toolbox_703/‘))

Once PLS_Toolbox is added to the path, return to the matlab commandline and type ‘evriinstall’ to run the installation wizard. 

>> evriinstall

Here you will be asked for a license code - obtainable from the download page on the Eigenvector website. Accepting all the defaults should be fine although if you’re very keen to have exactly version 7.0.3 to copy the original GalaxyM installation, you may wish to uncheck the box (underneath where you enter the license code) that offers to look for updates and newer versions of code etc.


Step 8. Install R
The LC-MS pipeline makes use of XCMS for an initial peak picking and alignment processing. This requires ‘R’ and various packages. R can be installed using:

$> sudo apt-get install r-base
$> sudo apt-get install r-cran-rcpp
$> sudo apt-get install libnetcdf-dev

Then you need to run R and install Hmisc,XCMS and CAMERA 

From the commandline, start R by typing simply:

$> sudo R

NB: it is advisable to run with sudo because various packages complain about not being able to write to various directories. Then within R, use ‘install.packages’ to install Hmisc. 

> install.packages(“Hmisc”, dependencies=TRUE )

XCMS and CAMERA are installed from Bioconductor so the command is slightly different

> source("http://bioconductor.org/biocLite.R")
> biocLite(“BiocUpgrade”)
> biocLite("xcms")
> biocLite("CAMERA")


Step 9. Get the data

The test data sets (both LCMS and DIMS) can be obtained from their larger data sets in Metabolights or from GigaDB. Download and if you want to copy the GalaxyM paper, move to ~/GalaxyM-TestData (with subfolders DIMS_DATA and LCMS_DATA for each modality). 

Step 10. Run Galaxy

from the Linux commandline, start Galaxy if it’s not already running by typing:

$> sh ~/galaxy-dist/run.sh

To interact with Galaxy, open up a web-browser and point it at your server. If you have access to a browser on the same system as Galaxy, you can load that up and enter 127.0.0.1:8080 in the address bar. 127.0.0.1 means ‘localhost’ and the browser speaks to the system it is running on. Alternatively, if you are running the server on a remote system (in the cloud perhaps) then you’ll need to ensure that the ‘galaxy.ini’ file (galaxy-dist/config/galaxy.ini) has the line ‘host = 0.0.0.0’ rather than the default ‘host = 127.0.0.1’. This line tells Galaxy to expect traffic from the internet rather than just local requests. That setting has already been changed in the GalaxyM file (you may wish to change it back to 127.0.0.1 if you want to protect your Galaxy server from the wide web). 

Hopefully, you should be presented with  the main Galaxy landing page. There should be a list of tools in the left hand panel, a welcome page in the middle and a history in the right hand side. 

Step 11. Log into Galaxy if you want to make use of workflows

The GalaxyM installation has a user already registered and two workflows have been stored. If you wish to access these then the login details are:

user: galaxym@galaxm.org
password: galaxym


Location of data in published work:
When logged in as galaxym user on the Ubuntu system,
~/GalaxyM-TestData/LCMS_DATA
~/GalaxyM-TestData/DIMS_DATA


 


