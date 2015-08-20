Galaxy-M
========

Metabolomics Tools for [Galaxy](http://galaxyproject.org/)

The folders stored here corresond to folders found in a standard Galaxy installation. Included are the tool files (both original code and .xml wrappers) for Metabolomics analysis and the Galaxy config files from a working installation of Galaxy. Not all of the config files have been altered, but all are included here for completeness. 

The Galaxy version that these files have been tested on is from the Master branch on Github: commit c429777c93680dcee449fe410f5360afbe673758.

The MI-Pack version that the GigaScience publication utilises is Github commit: 06ce0ace643ee1cf1d27550769a5272b7ea50825

Availability and Requirements
=============================

###Programming languages (versions given were used for development, other versions may also be compatible): 
* [Python (version 2.7)](https://www.python.org/download/releases/2.7/)
* [R programming language (version 3.0.1, x86 64bit)](http://cran.r-project.org/bin/windows/base/)
* [MATLAB Compiler Runtime (MCR) (version 8.1)](http://uk.mathworks.com/supportfiles/MCR_Runtime/R2013a/MCR_R2013a_glnxa64_installer.zip) - (Note: The MATLAB Runtime is a standalone set of shared libraries that enables the execution of compiled MATLAB applications or components on computers that do not have MATLAB installed. Each [release](https://github.com/Viant-Metabolomics/Galaxy-M/releases) includes executables for most of the tools)
* [Matlab (version 2012a)](http://uk.mathworks.com/products/matlab/) (Optional: tool development and use of [most recent source code](https://github.com/Viant-Metabolomics/Galaxy-M))
* [PLS-Toolbox for PCA tools (version 7.0.3)](http://www.eigenvector.com/software/pls_toolbox.htm) (Optional: few of the tools are PLS-toolbox dependent)

###Other requirements: 
* [Galaxy (Master branch on Github: commit c429777c93680dcee449fe410f5360afbe673758)](http://www.getgalaxy.org)
* [MI-Pack (Metabolite Identification Package, Master branch on Github: commit 06ce0ace643ee1cf1d27550769a5272b7ea50825)](https://github.com/Viant-Metabolomics/MI-Pack)
* [XCMS](https://metlin.scripps.edu/xcms/)
* [Thermo Scientific MSFileReader package](http://sjsupport.thermofinnigan.com/public/detail.asp?id=703)
* [WineHQ (version 1.6.2)](winehq.org)

###License: 
* [GNU General Public License version 3.0 (GPLv3)](https://www.gnu.org/licenses/quick-guide-gplv3.html)

###Any restrictions to use by non-academics: 
* no restrictions on use

###Pre-Installed Virtual Machine 
A virtual machine implementation has been created to support reuse by the community as well as to act as a snapshot of reproducibility for the work in the GigaScience publication. This will be available via the GigaDB repository (accession details to follow). You can use the VM in VMWare or VirtualBox, there is one virtual hard disk (VMDK file) that suits both systems however we've had to create separate metadata files for each (MF and OVF files). To import into either system you need the VMDK to be in the same folder as the relevant MF and OVF file. 

Installation instructions to recreate the virtual machine can be found below. They are quite complex and can involve a lot of troubleshooting, so it is advised to work direct from the virtual machine where possible. We will aim to share the VM via Amazon AWS for access without download. The GalaxyM installation has a user already registered and two workflows have been stored/published. If you wish to access the original user account then the login details are:

Ubuntu logins:

    user: galaxym
    password: galaxym
    
Galaxy logins:

    email: galaxym@galaxym.org
    user: galaxym
    password: galaxym


Location of data in published work:

When logged in as ``galaxym`` user on the Ubuntu system, the Mass Spectra can be found in :

* `~/GalaxyM-TestData/LCMS_DATA`
* `~/GalaxyM-TestData/DIMS_DATA`

The pipelines/workflows also make use of CSV files for certain instructions and these can be found in the GalaxyM-TestData folder mentioned above. See the published histories or workflows for where to use these files. 


Installation instructions for Ubuntu 14.04LTS 64bit 
================================================


###Step 0. Update the package manager

You may need to run this with superuser privileges, I’ll assume that you do so my commands will begin ‘sudo’. This requires you to provide your password and Ubuntu will check that you have superuser privileges

If you are running a 64bit (amd) version of Ubuntu (14.04) then it may be problematic to install WINE because it expects an i386 architecture. A workaround is:

    $> sudo dpkg --add-architecture i386
    $> sudo apt-get update

In any event, you should start with 

    $> sudo apt-get update

###Step 1. Install WINE

You may wish to see the official instructions [here](https://www.winehq.org/download/ubuntu).

If using Ubuntu 64bit, there are problems with WINE and Python2.7 *sometimes*. WINE is notoriously difficult. Some nice people have set up an Ubuntu WINE repository and so it’s advised to set that up first: 

To set up the WINE/Ubuntu repository:

    $> sudo add-apt-repository ppa:ubuntu-wine/ppa

Then update

    $> sudo apt-get update

It is then advised to install the latest version of WINE (made available by the new repository). I’ll be installing WINE1.7 but you may wish to check which versions are available to you using:

    $> sudo apt-cache search wine
    
We then perform the standard ‘apt-get install’ on our chosen WINE version. On running the command, the package manager will ask you if you are sure you want to install the package - type Y to confirm.

    $> sudo apt-get install wine1.7

Once it’s finished downloading all the packages, it might present you with an End User Agreement for Microsoft software that you can scroll through with the arrow keys and select the ‘yes’/‘no’ options by pressing Tab and then Enter or Space to select. Having agreed, more installation will proceed. On my most recent WINE1.7 install, none of this happened, but on previous attempts it did.

IMPERATIVE!: You must force wine to install everything in 32bit mode (because Thermo Fisher won't update their dlls). In order to get WINE to run in 32bit, remove the .wine folder in your home directory (or move it to e.g. ~/.wine_backup) which will remove all installed software and data from WINE. Then set the ``WINEARCH`` system variable by updating your ``~/.profile`` file with the line ``export WINEARCH=win32`` (place at end of file). Then try installing the software as follows... if you do not do this, you will likely get errors from the SimStitch tools along the lines of 'Invalid file identifier' - this is because they can't use the MSFileReader tool to open the DIMS data. 
Note that if this is the first time wine has been installed and it has yet to be configured the .wine folder will not have been created yet so will not need to be removed.

###Step 2. Install Windows packages in WINE

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
WINE can’t install .MSI files directly. You need to use the msiexec command (see below). When presented with Python installation GUI, just accept the defaults. Update: For Ubuntu 14.04, you may be asked to confirm installation of the Wine MONO installer, and Wine Gecko installer etc. Just agree to 'install':

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


Install visual C using winetricks (NB: This will not work in China. download.microsoft.com is not accessible, tested 7th Aug 2015)

At this stage, Wine may produce errors because it can not register the XRawFile2.dll required for MSFileReader. 
You may already have installed winetricks during the previous step when you installed wine 1.7 but in case you haven not, you can install it using apt-get. Then you need to install Visual C. 
These two steps can be done as follows (the vc install will ask you to accept various agreements etc, just say yes!):

    $> sudo apt-get install winetricks
    $> winetricks vcrun2008


###Step 3: Install MI-Pack python package
Download the MI-Pack python package from https://github.com/Viant-Metabolomics/MI-Pack (there is a version installed on the GalaxyM VM already, the original package files are at ``/home/galaxym/Galaxy_MI_Pack/MI_Pack_Python_package/``) and install as you would any standard Python package. That is to say:
    
    $> sudo apt-get install git
    $> cd ~
    $> git clone https://github.com/Viant-Metabolomics/MI-Pack.git
    $> cd MI-Pack
    $> sudo python setup.py build
    $> sudo python setup.py install

Then if you have multiple processing cores on your system and wish to make use of parallel processing to speed up e.g. Empirical Formula search, install Parallel Python. This can be downloaded from [parallelpython.com](http://www.parallelpython.com/) or there is a version within the MI-Pack folder on the VM.
Install this package in the same way as MI-Pack (``sudo python setup.py install``). 
It is NOT necessary to have parallel python installed for MI-Pack to work, but it is advised if you have multiple cores. 

    $> cd ~
    $> curl --remote-name http://www.parallelpython.com/downloads/pp/pp-1.6.4.tar.gz
    $> tar -xzvf pp-1.6.4.tar.gz
    $> cd pp-1.6.4/
    $> sudo python setup.py build
    $> sudo python setup.py install


###Step 4: Install Galaxy

You may wish to check out the instructions at [Get Galaxy](https://wiki.galaxyproject.org/Admin/GetGalaxy).

Upon download, the Galaxy distribution (using the following settings) will create a folder called ``galaxy`` in the directory that you call the command from. So if you are copying these instructions verbatim, make sure you're in a directory you’re happy to work from. I’ll be in my user’s home directory ``/home/galaxym``. 

    $> cd ~
    $> git clone https://github.com/galaxyproject/galaxy/
    $> cd galaxy
    $> git checkout -b master origin/master

###Step 5: Add GalaxyM to Galaxy

Obtain the latest version of GalaxyM from [github](https://github.com/Viant-Metabolomics/Galaxy-M).

    $> cd ~
    $> git clone https://github.com/Viant-Metabolomics/Galaxy-M.git

Merge or replace the:

``galaxy/tools/`` directory with ``GalaxyM/tools/``

    $> cp -r Galaxy-M/tools/ galaxy/

``galaxy/test-data/`` directory with ``GalaxyM/test-data/``

    $> cp -r Galaxy-M/tool-data/ galaxy/

``galaxy/tool-data/`` directory with ``GalaxyM/tool-data/``

    $> cp -r Galaxy-M/test-data/ galaxy/

Edit/replace/create ``galaxy/config/tool_conf.xml`` with [the contents of] ``GalaxyM/config/tool_conf.xml`` 

    $> cp Galaxy-M/config/tool_conf.xml galaxy/config/

Edit/replace/create ``galaxy/config/datatypes_conf.xml`` with the SQLite ``datatype extension`` tags from ``GalaxyM/config/datatypes_conf.xml`` (or simply copy the GalaxyM file to your galaxy/config/ directory).

    $> cp Galaxy-M/config/datatypes_conf.xml galaxy/config/

Copy the Galaxy-M binary type definitions file from ``GalaxyM/lib/galaxy/datatypes/galaxym.py`` into the ``galaxy/lib/galaxy/datatypes/`` directory.

    $> cp Galaxy-M/lib/galaxy/datatypes/galaxym.py galaxy/lib/galaxy/datatypes/

Copy the Galaxy-M welcome page files from ``GalaxyM/static/`` into the ``galaxy/static/`` directory.

    $> cp -r Galaxy-M/static/welcome.html galaxy/static/
    $> cp -r Galaxy-M/static/images/logos/ galaxy/static/images/

If you were already running Galaxy, restart now.

NB: to run Galaxy, you need to call the run.sh bash script within the main ``galaxy`` folder e.g.

    $> sh ~/galaxy/run.sh

If run in this way, Galaxy can be stopped with the command Ctrl-C.


###Step 6. MATLAB Compiler Runtime (MCR)

The MATLAB Runtime is a standalone set of shared libraries that enables the execution of compiled MATLAB applications or components on computers that do not have MATLAB installed. Each [release](https://github.com/Viant-Metabolomics/Galaxy-M/releases) includes executables for most of the tools

    $> cd ~
    $> curl --remote-name http://uk.mathworks.com/supportfiles/MCR_Runtime/R2013a/MCR_R2013a_glnxa64_installer.zip
    $> unzip MCR_R2013a_glnxa64_installer.zip -d MCR_R2013a_glnxa64_installer
    $> cd MCR_R2013a_glnxa64_installer
    $> sudo ./install
    
Follow the instructions on the screen to install MCR. With the following EXCEPTION:
The Matlab installation will tell you to add certain environment variables to your path. Do NOT do this - it will prevent Galaxy from starting and also your Ubuntu Unity desktop system from logging in/working. 

However, you must put something in so, on the basis that you've copied the instructions above, installed MCR R2013a, which is Version 81, edit your ~/.profile file to include the following:

    export LD_LIBRARY_PATH=/usr/local/MATLAB/MATLAB_Compiler_Runtime/v83/runtime/glnxa64:/usr/local/MATLAB/MATLAB_Compiler_Runtime/v83/bin/glnxa64:/usr/local/MATLAB/MATLAB_Compiler_Runtime/v83/sys/os/glnxa64:

Set the XAPPLRESDIR environment variable to the following value

    export XAPPLRESDIR=/usr/local/MATLAB/MATLAB_Compiler_Runtime/v83/X11/app-defaults

What we've done here is deliberately give the wrong version number (v83 instead of v81). Oddly, this seems to allow everything to work! 

Also: Point the MCR_CACHE_ROOT environment variable to a local temporary directory (failure to do so produces the error: "Could not access the MCR component cache.") 

    MCR_CACHE_ROOT=/tmp/mcr_cache_$USER/

Having added all these environment variables to your ~/.profile file, you need to refresh your paths by typing:

    $> source ~/.profile
    
Then manually create that temporary MCR_CACHE_ROOT directory:

    $> mkdir $MCR_CACHE_ROOT
    
For more info on Environment Variables see [here](https://help.ubuntu.com/community/EnvironmentVariables)  (Section Persistent environment variables).


###Step 7. Install R

The LC-MS pipeline makes use of XCMS for an initial peak picking and alignment processing. This requires ‘R’ and various packages. R can be installed using:

    $> sudo apt-get install r-base
    $> sudo apt-get install r-cran-rcpp
    $> sudo apt-get install libnetcdf-dev netcdf-bin

Then you need to run R and install [Hmisc](http://cran.r-project.org/web/packages/Hmisc/index.html), [XCMS](https://metlin.scripps.edu/xcms/) and [CAMERA](http://bioconductor.org/packages/release/bioc/html/CAMERA.html) 

From the commandline, start R by typing simply:

    $> sudo R

NB: it is advisable to run with sudo because various packages complain about not being able to write to various directories. Then within R, use ``install.packages`` to install Hmisc. 

    > install.packages('Hmisc', dependencies=TRUE )
    
If errors appear, saying that e.g. "Installation of packages 'ggplot2' had non-zero exit status", check your R version and if it's less than 3.1 (Ubuntu 14.04LTS apt-get base R is 3.0.2), then you'll need to upgrade. We followed these instructions for the Virtual Machine released to support this package (http://sysads.co.uk/2014/06/install-r-base-3-1-0-ubuntu-14-04/)

XCMS and CAMERA are installed from Bioconductor so the command is slightly different

    > source("http://bioconductor.org/biocLite.R")
    > biocLite()
    > biocLite("xcms")
    > biocLite("CAMERA")

###Step 8. Get the data

The test data sets (both LCMS and DIMS) can be obtained from their larger data sets in Metabolights (DIMS: Accession MTBLS79; LCMS: Accession MTBLS146) or from GigaDB (accession to follow). In the case of the DIMS dataset, the MetaboLights accession does not include the necessary .dat files that allow SimStitch processing, therefore the GigaDB repository is required (or contact Viant Lab). Download the datasets and, if you want to copy the GalaxyM paper, move to ``~/GalaxyM-TestData`` (with subfolders DIMS_DATA and LCMS_DATA for each modality). 

###Step 9. Run Galaxy

from the Linux commandline, start Galaxy if it’s not already running by typing:

    $> sh ~/galaxy/run.sh

To interact with Galaxy, open up a web-browser and point it at your server. If you have access to a browser on the same system as Galaxy, you can load that up and enter ``127.0.0.1:8080`` in the address bar. ``127.0.0.1`` means ``localhost`` and the browser speaks to the system it is running on. Alternatively, if you are running the server on a remote system (in the cloud perhaps) then you’ll need to ensure that the ``galaxy.ini`` file (``galaxy/config/galaxy.ini``) has the line ``host = 0.0.0.0`` rather than the default ``host = 127.0.0.1``. This line tells Galaxy to expect traffic from the internet rather than just local requests. That setting has already been changed in the GalaxyM file (you may wish to change it back to ``127.0.0.1`` if you want to protect your Galaxy server from the wide web). 

Hopefully, you should be presented with  the main Galaxy landing page. There should be a list of tools in the left hand panel, a welcome page in the middle and a history in the right hand side. 

###Install Matlab (Optional: tool development and use of most recent source code)

7a) install Java RunTime Environment and WebStart. 
This enables you to use the Matlab installers as downloaded from the Mathworks website. Alternatively you can just download the product files directly. To install Java JRE, first check which version is available:

    $> sudo apt-cache search openjdk-*

Within the results you should see packages with names like ``openjdk-6-jre``. The number indicates the Java version. I’ll be installing the latest one available in my package list, version 7.

    $> sudo apt-get install openjdk-7-jre

In order to use the Mathworks’ ‘Download Agent’, you’ll also need java webstart. On Ubuntu, you can download this from the apt-get repository as follows:

    $> sudo apt-get install icedtea-netx

7b) Download Matlab and Statistics Toolbox (version 2012a) using Download Agent
Within your Mathworks website account, navigate to your licenses and start the download products process. The site intends you to use their Java ‘webstart’ based Download Agent but there are alternatives. The website will try to see that you have Java installed before download. If it can’t find your Java JRE but you installed it in the previous step, just click ‘i have java’ and continue. A file called Download Agent will start to download. 

Your internet browser may ask if you want to open with IcedTea - feel free to do so. Alternatively download the file and then open with IcedTea. In theory you ought to be able to do this via the commandline using ‘javaws /path/to/my/file’ but I found this produced an error. However, simply opening up a file browser and double clicking on the file automatically opened it with the Iced Tea tool. 

Simply follow all defaults to install Matlab. The only exception will be that Matlab will try to create a folder called /usr/local/MATLAB/R2012a but you may not have started the process with superuser privileges and so the installer won’t be allowed to do this. If there is an error, simply use the following commands (at the commandline) and then tell it to try again:

    $> sudo mkdir -p /usr/local/MATLAB/R2012a
    $> sudo chmod 777 /usr/local/MATLAB/R2012a

Of course, you should make sure you’ve created the same path that the installer wanted to make or that you change the path when running the installer. Assuming that this allows you to progress, and you can activate your software etc. it is a good idea to go back and change the permissions on that folder, undoing the ``chmod 777`` that we just did to allow anyone (including the installer) to write to that location. 

    $> sudo chmod 755 /usr/local/MATLAB/R2012a

Now, for ease, you should put a symbolic link between the matlab binary that you’ll want to run and the ``/usr/local/bin`` directory that holds most of these programs. This is optional but it means you can simply type ``matlab`` at the command line and Linux will know what you mean.

    $> sudo ln -s /usr/local/MATLAB/R2012a/bin/matlab /usr/local/bin/matlab

NB: for some versions of 64 bit Ubuntu (13.10 for example), Matlab declares an error when it runs - ‘blah blah blah… /lib64/libc.so.6: not found’. Technically it will run and seems to perform the GalaxyM tasks but this warning is seen as an error by Galaxy and so despite storing the correct output in the history, it appears in red (i.e. failed job) and that’s not very helpful. A workaround is to find the ``libc.so.6`` file and place a symbolic link in the location that Matlab is looking. 

First find the ``libc.so.6`` file:

    $> locate libc.so.6

In my case this produces the result: ``/lib/x86_64-linux-gnu/libc.so.6``
Next step is to create the symbolic link between that result (where the file actually is) and the directory that Matlab is looking for it (``/lib64/libc.so.6``)

    $> sudo ln -s /lib/x86_64-linux-gnu/libc.so.6 /lib64/libc.so

###Install PLS Toolbox (Optional: tool development and use of most recent source code)

Log into the Eigenvector.com website and navigate to software downloads. Choose your version (we’re using 7.0.3) and download the .zip file. Assuming you’ve downloaded that to the ``/home/your_user/Downloads`` folder (and that it’s version 7.0.3 and that you put your Matlab install in the default location), you can unzip it and place it in the Matlab toolbox folder as follows:

    $> unzip ~/Downloads/PLS_Toolbox_703.zip -d /usr/local/MATLAB/R2012a/toolbox/

We now need to run the install program which is a matlab script. Open up Matlab from the commandline:

    $> matlab

Once inside the Matlab environment, you’ll have to specify the location of your PLS_Toolbox folder (i.e. ``/usr/local/MATLAB/R2012a/toolbox/PLS_Toolbox_703/``) so Matlab knows where the PLS scripts are. Using the GUI, this can be done from drop down menus: File-> Set Path -> Add with subfolders. Or it can be done from the Matlab command line as follows:

    >> addpath(genpath(‘/usr/local/MATLAB/R2012a/toolbox/PLS_Toolbox_703/‘))

Once PLS_Toolbox is added to the path, return to the matlab commandline and type ``evriinstall`` to run the installation wizard. 

    >> evriinstall

Here you will be asked for a license code - obtainable from the download page on the Eigenvector website. Accepting all the defaults should be fine although if you’re very keen to have exactly version 7.0.3 to copy the original GalaxyM installation, you may wish to uncheck the box (underneath where you enter the license code) that offers to look for updates and newer versions of code etc.


