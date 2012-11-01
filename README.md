# Evident

E-vident is a web based application developed at the [Knight Lab][] that allows direct effect size comparison and elucidates sampling effort for microbial analysis studies.

A critical consideration in any clinical study is power analysis, yet it has been difficult to perform such analyses for microbiome studies because the effect sizes of different disorders are unknown. Fortunately, several larger cohort studies, including but not limited to the Human Microbiome Project, now allow us to identify effect sizes for differences among ages and populations, and differences associated with obesity, IBD, and other disorders.

--------------------------------------------------------------------------------

# Installation

E-vident is written in [Python][], using [QIIME][] as a base framework and the [Apache][] HTTP Server. The following installation instructions assume both of these dependencies to be correctly installed.

Clone the E-vident repository from the github site:

    git clone git@github.com/qiime/evident.git ~/git_sw/evident

To access all th sample data you will need to create a symbolic link from the repository you just downloaded, like this:

    sudo ln -s /Users/yoshiki/git_sw/evident/data /evident

The following steps will guide you through the installation and configuration of mod_python.

Modpython is an Apache module that embeds the Python interpreter within the server and this is what powers E-vident. To get the package you will need to download the most recent SVN version by running the following command:

    svn co https://svn.apache.org/repos/asf/quetzalcoatl/mod_python/trunk modpython

The installation of mod_python varies depending on the OS version you are using.

## Mac OS X 10.6.x Snow Leopard

Before proceeding, you will need to have the Apple Developer tools installed, this tools can be found in the installation CDs that ship with you computer. If you have Xcode installed in your computer it is very likely that you have the tools already.

    cd modpython
    ./configure --with-apxs=/usr/sbin/apxs
    make

    make clean
    make 
    sudo make install

## Mac OS X 10.7.x Lion and Mountain Lion

To proceed with this installation you will need to have installed [Xcode][] (which you can download from the Mac App Store) after you have done this, the next step is to install the [command line tools for OS X][]. Once this is done you can proceed to build the package.

    cd modpython
    ./configure --with-apxs=/usr/sbin/apxs
    make

## Ubuntu

If you don't have Apache installed already proceed to install it and then install mod_python:

    sudo apt-get install apache2
    sudo apt-get install libapache2-mod-python libapache2-mod-python-doc

# Configuration

The last few steps require you to modify your Apache configuration files located in:

- Mac OS X `/etc/apache2/httpd.conf`
- Ubuntu `/etc/apache2/apache.conf`

You will need to add all of the following lines:

    LoadModule python_module libexec/apache2/mod_python.so
    AddHandler mod_python .psp
    PythonHandler mod_python.psp
    PythonDebug On
    PYTHONPATH "['/Users/yoshiki/Applications/sw/pynast','/Users/yoshiki/Applications/sw/qiime/','/Users/yoshiki/Applications/sw/pycogent/','/Users/yoshiki/Applications/sw/numpy/','/Users/yoshiki/Applications/sw/biom-format','/Users/yoshiki/Applications/sw/matplotlib/']+sys.path" 

Notice that on the last line each path corresponds to a python module that you should have installed in your computer, as an example here all our python modules are living under `/Users/yoshiki/Applications/sw/` hence we have to add the location of each module.

It is possible however that if you have all this modules installed in the site-packages directory of your python installation folder, if this is the case you won't have individual folders per module and you will not need to add all these paths as python should automatically be able to see the information.

Now edit `/private/etc/apache2/users/[username].conf` to add the following lines:

    <Directory "/Users/[username]/Sites/">
        Options Indexes MultiViews FollowSymLinks
        AllowOverride None
        Order allow,deny
        Allow from all
    </Directory>

Restart Apache and you should be all set.

### Installation FAQs & Issues

A known issue when building this package is the following error message:

    llvm-gcc-4.2: error trying to exec '/Developer/usr/bin/../llvm-gcc-4.2/bin/powerpc-apple-darwin10-llvm-gcc-4.2': execvp: No such file or directory
    lipo: can't figure out the architecture type of:

To overcome this issue you will need to install [Xcode 3.2.1][]. Once this is installed, go to the mod-python direcotry and try to re-build the package:

--------------------------------------------------------------------------------

`error: command gcc-4.2 failed with exit status 1`

Locate GCC by running `which gcc` from the terminal and create a symbolic link, for example:

    $ which gcc
    /usr/bin/gcc
    $ sudo ln -s /usr/bin/gcc /usr/bin/gcc-4.2

--------------------------------------------------------------------------------

[Knight Lab]: http://chem.colorado.edu/knightgroup/
[Python]: http://www.python.org/
[QIIME]: https://github.com/qiime/qiime
[Apache]: httpd.apache.org/
[Xcode 3.2.1]: http://instaar.colorado.edu/KingdomProkarya/software/xcode321_10m2003_developerdvd.dmg
[command line tools for OS X]: http://adcdownload.apple.com/Developer_Tools/command_line_tools_os_x_lion_for_xcode__november_2012/xcode452cltools10_76938212a.dmg
[Xcode]: http://itunes.apple.com/us/app/xcode/id497799835?ls=1&mt=12