# Installing Virtual Machine

These instructions are to set up a virtual machine as a web server where Evident will be running. To access Evident you need to use the [Google Chrome](https://www.google.com/chrome) web browser in your real machine. Basically, you only need to start the virtual machine when you will use Evident and you access it using Chrome from your machine.

## Downloads and Installing
Make sure you have VirtualBox installed in your machine, if not download it and install it from here: [VirtualBox](https://www.virtualbox.org/wiki/Downloads); then download and uncompress this file: [Evident Virtual Machine](ftp://thebeast.colorado.edu/pub/evident-snapshop-VMs/Evident.vdi_100312.tgz) [~ 2GB]

To install the virtual hard drive:

1. Open VirtualBox

2. In the top menu go to VirtualBox -> Preferences -> Network: Hit plus sign on the right and then accept all the changes, this should create a device called vboxnet0, hit OK.

3. Create a new machine (press the New button). A new window will show up. Click ‘Next’. In this screen type Evident as the name for the virtual machine. Then select Linux as the Operating System, and Ubuntu (64 bit) as the version. Click Next. Select the amount of RAM (memory). You will need at least 3G, but the best option is based on your machine. Note that we recommend to use as much as possible. After selecting the amount of RAM, click Next. Select “Use existing hard drive”, and click the folder icon next to the selector (it has a green up arrow). In the new window click ‘Add’, and locate the virtual hard drive that was downloaded before. Click Select and then click Next. In the new window click Finish. Do not start the machine.

4. Go to your newly created VirtualMachine -> Settings -> Network, then

	Adapter 1 -> set "Attached to" menu to HostOnly -> select the device created in step 2

	Adapter 2 -> select "Enable Network Adapter" -> set "Attached to" menu to NAT

5. Now start your machine and wait for it to boot, it will take a while because it's configuring all the hardware. The username and password are: `evident`. To connect to the Evident frontpage/website, try opening this address in Google Chrome in your real machine: [http://192.168.56.101/evident/](http://192.168.56.101/evident/), if this doesn't work, it means the configuration changed the IP address. To find the new IP: start a terminal session in the virtual machine, run: `ifconfig | grep eth -A 2`, find `eth0`, and use the IP address listed there. Note, that you should only have to do this once.