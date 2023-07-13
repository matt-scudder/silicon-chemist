###SiC³

##Aims
- Build a solid foundation for examining reactivity of a particular set of chemical species
- Use a more sophisticated algorithm for mechanism decision-making than in the original SiC
- Create a Web-based graphical user interface that will be friendlier to users than the original SiC

##Dependencies

#Python (installed through pip)
- flask
- openbabel (requires the libopenbabel-dev package)
- sortedcontainers (to automatically keep choices sorted)

#Java
- ReactionDecoder (https://github.com/asad/ReactionDecoder), included in install.

#Linux Packages for basic install
- python2-dev
- apache2
- libopenbabel-dev
- swig
- default-jre

#Optional Linux Packages required for server
- apache2

##How to run this program
To run the server, run "runserver.py" in the root of the project. Do not attempt to run sic/sigc/server.py directly;
due to the use of relative path imports, this is NOT viable, especially if you want to run the server with debug mode on.
Always use runserver.py and update the settings for the server application there (see Flask docs for details).
