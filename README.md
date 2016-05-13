###SiC³

##Aims
- Build a solid foundation for examining reactivity of a particular set of chemical species
- Use a more sophisticated algorithm for mechanism decision-making than in the original SiC
- Create a Web-based graphical user interface that will be friendlier to users than the original SiC

##Dependencies

#Python
- Flask
- openbabel (which requires the OpenBabel package - see http://openbabel.org/docs/2.3.1/Installation/install.html#compile-bindings for details)
- sortedcontainers (to automatically keep choices sorted)

#Java
- ReactionDecoder (https://github.com/asad/ReactionDecoder), included in install, requires modified version.

#Linux Packages for basic install
- python-dev
- apache2
- openbabel
- python-openbabel (so that the Pybel bindings get compiled.)

#Optional Linux Packages required for server
- apache2

##How to run this program
To run the server, run "runserver.py" in the root of the project. Do not attempt to run sic/sigc/server.py directly;
due to the use of relative path imports, this is NOT viable, especially if you want to run the server with debug mode on.
Always use runserver.py and update the settings for the server application there (see Flask docs for details).
