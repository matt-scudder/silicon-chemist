###SiC³

##Aims
- Build a solid foundation for examining reactivity of a particular set of chemical species
- Use a more sophisticated algorithm for mechanism decision-making than in the original SiC
- Create a Web-based graphical user interface that will be friendlier to users than the original SiC

##Dependencies

#Python
- Flask

#Linux Packages (for full WSGI setup; not needed for basic install)
- python-dev
- apache2


##How to run this program
Because SiC³ is set up as a set of separate modules, relative path imports are used, so as to not assume that
all of the modules are installed in python's local libs. However, this leads to the problem that, for example,
neither unit tests nor SiGC will run on their own without throwing an error. To make them work, run them as modules
with an example syntax like this:

python -m sigc.server

This will tell Python to think in terms of "modules" and take advantage of the \_\_init\_\_.py files present in each directory
to organize itself.
