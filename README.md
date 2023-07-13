# SiC³

## Aims

- Build a solid foundation for examining reactivity of a particular set of chemical species
- Use a more sophisticated algorithm for mechanism decision-making than in the original SiC
- Create a Web-based graphical user interface that will be friendlier to users than the original SiC

## Dependencies

### Linux Packages for basic install
The following dependencies can all be installed via the linux package manager

- `python2-dev`
- `apache2`
- `openbabel`
- `libopenbabel-dev`
- `swig`
- `default-jre`

### Python
The following dependencies are available from [PyPI](https://pypi.org/) and installable via `sudo pip2 install` on linux

- `flask`
- `sortedcontainers`
- `openbabel` (requires the `libopenbabel-dev` and `swig` packages to be installed)

### Java

- [Reaction Decoder Tool](https://github.com/asad/ReactionDecoder) with dependencies, included in this repository.


## How to run this program
To run the server, run "python2 runserver.py" in the directory which contains it.
Do not attempt to run sic/sigc/server.py directly; due to the use of relative path imports, this is NOT viable, especially if you want to run the server with debug mode on.
Always use runserver.py and update the settings for the server application there (see Flask docs for details).