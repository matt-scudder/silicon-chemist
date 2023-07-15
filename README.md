# SiC<sup>3</sup>

## Aims

- Build a solid foundation for examining reactivity of a particular set of chemical species
- Use a more sophisticated algorithm for mechanism decision-making than in the original SiC
- Create a Web-based graphical user interface that will be friendlier to users than the original SiC

## Dependencies

### Linux Packages
The following dependencies can all be installed via the linux package manager

- `apache2`
- `default-jre`
- `python3-openbabel`
- `python3-pip`

### Python
The following dependencies are available from [PyPI](https://pypi.org/) and installable via `python3 -m pip install` on linux

- `flask`
- `sortedcontainers`

### Java

- [Reaction Decoder Tool](https://github.com/asad/ReactionDecoder) with dependencies, included in this repository.


## How to run this program
To run the server, run `python3 runserver.py` in the root of the repository.
You can also mark it as executable with `chmod +x runserver.py` and invoke it with `./runserver.py` 
Always use runserver.py and update the settings for the server application there (see Flask docs for details).