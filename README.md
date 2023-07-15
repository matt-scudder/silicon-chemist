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
- `libopenbabel-dev`
- `python3-pip`

### Python
The following dependencies are available from [PyPI](https://pypi.org/) and on linux can be installed by running the following from the repository root:  
```bash
python3 -m pip install -r requirements.txt
```

- `flask`
- `sortedcontainers`
- `openbabel` 
  - This requires the `libopenbabel-dev` package to be installed first.
  - If you are installing this manually, you may need to provide the path to openbabel's header files to the pip installer.
    - Check `requirements.txt` for the syntax, as it adds the path `libopenbabel-dev` unpacks them to on Ubuntu/Debian

### Java

- [Reaction Decoder Tool](https://github.com/asad/ReactionDecoder) with dependencies, included in this repository.


## How to run this program
To run the server, run `python3 runserver.py` in the root of the repository.
You can also mark it as executable with `chmod +x runserver.py` and invoke it with `./runserver.py` 
Always use `runserver.py` and update the settings for the server application there (see Flask docs for details).