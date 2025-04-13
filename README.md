# SiC<sup>3</sup>

## Aims

- Build a solid foundation for examining reactivity of a particular set of chemical species
- Use a more sophisticated algorithm for mechanism decision-making than in the original SiC
- Create a Web-based graphical user interface that will be friendlier to users than the original SiC

## Development

The fastest way to set up the enviroment will all dependencies is to have [Docker](https://www.docker.com/get-started/) installed and use the included Dockerfile

#### VS Code
This repository contains devcontainer configs, so opening this directory in VS Code will prompt you to install the [Dev Containers](https://marketplace.visualstudio.com/items/?itemName=ms-vscode-remote.remote-containers) extension, and then to reopen the workspace in a container.  The workspace is mounted to `/usr/src/app` so changes made inside the container will persist.  
VS Code will automatically share the Git Credential Manager with the container if you cloned with HTTPS.  
There are debug profiles in `.vscode/launch.json` for the Flask development server and the SiC command line.  
Testing uses `unittest` with the `test_*` prefix.  

#### Docker Compose
Alternatively, there is an included `compose.yml` in the workspace set up to support hot reloading if you're not using VS Code.  
Simply run the following to start the development server:
```bash
docker compose up --watch
```

## Dependencies

### Linux Packages
The following dependencies can all be installed via the linux package manager

- `openjdk-17-jre`
- `libopenbabel-dev`
- `python3-pip`
- `swig`

### Python
The following dependencies are available from [PyPI](https://pypi.org/):

- `flask`
- `sortedcontainers`
- `openbabel`

On many linux distros, libopenbabel-dev does not place its files where the openbabel pip installer espects them, but we can link the folder to the expected location:
```bash
ln -s /usr/include/openbabel3 /usr/local/include/openbabel3
```
  - This requires the `libopenbabel-dev` package to be installed first.

The rest of the dependencies should install with:
```bash
python3 -m pip install -r requirements.txt
```

### Java

- [Reaction Decoder Tool](https://github.com/asad/ReactionDecoder) with dependencies, included in this repository.


## How to run this program
To run the server, run `python3 runserver.py` in the root of the repository.
You can also mark it as executable with `chmod +x runserver.py` and invoke it with `./runserver.py` 
Always use `runserver.py` and update the settings for the server application there (see Flask docs for details).