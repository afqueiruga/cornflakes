# Superbuild and Docker image

## Building the docker image


The included `dockerfile` sets up a Debian based docker image. Build
it inside of the cornflakes directory:
```bash
docker build -f ./superbuild/Dockerfile -t cornflakes .
```

Once the docker image is in your local registry, the environment can be loaded with
```bash
docker run -ti -v `pwd`:/home/user cornflakes
```
You will have acess to the current directory of your shell inside of
the docker image.

The docker image is also hosted in the public repositories where it is automatically built agains the github repository, so you can just directly
```bash
docker pull afqu/cornflakes
```
without having to deal with the source code at all! The name is `afqu/cornflakes` instead of just `cornflakes`.
