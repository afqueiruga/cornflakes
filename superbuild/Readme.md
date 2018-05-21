# Superbuild and Docker image

## Building the docker image


The included `dockerfile` sets up a Debian based docker image. Build
it just outside of the cornflakes directory:
```bash
docker build -f cornflakes/superbuild/Dockerfile -t cornflakes .
```


Once the docker image is in your local registry, the environment can be loaded with
```bash
docker run -ti -v `pwd`:/home/user cornflakes
```
You will have acess to the current directory of your shell inside of
the docker image.
