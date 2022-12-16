# How to use the SIAC Docker container

## Modified from ChatGPT

### Feng Yin
### Department of Geography, UCL
### ucfafyi@ucl.ac.uk

![Docker Image Size (latest by date)](https://img.shields.io/docker/image-size/marcyin/siac)

![](siac_docker/siac_docker_diag.png)

The [`marcyin/siac`](https://hub.docker.com/r/marcyin/siac) Docker container provides a convenient way to run the SIAC software on any system that has Docker installed. This can be especially useful if you don't want to install the dependencies required to run SIAC on your local system.

To use the SIAC Docker container, follow these steps:

1. Make sure you have Docker installed on your system. If you don't have Docker, you can [download and install it here](https://www.docker.com/get-started/).
   
2. Pull the `marcyin/siac` Docker image from Docker Hub by running the following command:
```bash
docker pull marcyin/siac
```

1. Once the image is downloaded, you can run the SIAC Docker container with the following command:
```bash
docker run --rm -it marcyin/siac bash
```
This will start the SIAC Docker container and drop you into a shell inside the container.

1. To mount external drive folders to the container, you can use the `-v` flag to bind mount the host directories to the container. For example, to mount the `MCD43`, `DEM`, `water_mask`, `S2_L1C`, and `S2_L2A` directories, you can use the following command:
```bash
docker run --rm --name SIAC -it -v /path/to/MCD43:/MCD43 -v /path/to/DEM:/DEM -v /path/to/water_mask:/water_mask -v /path/to/S2_L1C:/S2_L1C -v /path/to/S2_L2A:/S2_L2A marcyin/siac
```

Replace `/path/to/MCD43`, `/path/to/DEM`, `/path/to/water_mask`, `/path/to/S2_L1C`, and `/path/to/S2_L2A` with the actual paths to the directories on your host system.

The `--name` flag allows you to specify a custom name for the Docker container. In this case, we've named it `SIAC`. This can be useful if you need to reference the container by name later on, for example when you want to stop or delete it.


5. Once you're inside the SIAC Docker container, you can run the SIAC commands as you would normally do. When you're done using the container, you can simply exit the shell to stop the container. The `--rm` flag specified in the `docker run` command will automatically delete the container when it's stopped, so you don't have to worry about cleaning up after yourself.