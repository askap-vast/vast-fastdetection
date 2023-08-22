# How to build and run the container.


## Build the docker image:
To run the containerised version you need to create a container image first. You can do this by the following command:
```bash
docker build --rm --no-cache -t vast_fast .
```
## Run the docker image
Once you have got the image built, then you can run the container easily in your local environment. The container 
expects that you will bind the volume (input data directory) while running the image. A sample command will be like the 
following:

```bash
docker run -v </path/to/your/input/data>:/input vast_fast
```

You must have a ```app_settings.ini``` file to describe the settings for this run in the data directory. You can change 
the settings and run the image again to get different processing. The sample content of the ```app_settings.ini``` file 
is given below:

```ini
[PROCESS]
beam = 00
input_directory_name = SB41912_cutout_beam00

[INSTRUMENT]
n_beam = 36

[RESOURCE]
mem = 10g
n_cpu = 4
time = 02:00:00

[EXE]
setup_file = /fred/oz999/tiger/SS2022B-TMurphy/ozstar-dev.sh
python_file = process_beam.py

[RUN_SETTINGS]
OUT_PREFIX = "output"

# list of maps
KTYPELIST = chisquare, peak, std


# Gaussian width
G_WIDTH = 4

# chunksize used in dask.array when generating Gaussian map
# chunksize = (CHUNK_0, CHUNK_1, time_dim), where time_dim (the number of time slices) is determined at run time
# the default value is (100, 100, time_dim), which has been verified during the course of devlopment at ADACS
# NOTICE: if the chunksize is too big, it will result in out-of-memory issue or kill the process
CHUNK_0 = 100
CHUNK_1 = 100

```

# Publish the docker image to a docker registry (Docker Hub example)

To push the Docker image you've built to Docker Hub, you need to follow these steps:

1. **Tag the Image**: Before pushing the image to Docker Hub, you need to tag it with the appropriate Docker Hub repository name. The format for tagging is generally `username/repository:tag`.

2. **Login to Docker Hub**: You need to log in to your Docker Hub account using the `docker login` command.

3. **Push the Image**: After logging in, you can push the tagged image to Docker Hub using the `docker push` command.

Here's how to do it:

1. **Tag the Image**:
   You need to retag your image with the Docker Hub repository name. Replace `username` with your Docker Hub username and `repository` with the name you want for your repository.

   ```bash
   docker tag vast_fast dockerhub_username/vast_fast:latest
   ```

2. **Login to Docker Hub**:
   Run the following command and provide your Docker Hub username and password when prompted:

   ```bash
   docker login
   ```

3. **Push the Image**:
   Push the tagged image to Docker Hub:

   ```bash
   docker push dockerhub_username/vast_fast:latest
   ```

After a successful push, your Docker image will be available on Docker Hub. You can access it using the URL `https://hub.docker.com/r/dockerhub_username/vast_fast`.