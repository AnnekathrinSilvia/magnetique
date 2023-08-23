To build the image:

    ./deploy/build.sh

To run the image:

    docker run --privileged -it -e "RENV_PATHS_CACHE=${RENV_PATHS_CACHE_CONTAINER}" -v "${RENV_PATHS_CACHE_HOST}:${RENV_PATHS_CACHE_CONTAINER}" magnetique "./restore.sh"

After updating a R package (with renv::install or renv::record):

    docker cp $container:/root/magnetique/renv.lock renv.lock

Dev run:
    ```
    docker  run -p 68:3838 --rm --privileged -it -e "RENV_PATHS_CACHE=${RENV_PATHS_CACHE_CONTAINER}" -v "${RENV_PATHS_CACHE_HOST}:${RENV_PATHS_CACHE_CONTAINER}" magnetique ./restore.sh
    ```
Logs location:
    R log  logs: /criu_dumps/rlogs.log
    restore log: /criu_dumps/restore.log



