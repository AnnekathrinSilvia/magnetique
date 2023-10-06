#!/bin/bash

export RENV_PATHS_CACHE_HOST=~/renv/cache
export RENV_PATHS_CACHE_CONTAINER=/renv/cache


mkdir -p $RENV_PATHS_CACHE_HOST

docker build -f deploy/Dockerfile_base --progress=plain -t shinyverse_criu_base .
docker build -f deploy/Dockerfile --progress=plain -t magnetique_prefreeze:latest .

container=$(docker run -d --rm -i --privileged \
    -e "RENV_PATHS_CACHE=${RENV_PATHS_CACHE_CONTAINER}" \
    -v "${RENV_PATHS_CACHE_HOST}:${RENV_PATHS_CACHE_CONTAINER}" \
    -v "$(pwd)/criu_dumps/:/criu_dumps/" \
    -v "$(pwd)/data/:/data/" \
    magnetique_prefreeze:latest /bin/bash launch_shiny_app.sh
)
dump_done=''
while [[ $dump_done != "/criu_dumps/dump.done" ]]
do
    dump_done=$(docker exec "$container" ls /criu_dumps/dump.done)
    echo "Waiting for dump..."
    sleep 5
done
docker commit "$container" magnetique:latest
docker container stop "$container"
