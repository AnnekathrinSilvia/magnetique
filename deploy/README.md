To build the image:

    ./deploy/build.sh

To run the image:

    docker run --privileged  -it  magnetique "./restore.sh"

After updating a R package:

    docker cp $container:/root/magnetique/renv.lock renv.lock


Logs location:
    R log  logs: /criu_dumps/rlogs.log
    restore log: /criu_dumps/restore.log
