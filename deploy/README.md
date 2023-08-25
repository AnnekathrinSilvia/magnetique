To build the image:
    
    rm -f criu_dumps/  
    ./deploy/build.sh
    sudo cp criu_dumps/* /tmp_shiny/magnetique_criu_dump

After updating a R package (with renv::install or renv::record):

    docker cp $container:/root/magnetique/renv.lock renv.lock

Check logs:
    docker run magnetique cat rlogs.log
    
