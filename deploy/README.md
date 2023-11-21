On the server, fetch the newest version from github:
  
  git pull
  
If required, build the new database with `database.R` or copy it to the server.

To build the image:
    
    rm -f criu_dumps/  
    ./deploy/build.sh
    sudo cp criu_dumps/* /tmp_shiny/magnetique_criu_dump


and then restart the server
  sudo systemctl restart shinyproxy-prod

# 
After updating a R package (with renv::install or renv::record):

    docker cp $container:/root/magnetique/renv.lock renv.lock

Check logs:
    docker run magnetique cat rlogs.log
    
