FROM shinyverse_criu_base:latest

WORKDIR /root/magnetique/
# setup app
COPY MAGNetApp.rdata /root/magnetique/
COPY app.R /root/magnetique/
COPY utils.R /root/magnetique/
COPY database.R /root/magnetique/
COPY www/ /root/magnetique/www/
COPY tours/ /root/magnetique/tours/

COPY deploy/rstart.R rstart.R
COPY deploy/restore.sh restore.sh
COPY deploy/launch_shiny_app.sh launch_shiny_app.sh
RUN mkdir -p /criu_dumps \
    && chmod 755 restore.sh launch_shiny_app.sh
EXPOSE 3838
CMD /bin/bash /restore.sh
