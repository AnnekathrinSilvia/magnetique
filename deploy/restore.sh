#!/bin/bash
tail -f --retry  /app/rlogs_container.log &
echo $SHINYPROXY_USERNAME > /shiny_user.txt
echo $SHINYPROXY_USERGROUPS > /shiny_usergroups.txt
touch /restore.marker
criu-ns restore -vvvv -D /criu_dumps -o /app/restore.log