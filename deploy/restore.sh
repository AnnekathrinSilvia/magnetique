#!/bin/bash
tail -f --retry  /criu_dumps/rlogs_container.log &
echo $SHINYPROXY_USERNAME > /shiny_user.txt
echo $SHINYPROXY_USERGROUPS > /shiny_usergroups.txt
touch /restore.marker
criu-ns restore -vvvv -D /criu_dumps -o /criu_dumps/restore.log