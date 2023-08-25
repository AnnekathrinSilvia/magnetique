#!/bin/bash
export RENV_WATCHDOG_ENABLED=0
R -e 'renv::init(bioconductor=TRUE); renv::restore(prompt = FALSE); renv::isolate()'

echo "Running renv $(pidof  R)"
setsid R -f rstart.R > /app/rlogs.log 2>&1 &
while ! grep -q -F 'Library loaded' /app/rlogs.log
do
  echo 'Waiting for library load ...'
  sleep 5
done
rpid=$(pidof  R | sed -e "s/ / -t /g")
echo "Running rstar $(pidof  R)"

rm -rf /criu_dumps/*

criu-ns dump -t $rpid -vvv -o /criu_dumps/dump.log -D /criu_dumps/ && echo OK
touch /criu_dumps/dump.done
read -p 'Waiting for input'
