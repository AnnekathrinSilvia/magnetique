#!/bin/sh
curl --output data.zip 'https://data.dieterichlab.org/s/3gCTLGT4DaAqaqb/download' -H 'Accept: text/html' --compressed  -H 'Connection: keep-alive'
unzip data.zip