#!/bin/sh
wget --output-document data.zip 'https://data.dieterichlab.org/s/3gCTLGT4DaAqaqb/download' --header 'Accept: text/html' --header 'Connection: keep-alive'
unzip -o data.zip
rm data.zip