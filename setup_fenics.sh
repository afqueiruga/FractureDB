#!/bin/sh
sudo apt-get update
sudo apt-get install -y gmsh
sudo pip3 install Image
cd oss/SimDataDB ; sudo python3 setup.py install ; cd ../..
