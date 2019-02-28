#!/bin/sh
sudo apt-get update
sudo apt-get install -y gmsh
sudo pip3 install Image
sudo python3 SimDataDB/setup.py install
