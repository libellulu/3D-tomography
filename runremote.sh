#!/usr/bin/env bash
sshpass -p '4fastdata' rsync -avz ./ lbechtet@betelgeuse.ipfn.tecnico.ulisboa.pt:/home/lbechtet/3D-tomography/

sshpass -p '4fastdata' ssh lbechtet@betelgeuse.ipfn.tecnico.ulisboa.pt -C "/bin/bash -i -c 'cd /home/lbechtet/3D-tomography && python3.7 $1'"
