#!/bin/bash

set -e

# changing permissions to reflect singularityhub
sudo chmod go-rwx -R .
sudo chown root.root -R .

# build image
sudo singularity build travis_test_build.simg Singularity

# revert permissons to previous state (sortof)
sudo chmod 777 -R .
sudo chmod travis.travis -R .

# move image to correct subdirectory
mkdir dynverse
mv travis_test_build.simg dynverse
