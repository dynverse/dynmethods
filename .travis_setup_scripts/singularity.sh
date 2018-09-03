#!/bin/bash

set -e

export GOPATH=${HOME}/go
export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin

export SINGULARITY_DIR=$GOPATH/src/github.com/singularityware/singularity
echo $SINGULARITY_DIR

if [ "$TRAVIS_OS_NAME" == "osx" ]; then # use homebrew version
	echo "Panic!"
else
	if [ -f $SINGULARITY_DIR/builddir/singularity ]; then
	    echo "using cached build"
    else
		sudo apt-get update
		sudo apt-get install -y build-essential libssl-dev uuid-dev libgpgme11-dev

		# Install go
		export VERSION=1.10.3 OS=linux ARCH=amd64
		wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz
		sudo tar -C /usr/local -xzf go$VERSION.$OS-$ARCH.tar.gz

		# Not sure why this is here
		# sudo sed -i -e 's/^Defaults\tsecure_path.*$//' /etc/sudoers
		
		# Clone singularity repo
		mkdir -p $GOPATH/src/github.com/singularityware
		pushd $GOPATH/src/github.com/singularityware
		git clone https://github.com/singularityware/singularity.git
		cd singularity

		# Install golang dependencies
		go get -u -v github.com/golang/dep/cmd/dep

		# Build Singularity (finally)
		pushd $SINGULARITY_DIR
		./mconfig

		popd
		popd
	fi

	# Install Singularity
	pushd $SINGULARITY_DIR/builddir
	sudo make install
	popd
fi

