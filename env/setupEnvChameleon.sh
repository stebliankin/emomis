#!/bin/bash


# Install basic tools for compiling
sudo yum groupinstall -y 'Development Tools'
# Ensure EPEL repository is available
sudo yum install -y epel-release
# Install RPM packages for dependencies
sudo yum install -y \
    libseccomp-devel \
    squashfs-tools \
    cryptsetup

export VERSION=1.16.7 OS=linux ARCH=amd64  # change this as you need

# Install GO
wget -O /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz \
  https://dl.google.com/go/go${VERSION}.${OS}-${ARCH}.tar.gz
sudo tar -C /usr/local -xzf /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz
echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
source ~/.bashrc

# Install Singularity
git clone https://github.com/sylabs/singularity.git
cd singularity
git checkout v3.8.1

# Compile
./mconfig
make -C builddir
sudo make -C builddir install

sudo yum install nano

cd /
sudo git clone https://github.com/stebliankin/emomis-dev
cd /emomis-dev/env
sudo singularity build emomis.sif emomis.def

cp emomis.sif $HOME/my_mounting_point/singularity_containers/