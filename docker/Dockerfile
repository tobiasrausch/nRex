# nRex Dockerfile
#
# Version 0.0.1
# build with: sudo docker build --rm --no-cache -t nrex - < Dockerfile
# run with: sudo docker run -i -t nrex /bin/bash

# use the ubuntu base image
# FROM ubuntu:13.10
FROM ubuntu

MAINTAINER Tobias Rausch rausch@embl.de

# update package repository
RUN apt-get update

# install g++, git, zlib, cmake, boost, ...
RUN apt-get install -y build-essential g++ git cmake zlib1g-dev ncurses-dev openjdk-7-jdk ant wget

# install delly
RUN git clone --recursive https://github.com/tobiasrausch/nRex.git
RUN cd nRex/ ; make all
