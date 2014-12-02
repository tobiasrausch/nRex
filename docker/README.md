Dockerized nRex
===============

This folder is used to create a Docker image of nRex. Installation and usage instruction can be found on the [docker](https://www.docker.io) website.
Once you have docker installed you can just pull the nRex image:

`sudo docker pull trausch/nrex`

You can then run nRex from that image. Below we assume your bam files are in /var/data which is mounted as /data in the docker container:

`sudo docker run -i -t -v /var/data:/data trausch/nrex /bin/bash`

`./nRex/nRex.sh /data/ref.fa /data/tumor.bam /data/normal.bam`

Once nRex is finished you can just exit the Docker image:

`exit`
