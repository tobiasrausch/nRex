#!/bin/bash

wget 'https://github.com/odelaneau/shapeit4/archive/refs/tags/v4.2.2.tar.gz'
tar -xzf v4.2.2.tar.gz
mv shapeit4-4.2.2/maps/* .
tar -xzf genetic_maps.b37.tar.gz 
tar -xzf genetic_maps.b38.tar.gz
rm -rf genetic_maps.b37.tar.gz genetic_maps.b38.tar.gz v4.2.2.tar.gz shapeit4-4.2.2/
