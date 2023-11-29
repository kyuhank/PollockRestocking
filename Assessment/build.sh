#!/bin/bash

set -ex

make

cp *.RData /output/
#cp *.html /output/
