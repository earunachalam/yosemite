#!/bin/bash --login

cd $(dirname $1)
qsub run.sh
