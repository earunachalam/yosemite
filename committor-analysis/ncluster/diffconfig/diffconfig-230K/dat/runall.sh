#!/bin/bash --login

find . -name run.sh -exec ./runnative.sh {} \; > runall.log
