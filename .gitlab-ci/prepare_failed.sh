#!/bin/bash

rm -rf nektar/build/dist
tar czf nektar.tar.gz nektar
rm -rf nektar
