#!/bin/bash

set -e

if [ -d /docker-entrypoint ]; then
    mkdir -p $HOME/data
    cp -r /docker-entrypoint/* $HOME/data
fi

exec "$@"
