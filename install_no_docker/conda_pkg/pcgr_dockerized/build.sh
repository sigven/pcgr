#!/usr/bin/env bash

mkdir -p ${PREFIX}/bin
chmod +x ${SRC_DIR}/*.py
mv ${SRC_DIR}/*.py ${PREFIX}/bin/
