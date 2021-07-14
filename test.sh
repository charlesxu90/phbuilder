#!/bin/bash

rm -f universe phprocessed.pdb

phbuilder gentopol -f proteins/1cvo.pdb -m all -v 3

rm -rf __pychache__
