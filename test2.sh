#!/bin/bash

rm -f universe phprocessed.pdb

phbuilder gentopol -f proteins/1cvo.pdb -m all -v 2

rm -rf \_\_py*
