#!/bin/bash

convert 4HFI_4_0*.png +append A.png
convert 4HFI_7_0*.png +append B.png
convert 6ZGD_4_0*.png +append C.png
convert 6ZGD_7_0*.png +append D.png
# convert 6ZGD_4_long_*.png +append E.png
# convert 6ZGD_7_long_*.png +append F.png

convert A.png B.png C.png D.png -append total_1.png
# convert C.png D.png E.png F.png -append total_2.png

rm -f A.png B.png C.png D.png E.png F.png
