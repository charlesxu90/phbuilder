#!/bin/bash

# E35, E35-Na+, E35-T158

convert proto_6ZGD_7_E35.png       proto_6ZGD_4_E35.png       proto_4HFI_7_E35.png       proto_4HFI_4_E35.png       +append A.png
convert mindst_6ZGD_7_E35-Na+.png  mindst_6ZGD_4_E35-Na+.png  mindst_4HFI_7_E35-Na+.png  mindst_4HFI_4_E35-Na+.png  +append B.png
convert mindst_6ZGD_7_E35-T158.png mindst_6ZGD_4_E35-T158.png mindst_4HFI_7_E35-T158.png mindst_4HFI_4_E35-T158.png +append C.png

convert A.png B.png C.png -append panel_E35.png; rm A.png B.png C.png

convert proto_6ZGD_7_E243.png       proto_6ZGD_4_E243.png       proto_4HFI_7_E243.png       proto_4HFI_4_E243.png       +append A.png
convert mindst_6ZGD_7_E243-K248.png mindst_6ZGD_4_E243-K248.png mindst_4HFI_7_E243-K248.png mindst_4HFI_4_E243-K248.png +append B.png
convert rmsd_6ZGD_7_246-252.png     rmsd_6ZGD_4_246-252.png     rmsd_4HFI_7_246-252.png     rmsd_4HFI_4_246-252.png     +append C.png

convert A.png B.png C.png -append panel_E243.png; rm A.png B.png C.png
