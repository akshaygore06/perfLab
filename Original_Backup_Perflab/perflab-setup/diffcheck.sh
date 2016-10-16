echo "Difference Checking Script by Akshay Gore"

echo "Gauss Filter"

python picturediff.py  filtered-gauss-blocks-small.bmp ~/640/PerformanceLab_OriginalBackup/CSCI540-PerfLab/perflab-setup/filtered-gauss-blocks-small.bmp


echo "Average Filter"

python picturediff.py  filtered-avg-blocks-small.bmp ~/640/PerformanceLab_OriginalBackup/CSCI540-PerfLab/perflab-setup/filtered-avg-blocks-small.bmp


echo "Hline Filter"

python picturediff.py  filtered-hline-blocks-small.bmp ~/640/PerformanceLab_OriginalBackup/CSCI540-PerfLab/perflab-setup/filtered-hline-blocks-small.bmp

echo "Emboss Filter"

python picturediff.py  filtered-emboss-blocks-small.bmp ~/640/PerformanceLab_OriginalBackup/CSCI540-PerfLab/perflab-setup/filtered-emboss-blocks-small.bmp



