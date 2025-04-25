#!/bin/bash

cd PPARg
cp ../min_grompp_mdrun.sh .

echo "#!/bin/bash" > temp
echo "#SBATCH --job-name=PPARg_min" >> temp
cat min_grompp_mdrun.sh >> temp
mv temp min_grompp_mdrun.sh

sbatch min_grompp_mdrun.sh
cd ..
