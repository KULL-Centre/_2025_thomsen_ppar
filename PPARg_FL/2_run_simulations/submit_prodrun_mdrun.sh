#!/bin/bash

for rep in 1 2
do

cd PPARg_rep${rep}
cp ../prodrun_mdrun.sh .

echo "#!/bin/bash" > temp
echo "#SBATCH --job-name=PPARg_${rep}_md" >> temp
cat prodrun_mdrun.sh >> temp
mv temp prodrun_mdrun.sh

sbatch prodrun_mdrun.sh
cd ..

done
