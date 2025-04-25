#!/bin/bash

cp -r PPARg PPARg_rep2
mv PPARg PPARg_rep1

for rep in 1 2
do

cd PPARg_rep${rep}
cp ../relax_grompp_mdrun.sh .

echo "#!/bin/bash" > temp
echo "#SBATCH --job-name=PPARg_${rep}_eq" >> temp
cat relax_grompp_mdrun.sh >> temp
mv temp relax_grompp_mdrun.sh

sbatch relax_grompp_mdrun.sh
cd ..

done
