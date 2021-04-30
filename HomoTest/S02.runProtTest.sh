#!/bin/bash

phy=($(find $1 -maxdepth 1 -type f -name "*trimAl"))

for input in "${phy[@]}"
do
	sub=$(echo $input | sed 's/.*\///g' | sed 's/trimAl/sh/g')
	out=$(echo $input | sed 's/trimAl/tab/g')
	echo -e "#!/bin/bash -l\n\n" > $sub
	echo -e "time java -jar /home-user/software/prottest-3.4.2/prottest-3.4.2.jar -i $input -o $out -all-distributions -F -BIC -tc 0.5 -threads 4\n" >> $sub

	chmod +x $sub
	sbatch -n4 $sub
done

