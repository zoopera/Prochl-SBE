#!/bin/bash

pys=($(find $1 -name "*py"))

for py in "${pys[@]}"
do
	sub=$(echo $py|sed 's/.*\///g'|sed 's/py/lsf/g')
	out=$(echo $py|sed 's/.*\///g'|sed 's/py/p4/g')

	echo -e "#!/bin/bash -l\n\n#SBATCH -n 1\n" > $sub
#	echo -e "source activate p4_env\n" >> $sub
	echo -e "time p4 $py &>$out" >> $sub
#	echo -e "source deactivate\n" >> $sub

	chmod +x $sub
	sbatch $sub
done
