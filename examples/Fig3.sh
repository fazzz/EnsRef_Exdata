#!/bin/sh

if [ ! -f J3_HN_HA_ave.txt ]; then
    python J3_HN_HA.py ../../production_abci/production.xtc \
                       ../../production_abci/production.pdb \
		       J3_HN_HA list_pairs_6LYT.txt
fi


p=4.5
if [ ! -f tmp/${p}.txt -o ! -f tmp/w_${p}.txt -o ! -f tmp/o_${p}.txt ]; then
    ./ensRefEx --Del ${p} \
	       exp_6LYT.txt \
	       J3_HN_HA_trj.txt \
	       tmp/w_${p}.txt tmp/o_${p}.txt > tmp/${p}.txt
fi

diff_tgt=$( grep "Initial Value" tmp/${p}.txt | awk '{printf("%8.3f ",$6)}' )

diff_tol=$( ./merge_sim_exp.sh J3_HN_HA_ave.txt J3_HN_HA_experment_1991_Biochem.txt \
     | awk 'BEGIN{sum=0.0}{sum=sum+$3;}END{printf("%8.3f\n",sum)}' )
printf  "%8.3f %8.3f %8.3f %8.3f %8.3f\n" 0.0 ${diff_tgt} ${diff_tol} 0.0 0.0 \
	> div_vs_chi.txt

if [ ! -d tmp ]; then
    mkdir tmp
fi

for p in 4.0 3.5 3.0 2.5 2.0 1.5 1.0; do
    
    if [ ! -f tmp/${p}.txt -o ! -f tmp/w_${p}.txt -o ! -f tmp/o_${p}.txt ]; then
	./ensRefEx --Del ${p} \
		   exp_6LYT.txt \
		   J3_HN_HA_trj.txt \
		   tmp/w_${p}.txt tmp/o_${p}.txt > tmp/${p}.txt
    fi

    kld=$( cat tmp/${p}.txt | grep Final | awk '{printf("%8.3f ",$9)}' )
    chis_tgt=$( cat tmp/${p}.txt | grep Final | awk '{printf("%8.3f ",$6)}' )
    
    python J3_HN_HA_weight.py ../../production_abci/production.xtc \
                              ../../production_abci/production.pdb \
			      J list_pairs_6LYT.txt tmp/w_${p}.txt

    chis_tol=$( ./merge_sim_exp.sh J_wave.txt J3_HN_HA_experment_1991_Biochem.txt \
		    | awk 'BEGIN{sum=0.0}{sum=sum+$3;}END{print sum}' )

    delta_tgt=$( echo "scale=2.0;${diff_tgt}-${chis_tgt}" | bc )
    delta_tol=$( echo "scale=2.0;${diff_tol}-${chis_tol}" | bc )

    printf  "%8.3f %8.3f %8.3f %8.3f %8.3f\n" \
	    ${kld} ${chis_tgt} ${chis_tol} ${delta_tgt} ${delta_tol}
    
done >> div_vs_chi.txt

#gnuplot div_vs_chi.gp
