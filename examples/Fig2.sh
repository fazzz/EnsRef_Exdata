#!/bin/sh

if [ ! -f J3_HN_HA_ave.txt ]; then
    python J3_HN_HA.py ../../production_abci/production.xtc \
                       ../../production_abci/production.pdb \
		       J3_HN_HA list_pairs_6LYT.txt
fi

./merge_sim_exp.sh J3_HN_HA_ave.txt \
		   J3_HN_HA_experment_1991_Biochem.txt \
		   > J3_HN_HA_6LYT_100nsave_vs_exp.txt

cat list_pairs_6LYT.txt | while read line; do
    grep ${line} J3_HN_HA_6LYT_100nsave_vs_exp.txt ;
done > J3_HN_HA_6LYT_100nsave_vs_exp_select.txt

if [ ! -f tmp/2.5.txt -o ! -f tmp/w_2.5.txt -o ! -f tmp/o_2.5.txt ]; then
	./ensRefEx --Del 2.5 \
		   exp_6LYT.txt \
		   J3_HN_HA_trj.txt \
		   tmp/w_2.5.txt tmp/o_2.5.txt > tmp/2.5.txt
fi

python J3_HN_HA_weight.py ../../production_abci/production.xtc \
                          ../../production_abci/production.pdb \
			  J list_pairs_6LYT.txt tmp/w_2.5.txt

./merge_sim_exp.sh J_wave.txt J3_HN_HA_experment_1991_Biochem.txt \
    > J3_HN_HA_6LYT_100nswave_vs_exp.txt

cat list_pairs_6LYT.txt | while read line; do
    grep ${line} J3_HN_HA_6LYT_100nswave_vs_exp.txt ;
done > J3_HN_HA_6LYT_100nswave_vs_exp_select.txt

gnuplot J3_HN_HA_6LYT_100nsopt_vs_exp.gp
