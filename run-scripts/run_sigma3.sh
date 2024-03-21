#!/bin/bash

#
#

echo "#############################################"
echo "\n\n\n\n"
echo "SIGMA_3 NOT IMPLEMENTED FOR ANYTHING ELSE THAN RESUM BK"
echo "\n\n\n\n"
echo "#############################################"
#
#
./plottool sub resumbk fixedrc z2simple unboundloop | tee -a out/sigma3_comparison/sub_resumbk_fc_z2sim_ubloop.dat
./plottool unsub resumbk fixedrc z2simple unboundloop | tee -a out/sigma3_comparison/unsub_resumbk_fc_z2sim_ubloop.dat
./plottool unsub+ resumbk fixedrc z2simple unboundloop | tee -a out/sigma3_comparison/unsub+_resumbk_fc_z2sim_ubloop.dat
#
#
#./plottool unsub+ lobk fixedrc z2simple unboundloop | tee -a out/unsub+_lobk_fc_z2sim_ubloop.dat
##
#./plottool unsub+ lobk fixedrc z2improved z2boundloop | tee -a out/unsub+_lobk_fc_z2imp_z2loop.dat
#./plottool unsub+ lobk fixedrc z2improved unboundloop | tee -a out/unsub+_lobk_fc_z2imp_ubloop.dat
#./plottool unsub+ lobk fixedrc z2simple z2boundloop | tee -a out/unsub+_lobk_fc_z2sim_z2loop.dat
##
#./plottool unsub+ trbk fixedrc z2improved z2boundloop | tee -a out/unsub+_trbk_fc_z2imp_z2loop.dat
#./plottool unsub+ trbk fixedrc z2improved unboundloop | tee -a out/unsub+_trbk_fc_z2imp_ubloop.dat
#./plottool unsub+ trbk fixedrc z2simple z2boundloop | tee -a out/unsub+_trbk_fc_z2sim_z2loop.dat
##
##
#./plottool unsub+ lobk parentrc z2simple unboundloop | tee -a out/unsub+_lobk_pdrc_z2sim_ubloop.dat
##
#./plottool unsub+ lobk parentrc z2improved z2boundloop | tee -a out/unsub+_lobk_pdrc_z2imp_z2loop.dat
#./plottool unsub+ lobk parentrc z2improved unboundloop | tee -a out/unsub+_lobk_pdrc_z2imp_ubloop.dat
#./plottool unsub+ lobk parentrc z2simple z2boundloop | tee -a out/unsub+_lobk_pdrc_z2sim_z2loop.dat
##
#./plottool unsub+ trbk parentrc z2improved z2boundloop | tee -a out/unsub+_trbk_pdrc_z2imp_z2loop.dat
#./plottool unsub+ trbk parentrc z2improved unboundloop | tee -a out/unsub+_trbk_pdrc_z2imp_ubloop.dat
#./plottool unsub+ trbk parentrc z2simple z2boundloop | tee -a out/unsub+_trbk_pdrc_z2sim_z2loop.dat
##
#./plottool unsub+ lobk guillaumerc z2simple unboundloop | tee -a out/unsub+_lobk_grc_z2sim_ubloop.dat
#
#./plottool unsub+ lobk guillaumerc z2improved z2boundloop | tee -a out/unsub+_lobk_grc_z2imp_z2loop.dat
#./plottool unsub+ lobk guillaumerc z2improved unboundloop | tee -a out/unsub+_lobk_grc_z2imp_ubloop.dat
#./plottool unsub+ lobk guillaumerc z2simple z2boundloop | tee -a out/unsub+_lobk_grc_z2sim_z2loop.dat
#
#./plottool unsub+ trbk guillaumerc z2improved z2boundloop | tee -a out/unsub+_trbk_grc_z2imp_z2loop.dat
#./plottool unsub+ trbk guillaumerc z2improved unboundloop | tee -a out/unsub+_trbk_grc_z2imp_ubloop.dat
#./plottool unsub+ trbk guillaumerc z2simple z2boundloop | tee -a out/unsub+_trbk_grc_z2sim_z2loop.dat
#
