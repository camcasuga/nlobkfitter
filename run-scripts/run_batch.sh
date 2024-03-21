#!/bin/bash

./plottool sub lobk fixedrc z2simple unboundloop | tee -a out/sub_lobk_fc_z2sim_ubloop.dat
./plottool unsub lobk fixedrc z2simple unboundloop | tee -a out/unsub_lobk_fc_z2sim_ubloop.dat
./plottool sub lobk fixedrc z2improved z2boundloop | tee -a out/sub_lobk_fc_z2imp_z2loop.dat
./plottool sub lobk fixedrc z2improved unboundloop | tee -a out/sub_lobk_fc_z2imp_ubloop.dat
./plottool sub lobk fixedrc z2simple z2boundloop | tee -a out/sub_lobk_fc_z2sim_z2loop.dat
#
./plottool unsub lobk fixedrc z2improved z2boundloop | tee -a out/unsub_lobk_fc_z2imp_z2loop.dat
./plottool unsub lobk fixedrc z2improved unboundloop | tee -a out/unsub_lobk_fc_z2imp_ubloop.dat
./plottool unsub lobk fixedrc z2simple z2boundloop | tee -a out/unsub_lobk_fc_z2sim_z2loop.dat
#
./plottool sub trbk fixedrc z2improved z2boundloop | tee -a out/sub_trbk_fc_z2imp_z2loop.dat
./plottool sub trbk fixedrc z2improved unboundloop | tee -a out/sub_trbk_fc_z2imp_ubloop.dat
./plottool sub trbk fixedrc z2simple z2boundloop | tee -a out/sub_trbk_fc_z2sim_z2loop.dat
#
./plottool unsub trbk fixedrc z2improved z2boundloop | tee -a out/unsub_trbk_fc_z2imp_z2loop.dat
./plottool unsub trbk fixedrc z2improved unboundloop | tee -a out/unsub_trbk_fc_z2imp_ubloop.dat
./plottool unsub trbk fixedrc z2simple z2boundloop | tee -a out/unsub_trbk_fc_z2sim_z2loop.dat
#
#
./plottool sub lobk parentrc z2simple unboundloop | tee -a out/sub_lobk_pdrc_z2sim_ubloop.dat
./plottool unsub lobk parentrc z2simple unboundloop | tee -a out/unsub_lobk_pdrc_z2sim_ubloop.dat
./plottool sub lobk parentrc z2improved z2boundloop | tee -a out/sub_lobk_pdrc_z2imp_z2loop.dat
./plottool sub lobk parentrc z2improved unboundloop | tee -a out/sub_lobk_pdrc_z2imp_ubloop.dat
./plottool sub lobk parentrc z2simple z2boundloop | tee -a out/sub_lobk_pdrc_z2sim_z2loop.dat
#
./plottool unsub lobk parentrc z2improved z2boundloop | tee -a out/unsub_lobk_pdrc_z2imp_z2loop.dat
./plottool unsub lobk parentrc z2improved unboundloop | tee -a out/unsub_lobk_pdrc_z2imp_ubloop.dat
./plottool unsub lobk parentrc z2simple z2boundloop | tee -a out/unsub_lobk_pdrc_z2sim_z2loop.dat
#
./plottool sub trbk parentrc z2improved z2boundloop | tee -a out/sub_trbk_pdrc_z2imp_z2loop.dat
./plottool sub trbk parentrc z2improved unboundloop | tee -a out/sub_trbk_pdrc_z2imp_ubloop.dat
./plottool sub trbk parentrc z2simple z2boundloop | tee -a out/sub_trbk_pdrc_z2sim_z2loop.dat
#
./plottool unsub trbk parentrc z2improved z2boundloop | tee -a out/unsub_trbk_pdrc_z2imp_z2loop.dat
./plottool unsub trbk parentrc z2improved unboundloop | tee -a out/unsub_trbk_pdrc_z2imp_ubloop.dat
./plottool unsub trbk parentrc z2simple z2boundloop | tee -a out/unsub_trbk_pdrc_z2sim_z2loop.dat
#
#
./plottool sub lobk guillaumerc z2simple unboundloop | tee -a out/sub_lobk_grc_z2sim_ubloop.dat
./plottool unsub lobk guillaumerc z2simple unboundloop | tee -a out/unsub_lobk_grc_z2sim_ubloop.dat
./plottool sub lobk guillaumerc z2improved z2boundloop | tee -a out/sub_lobk_grc_z2imp_z2loop.dat
./plottool sub lobk guillaumerc z2improved unboundloop | tee -a out/sub_lobk_grc_z2imp_ubloop.dat
./plottool sub lobk guillaumerc z2simple z2boundloop | tee -a out/sub_lobk_grc_z2sim_z2loop.dat
#
./plottool unsub lobk guillaumerc z2improved z2boundloop | tee -a out/unsub_lobk_grc_z2imp_z2loop.dat
./plottool unsub lobk guillaumerc z2improved unboundloop | tee -a out/unsub_lobk_grc_z2imp_ubloop.dat
./plottool unsub lobk guillaumerc z2simple z2boundloop | tee -a out/unsub_lobk_grc_z2sim_z2loop.dat
#
./plottool sub trbk guillaumerc z2improved z2boundloop | tee -a out/sub_trbk_grc_z2imp_z2loop.dat
./plottool sub trbk guillaumerc z2improved unboundloop | tee -a out/sub_trbk_grc_z2imp_ubloop.dat
./plottool sub trbk guillaumerc z2simple z2boundloop | tee -a out/sub_trbk_grc_z2sim_z2loop.dat
#
./plottool unsub trbk guillaumerc z2improved z2boundloop | tee -a out/unsub_trbk_grc_z2imp_z2loop.dat
./plottool unsub trbk guillaumerc z2improved unboundloop | tee -a out/unsub_trbk_grc_z2imp_ubloop.dat
./plottool unsub trbk guillaumerc z2simple z2boundloop | tee -a out/unsub_trbk_grc_z2sim_z2loop.dat
#
