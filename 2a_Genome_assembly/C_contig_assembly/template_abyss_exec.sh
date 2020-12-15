#!/bin/bash

abyss-pe name=SP_abyss_kKMER j=32 k=KMER \
 lib='pea peb pec ped pef' mp='mpa mpb' \
 se='pse pse_run2 mse' \
 pea='SP_R1t_is_225.fq.gz SP_R2t_is_225.fq.gz' \
 peb='SP_R1t_is_350.fq.gz SP_R2t_is_350.fq.gz' \
 pec='SP_R1t_is_350_run2.fq.gz SP_R2t_is_350_run2.fq.gz' \
 ped='SP_R1t_is_550.fq.gz SP_R2t_is_550.fq.gz' \
 pef='SP_R1t_is_700.fq.gz SP_R2t_is_700.fq.gz' \
 mpa='mp_nxtrim_FR/SP_R1t_is_3000.fq.gz mp_nxtrim_FR/SP_R2t_is_3000.fq.gz' \
 mpb='mp_nxtrim_FR/SP_R1t_is_5000.fq.gz mp_nxtrim_FR/SP_R2t_is_5000.fq.gz' \
 pse='SP_R1np_is_350.fq.gz SP_R1np_is_550.fq.gz SP_R1np_is_700.fq.gz SP_R2np_is_350.fq.gz SP_R2np_is_550.fq.gz SP_R2np_is_700.fq.gz' \
 pse_run2='SP_R1np_is_350_run2.fq.gz SP_R2np_is_350_run2.fq.gz' \
 mse='SP_se_mp.fq.gz' \
