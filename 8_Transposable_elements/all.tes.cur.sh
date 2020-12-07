#!/bin/bash

#echo "###script to reannotate using PASTEC and Tce balst evidence using args 1 = input and args 2 = output and args 3 = species abbreviation (eg Tce)###"

echo "$1 $2"

cat $1 | \
sed 's/DHX-comp$/DHH#RC\/Helitron/' | \
sed 's/DHX-comp\sDNA\/Helitron$/DHH#RC\/Helitron/' | \
sed 's/DHX-comp\sLINE\/Jockey$/DHH#RC\/Helitron/' | \
sed 's/DHX-comp-chim$/DHH#RC\/Helitron_chim/' | \
sed 's/DHX-comp-chim\sDNA\/Helitron$/DHH#RC\/Helitron_chim/' | \
sed 's/DHX-incomp$/DHH#RC\/Helitron/' | \
sed 's/DHX-incomp\sDNA\/Helitron$/DHH#RC\/Helitron/' | \
sed 's/DHX-incomp\sDNA\/Mariner$/DHH#RC\/Helitron/' | \
sed 's/DHX-incomp\sDNA\/hAT$/DHH#RC\/Helitron/' | \
sed 's/DHX-incomp\sSINE\/SINE$/DHH#RC\/Helitron/' | \
sed 's/DHX-incomp-chim$/DHH#RC\/Helitron_chim/' | \
sed 's/DMX-incomp$/DMM#DNA\/Maverick/' | \
sed 's/DMX-incomp\sDNA\/Polinton$/DMM#DNA\/Maverick/' | \
sed 's/DMX-incomp-chim\sDNA\/Polinton$/DMM#DNA\/Maverick_chim/' | \
sed 's/DTX-comp$/DTX#DNA/' | \
sed 's/DTX-comp\sDNA\/Chapaev$/DTX#DNA\/Chapaev/' | \
sed 's/DTX-comp\sDNA\/Mariner$/DTT#DNA\/TcMar/' | \
sed 's/DTX-comp\sDNA\/MuDr$/DTX#DNA\/MuDr/' | \
sed 's/DTX-comp\sDNA\/P$/DTP#DNA\/P/' | \
sed 's/DTX-comp\sDNA\/PiggyBac$/DTB#DNA\/PiggyBac/' | \
sed 's/DTX-comp\sDNA\/hAT$/DTA#DNA\/hAT/' | \
sed 's/DTX-incomp$/DTX#DNA/' | \
sed 's/DTX-incomp\sDNA\/Academ$/DTX#DNA\/Academ/' | \
sed 's/DTX-incomp\sDNA\/Chapaev$/DTX#DNA\/Chapaev/' | \
sed 's/DTX-incomp\sDNA\/Ginger$/DTX#DNA\/Ginger/' | \
sed 's/DTX-incomp\sDNA\/Harbinger$/DTH#DNA\/Harbinger/' | \
sed 's/DTX-incomp\sDNA\/Helitron$/DTX#DNA/' | \
sed 's/DTX-incomp\sDNA\/Kolobok$/DTX#DNA\/Kolobok/' | \
sed 's/DTX-incomp\sDNA\/Mariner$/DTT#DNA\/TcMar/' | \
sed 's/DTX-incomp\sDNA\/MuDr$/DTX#DNA\/MuDr/' | \
sed 's/DTX-incomp\sDNA\/P$/DTP#DNA\/P/' | \
sed 's/DTX-incomp\sDNA\/PiggyBac$/DTB#DNA\/PiggyBac/' | \
sed 's/DTX-incomp\sDNA\/Polinton$/DMM#DNA\/Maverick/' | \
sed 's/DTX-incomp\sDNA\/Sola$/DTX#DNA\/Sola/' | \
sed 's/DTX-incomp\sDNA\/Transib$/DTR#DNA\/Transib/' | \
sed 's/DTX-incomp\sDNA\/hAT$/DTA#DNA\/hAT/' | \
sed 's/DTX-incomp\sLINE\/CR1$/DTX#DNA/' | \
sed 's/DTX-incomp\sLINE\/I$/DTX#DNA/' | \
sed 's/DTX-incomp\sLINE\/Jockey$/DTX#DNA/' | \
sed 's/DTX-incomp\sLINE\/L2$/DTX#DNA/' | \
sed 's/DTX-incomp\sLINE\/Penelope$/DTX#DNA/' | \
sed 's/DTX-incomp\sLINE\/R1$/DTX#DNA/' | \
sed 's/DTX-incomp\sLINE\/RTE$/DTX#DNA/' | \
sed 's/DTX-incomp\sLTR\/BEL$/DTX#DNA/' | \
sed 's/DTX-incomp\sLTR\/Copia$/DTX#DNA/' | \
sed 's/DTX-incomp\sLTR\/Gypsy$/DTX#DNA/' | \
sed 's/DTX-incomp\sOther\/EnSpm$/DTX#DNA/' | \
sed 's/DTX-incomp\sOther\/Homo$/DTX#DNA/' | \
sed 's/DTX-incomp-chim$/DTX#DNA\/DNA_chim/' | \
sed 's/DTX-incomp-chim\sDNA\/Chapaev$/DTX#DNA\/Chapaev_chim/' | \
sed 's/DTX-incomp-chim\sDNA\/Harbinger$/DTH#DNA\/Harbinger_chim/' | \
sed 's/DTX-incomp-chim\sDNA\/Helitron$/DTX#DNA\/DNA_chim/' | \
sed 's/DTX-incomp-chim\sDNA\/Mariner$/DTT#DNA\/TcMar_chim/' | \
sed 's/DTX-incomp-chim\sDNA\/PiggyBac$/DTB#DNA\/PiggyBac_chim/' | \
sed 's/DTX-incomp-chim\sDNA\/Polinton$/DMM#DNA\/Maverick_chim/' | \
sed 's/DTX-incomp-chim\sDNA\/hAT$/DTA#DNA\/hAT_chim/' | \
sed 's/DTX-incomp-chim\sLINE\/I$/DTX#DNA\/DNA_chim/' | \
sed 's/DTX-incomp-chim\sLINE\/Jockey$/DTX#DNA\/DNA_chim/' | \
sed 's/DTX-incomp-chim\sLINE\/RTE$/DTX#DNA\/DNA_chim/' | \
sed 's/DTX-incomp-chim\sLTR\/Gypsy$/DTX#DNA\/DNA_chim/' | \
sed 's/DXX$/DXX#DNA/' | \
sed 's/DXX\sDNA\/Academ$/DTX#DNA\/Academ/' | \
sed 's/DXX\sDNA\/Chapaev$/DTX#DNA\/Chapaev/' | \
sed 's/DXX\sDNA\/Ginger$/DTX#DNA\/Ginger/' | \
sed 's/DXX\sDNA\/Harbinger$/DTH#DNA\/Harbinger/' | \
sed 's/DXX\sDNA\/ISL2EU$/DTX#DNA\/ISL2EU/' | \
sed 's/DXX\sDNA\/Kolobok$/DTX#DNA\/Kolobok/' | \
sed 's/DXX\sDNA\/Mariner$/DTT#DNA\/TcMar/' | \
sed 's/DXX\sDNA\/Sola$/DTX#DNA\/Sola/' | \
sed 's/DXX\sDNA\/hAT$/DTA#DNA\/hAT/' | \
sed 's/DXX\sLINE\/Crack$/DXX#DNA/' | \
sed 's/DXX\sLINE\/RTE$/DXX#DNA/' | \
sed 's/DXX\sLTR\/Copia$/DXX#DNA/' | \
sed 's/DXX\sLTR\/Gypsy$/DXX#DNA/' | \
sed 's/DXX\sOther\/HERV$/DXX#DNA/' | \
sed 's/DXX\sSINE\/SINE$/DXX#DNA/' | \
sed 's/DXX-MITE$/DXX#DNA\/MITE/' | \
sed 's/DXX-MITE\sDNA\/Ginger$/DXX#DNA\/MITE/' | \
sed 's/DXX-MITE\sDNA\/Harbinger$/DXX#DNA\/MITE/' | \
sed 's/DXX-MITE\sDNA\/Helitron$/DXX#DNA\/MITE/' | \
sed 's/DXX-MITE\sDNA\/Mariner$/DXX#DNA\/MITE/' | \
sed 's/DXX-MITE\sDNA\/MuDr$/DXX#DNA\/MITE/' | \
sed 's/DXX-MITE\sDNA\/Polinton$/DXX#DNA\/MITE/' | \
sed 's/DXX-MITE\sDNA\/hAT$/DXX#DNA\/MITE/' | \
sed 's/DXX-MITE\sLINE\/Penelope$/DXX#DNA\/MITE/' | \
sed 's/DXX-MITE\sLTR\/Copia$/DXX#DNA\/MITE/' | \
sed 's/DXX-MITE\sOther\/EnSpm$/DXX#DNA\/MITE/' | \
sed 's/DXX-MITE\sSINE\/SINE$/DXX#DNA\/MITE/' | \
sed 's/DXX-MITE-chim$/DXX#DNA\/MITE_chim/' | \
sed 's/PotentialHostGene$/HG#HG\/HG/' | \
sed 's/PotentialHostGene\sDNA\/hAT$/HG#HG\/HG/' | \
sed 's/PotentialHostGene\sDNA\/Harbinger$/HG#HG\/HG/' | \
sed 's/PotentialHostGene\sDNA\/Helitron$/HG#HG\/HG/' | \
sed 's/PotentialHostGene\sDNA\/Mariner$/HG#HG\/HG/' | \
sed 's/PotentialHostGene\sDNA\/MuDr$/HG#HG\/HG/' | \
sed 's/PotentialHostGene\sDNA\/P$/HG#HG\/HG/' | \
sed 's/PotentialHostGene\sDNA\/PiggyBac$/HG#HG\/HG/' | \
sed 's/PotentialHostGene\sDNA\/Polinton$/HG#HG\/HG/' | \
sed 's/PotentialHostGene$/HG#HG\/HG/' | \
sed 's/PotentialHostGene\sLINE\/CR1$/HG#HG\/HG/' | \
sed 's/PotentialHostGene\sLINE\/Jockey$/HG#HG\/HG/' | \
sed 's/PotentialHostGene\sLINE\/RTE$/HG#HG\/HG/' | \
sed 's/PotentialHostGene\sLTR\/Copia$/HG#HG\/HG/' | \
sed 's/PotentialHostGene\sLTR\/Gypsy$/HG#HG\/HG/' | \
sed 's/PotentialHostGene\sSINE\/SINE$/HG#HG\/HG/' | \
sed 's/RIX-comp$/RIX#LINE/' | \
sed 's/RIX-comp\sLINE\/CR1$/RIX#LINE\/CR1/' | \
sed 's/RIX-comp\sLINE\/Crack$/RIX#LINE\/Crack/' | \
sed 's/RIX-comp\sLINE\/I$/RII#LINE\/I/' | \
sed 's/RIX-comp\sLINE\/Jockey$/RIJ#LINE\/Jockey/' | \
sed 's/RIX-comp\sLINE\/LOA$/RIX#LINE\/LOA/' | \
sed 's/RIX-comp\sLINE\/R1$/RIX#LINE\/R1/' | \
sed 's/RIX-comp\sLINE\/RTE$/RIR#LINE\/R2/' | \
sed 's/RIX-comp\sLINE\/Tx1$/RIX#LINE\/Tx1/' | \
sed 's/RIX-comp-chim\sLINE\/R1$/RIX#LINE\/R1/' | \
sed 's/RIX-comp-chim\sLINE\/RTE$/RIT#LINE\/RTE/' | \
sed 's/RIX-incomp$/RIX#LINE/' | \
sed 's/RIX-incomp\sDNA\/Academ$/RIX#LINE/' | \
sed 's/RIX-incomp\sDNA\/Chapaev$/RIX#LINE/' | \
sed 's/RIX-incomp\sDNA\/Ginger$/RIX#LINE/' | \
sed 's/RIX-incomp\sDNA\/Helitron$/RIX#LINE/' | \
sed 's/RIX-incomp\sDNA\/Kolobok$/RIX#LINE/' | \
sed 's/RIX-incomp\sDNA\/Polinton$/RIX#LINE/' | \
sed 's/RIX-incomp\sDNA\/Sola$/RIX#LINE/' | \
sed 's/RIX-incomp\sDNA\/hAT$/RIX#LINE/' | \
sed 's/RIX-incomp\sLINE\/CR1$/RIX#LINE\/CR1/' | \
sed 's/RIX-incomp\sLINE\/Crack$/RIX#LINE\/Crack/' | \
sed 's/RIX-incomp\sLINE\/I$/RII#LINE\/I/' | \
sed 's/RIX-incomp\sLINE\/Jockey$/RIJ#LINE\/Jockey/' | \
sed 's/RIX-incomp\sLINE\/L1$/RIL#LINE\/L1/' | \
sed 's/RIX-incomp\sLINE\/L2$/RIX#LINE\/L2/' | \
sed 's/RIX-incomp\sLINE\/LOA$/RIX#LINE\/LOA/' | \
sed 's/RIX-incomp\sLINE\/R1$/RIX#LINE\/R1/' | \
sed 's/RIX-incomp\sLINE\/RTE$/RIT#LINE\/RTE/' | \
sed 's/RIX-incomp\sLINE\/Tx1$/RIX#LINE\/Tx1/' | \
sed 's/RIX-incomp\sLTR\/Copia$/RIX#LINE/' | \
sed 's/RIX-incomp\sLTR\/Gypsy$/RIX#LINE/' | \
sed 's/RIX-incomp\sSINE\/SINE$/RIX#LINE/' | \
sed 's/RIX-incomp-chim$/RIX#LINE\/LINE_chim/' | \
sed 's/RIX-incomp-chim\sLINE\/CR1$/RIX#LINE\/CR1_chim/' | \
sed 's/RIX-incomp-chim\sLINE\/Jockey$/RIJ#LINE\/Jockey_chim/' | \
sed 's/RIX-incomp-chim\sLINE\/LOA$/RIX#LINE\/LOA_chim/' | \
sed 's/RIX-incomp-chim\sLINE\/RTE$/RIT#LINE\/RTE_chim/' | \
sed 's/RIX-incomp-chim\sLTR\/Gypsy$/RIX#LINE/' | \
sed 's/RLX-comp\sLTR\/Copia$/RLC#LTR\/Copia/' | \
sed 's/RLX-incomp$/RLX#LTR/' | \
sed 's/RLX-incomp\sDNA\/Helitron$/RLX#LTR/' | \
sed 's/RLX-incomp\sDNA\/Polinton$/RLX#LTR/' | \
sed 's/RLX-incomp\sDNA\/hAT$/RLX#LTR/' | \
sed 's/RLX-incomp\sLINE\/Jockey$/RLX#LTR/' | \
sed 's/RLX-incomp\sLINE\/RTE$/RLX#LTR/' | \
sed 's/RLX-incomp\sLTR\/BEL$/RLB#LTR\/BEL/' | \
sed 's/RLX-incomp\sLTR\/Copia$/RLC#LTR\/Copia/' | \
sed 's/RLX-incomp\sLTR\/Gypsy$/RLG#LTR\/Gypsy/' | \
sed 's/RLX-incomp\sOther\/EnSpm$/RLX#LTR/' | \
sed 's/RLX-incomp-chim$/RLX#LTR\/LTR_chim/' | \
sed 's/RLX-incomp-chim\sDNA\/Polinton$/RLX#LTR\/LTR_chim/' | \
sed 's/RLX-incomp-chim\sLTR\/BEL$/RLB#LTR\/BEL_chim/' | \
sed 's/RLX-incomp-chim\sLTR\/Gypsy$/RLG#LTR\/Gypsy_chim/' | \
sed 's/RPX-comp$/RPP#LTR\/Penelope/' | \
sed 's/RPX-comp\sLINE\/Penelope$/RPP#LTR\/Penelope/' | \
sed 's/RPX-comp\sLINE\/Poseidon$/RPP#LINE\/Poseidon/' | \
sed 's/RPX-comp-chim$/RPP#LTR\/Penelope_chim/' | \
sed 's/RPX-incomp$/RPP#LTR\/Penelope/' | \
sed 's/RPX-incomp\sLINE\/Penelope$/RPP#LTR\/Penelope/' | \
sed 's/RPX-incomp\sLINE\/Poseidon$/RPP#LINE\/Poseidon/' | \
sed 's/RPX-incomp-chim$/RPP#LTR\/Penelope_chim/' | \
sed 's/RPX-incomp-chim\sLINE\/I$/RPP#LTR\/Penelope/' | \
sed 's/RSX-comp$/RSX#SINE/' | \
sed 's/RSX-comp\sLTR\/Gypsy$/RSX#SINE/' | \
sed 's/RSX-comp\sSINE\/SINE$/RSX#SINE/' | \
sed 's/RSX-incomp$/RSX#SINE/' | \
sed 's/RSX-incomp\sDNA\/Harbinger$/RSX#SINE/' | \
sed 's/RSX-incomp\sDNA\/Kolobok$/RSX#SINE/' | \
sed 's/RSX-incomp\sDNA\/P$/RSX#SINE/' | \
sed 's/RSX-incomp\sDNA\/Polinton$/RSX#SINE/' | \
sed 's/RSX-incomp\sDNA\/Sola$/RSX#SINE/' | \
sed 's/RSX-incomp\sDNA\/hAT$/RSX#SINE/' | \
sed 's/RSX-incomp\sLINE\/CR1$/RSX#SINE/' | \
sed 's/RSX-incomp\sLINE\/Jockey$/RSX#SINE/' | \
sed 's/RSX-incomp\sLINE\/RTE$/RSX#SINE/' | \
sed 's/RSX-incomp\sLTR\/Copia$/RSX#SINE/' | \
sed 's/RSX-incomp\sLTR\/Gypsy$/RSX#SINE/' | \
sed 's/RSX-incomp\sSINE\/SINE$/RSX#SINE/' | \
sed 's/RSX-incomp-chim$/RSX#SINE\/SINE_chim/' | \
sed 's/RSX-incomp-chim\sDNA\/Mariner$/RSX#SINE\/SINE_chim/' | \
sed 's/RSX-incomp-chim\sDNA\/hAT$/RSX#SINE\/SINE_chim/' | \
sed 's/RSX-incomp-chim\sLINE\/Jockey$/RSX#SINE\/SINE_chim/' | \
sed 's/RSX-incomp-chim\sSINE\/SINE$/RSX#SINE\/SINE_chim/' | \
sed 's/RXX$/RXX#RETRO/' | \
sed 's/RXX\sDNA\/Polinton$/RXX#RETRO/' | \
sed 's/RXX\sDNA\/hAT$/RXX#RETRO/' | \
sed 's/RXX\sLINE\/CR1$/RIX#LINE\/CR1/' | \
sed 's/RXX\sLINE\/I$/RII#LINE\/I/' | \
sed 's/RXX\sLINE\/Jockey$/RIJ#LINE\/Jockey/' | \
sed 's/RXX\sLINE\/L2$/RIX#LINE\/L2/' | \
sed 's/RXX\sLINE\/R1$/RIX#LINE\/R1/' | \
sed 's/RXX\sLINE\/RTE$/RIT#LINE\/RTE/' | \
sed 's/RXX\sLINE\/Tx1$/RIX#LINE\/Tx1/' | \
sed 's/RXX\sLTR\/BEL$/RLB#LTR\/BEL/' | \
sed 's/RXX\sLTR\/Copia$/RLC#LTR\/Copia/' | \
sed 's/RXX\sLTR\/Gypsy$/RLG#LTR\/Gypsy/' | \
sed 's/RXX-LARD$/RXX#Other\/LARD/' | \
sed 's/RXX-TRIM$/RXX#Other\/TRIM/' | \
sed 's/RXX-chim$/RXX#RETRO/' | \
sed 's/RYX-comp$/RYX#DIRS\/Other/' | \
sed 's/RYX-comp-chim$/RYX#DIRS\/Other_chim/' | \
sed 's/RYX-incomp$/RYX#DIRS\/Other/' | \
sed 's/RYX-incomp\sLTR\/Gypsy$/RYX#DIRS\/Other/' | \
sed 's/RYX-incomp-chim$/RYX#DIRS\/Other_chim/' | \
sed 's/RYX-incomp-chim\sLTR\/Copia$/RYX#DIRS\/Other_chim/' | \
sed 's/RYX-incomp-chim\sLTR\/Gypsy$/RYX#DIRS\/Other_chim/' | \
sed 's/SSR$/SSR#Other/' | \
sed 's/SSR\sDNA\/Polinton$/SSR#Other/' | \
sed 's/SSR\sLINE\/I$/SSR#Other/' | \
sed 's/SSR\sLTR\/Copia$/SSR#Other/' | \
sed 's/SSR\sLTR\/Gypsy$/SSR#Other/' | \
sed 's/Unknown$/Unknown#Unknown/' | \
sed 's/noCat$/Unknown#Unknown/' | \
sed 's/noCat\sDNA\/Academ$/DTX#DNA\/Academ/' | \
sed 's/noCat\sDNA\/Chapaev$/DTX#DNA\/Chapaev/' | \
sed 's/noCat\sDNA\/Ginger$/DTX#DNA\/Ginger/' | \
sed 's/noCat\sDNA\/Harbinger$/DTX#DNA\/Harbinger/' | \
sed 's/noCat\sDNA\/Helitron$/DHH#RC\/Helitron/' | \
sed 's/noCat\sDNA\/ISL2EU$/DTX#DNA\/ISL2EU/' | \
sed 's/noCat\sDNA\/Kolobok$/DTX#DNA\/Kolobok/' | \
sed 's/noCat\sDNA\/Mariner$/DTT#DNA\/TcMar/' | \
sed 's/noCat\sDNA\/MuDr$/DTX#DNA\/MuDr/' | \
sed 's/noCat\sDNA\/P$/DTP#DNA\/P/' | \
sed 's/noCat\sDNA\/PiggyBac$/DTB#DNA\/PiggyBac/' | \
sed 's/noCat\sDNA\/Polinton$/DMM#DNA\/Maverick/' | \
sed 's/noCat\sDNA\/Sola$/DTX#DNA\/Sola/' | \
sed 's/noCat\sDNA\/hAT$/DTA#DNA\/hAT/' | \
sed 's/noCat\sLINE\/CR1$/RIX#LINE\/CR1/' | \
sed 's/noCat\sLINE\/Crack$/RIX#LINE\/Crack/' | \
sed 's/noCat\sLINE\/I$/RII#LINE\/I/' | \
sed 's/noCat\sLINE\/Jockey$/RIJ#LINE\/Jockey/' | \
sed 's/noCat\sLINE\/L1$/RIL#LINE\/L1/' | \
sed 's/noCat\sLINE\/L2$/RIX#LINE\/L2/' | \
sed 's/noCat\sLINE\/Penelope$/RIX#LTR\/Penelope/' | \
sed 's/noCat\sLINE\/Poseidon$/RIX#LINE\/Poseidon/' | \
sed 's/noCat\sLINE\/R1$/RIX#LINE\/R1/' | \
sed 's/noCat\sLINE\/RTE$/RIT#LINE\/RTE/' | \
sed 's/noCat\sLINE\/Tx1$/RIX#LINE\/Tx1/' | \
sed 's/noCat\sLTR\/BEL$/RLB#LTR\/BEL/' | \
sed 's/noCat\sLTR\/Copia$/RLC#LTR\/Copia/' | \
sed 's/noCat\sLTR\/Gypsy$/RLG#LTR\/Gypsy/' | \
sed 's/noCat\sOther\/EnSpm$/XXX#Other/' | \
sed 's/noCat\sOther\/HERV$/XXX#Other/' | \
sed 's/noCat\sOther\/Homo$/XXX#Other/' | \
sed 's/noCat\sSINE\/SINE$/XXX#SINE/' \
> $2
