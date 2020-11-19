# makefile for all variables needed in more than one make

# VARIABLES
# all species codes
SPECIES = 1_Tdi 1_Tps 2_Tcm 2_Tsi 3_Tce 3_Tms 4_Tbi 4_Tte 5_Tge 5_Tpa


# RECEPIES

# sanity check (read variable)
## print-VARIABLENAME print variable value
print-% :
	@echo $* = $($*)

.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<
