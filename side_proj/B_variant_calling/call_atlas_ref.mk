ifndef WINDOW_SIZE
WINDOW_SIZE := 100000
endif

THETACALLS := $(patsubst %, data/%/variant_calls/ref/atlas/ref_to_b3v05_w$(WINDOW_SIZE)_theta_estimates.txt, $(SPECIES))

.PHONY : call.all
call.all : $(THETACALLS)

data/%/variant_calls/ref/atlas/ref_to_b3v05_w$(WINDOW_SIZE)_theta_estimates.txt :
	cd data/$*/variant_calls/logs && \
		$$VROOT/B_variant_calling/atlas_estTheta_ref.sh $* b3v05 $(WINDOW_SIZE)
