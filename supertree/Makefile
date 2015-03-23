STEP_1_DIR=step-1-pruned-overlapping-tips
STEP_2_DIR=step-2-pruned-taxonomy
STEP_3_DIR=step-3-contested-tips
STEP_4_DIR=step-4-subproblems

# $(STEP_3_DIR)/contested-ott-ids.txt is not (currently) a part of the pipeline
# so we go from step2 to step4
#	$(STEP_3_DIR)/contested-ott-ids.txt

ARTIFACTS=$(STEP_2_DIR)/taxonomy.tre \
	$(STEP_4_DIR)/subproblem-ids.txt

# default is "all"
all: $(ARTIFACTS)
	


STEP_1_INPUT_FILENAMES=$(shell cat $(STEP_1_DIR)/tree-ranking.txt)
STEP_1_INPUTS=$(addprefix $(STEP_1_DIR)/, $(STEP_1_INPUT_FILENAMES))

$(STEP_2_DIR)/taxonomy.tre : input/taxonomy.tre $(STEP_1_INPUTS)
	cd $(STEP_1_DIR) ; pwd; \
	if otcprunetaxonomy ../input/taxonomy.tre -ftree-ranking.txt > .output ; \
	then mv .output ../$(STEP_2_DIR)/taxonomy.tre ; \
	fi

$(STEP_3_DIR)/contested-ott-ids.txt : $(STEP_2_DIR)/taxonomy.tre $(STEP_1_INPUTS)
	cd $(STEP_1_DIR) ; pwd; \
	if otcdetectcontested ../$(STEP_2_DIR)/taxonomy.tre -ftree-ranking.txt >.output 2>../$(STEP_3_DIR)/contested-errstream.txt ; \
	then mv .output ../$(STEP_3_DIR)/contested-ott-ids.txt ; \
	fi

$(STEP_4_DIR)/subproblem-ids.txt :$(STEP_2_DIR)/taxonomy.tre $(STEP_1_INPUTS) $(STEP_1_DIR)/tree-ranking.txt
	cd $(STEP_1_DIR) ; pwd; \
	bash ../run-subproblem-finder.sh ../$(STEP_4_DIR) ../$(STEP_2_DIR)/taxonomy.tre -ftree-ranking.txt
clean:
	rm $(STEP_2_DIR)/taxonomy.tre
	rm $(STEP_3_DIR)/contested-ott-ids.txt



