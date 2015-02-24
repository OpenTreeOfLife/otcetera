STEP_1_DIR=step-1-pruned-overlapping-tips
STEP_1_INPUT_FILENAMES=$(shell cat $(STEP_1_DIR)/tree-ranking.txt)
STEP_1_INPUTS=$(addprefix $(STEP_1_DIR)/, $(STEP_1_INPUT_FILENAMES))

STEP_2_DIR=step-2-pruned-taxonomy
$(STEP_2_DIR)/taxonomy.tre : input/taxonomy.tre $(STEP_1_INPUTS)
	cd $(STEP_1_DIR) ; pwd; \
	if otcpruntaxonomy ../input/taxonomy -ftree-ranking.txt > ../.output ; \
	then mv ../.output $(STEP_2_DIR)/taxonomy.tre ; \
	fi
