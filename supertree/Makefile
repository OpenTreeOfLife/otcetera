STEP_1_DIR=step_1
STEP_2_DIR=step_2
STEP_3_DIR=step_3
STEP_4_DIR=step_4
STEP_5_DIR=step_5
STEP_6_DIR=step_6

# $(STEP_3_DIR)/contested-ott-ids.txt is not (currently) a part of the pipeline
# so we go from step2 to step4
#	$(STEP_3_DIR)/contested-ott-ids.txt

ARTIFACTS=$(STEP_4_DIR)/step-5-phylo-args.txt \
	$(STEP_5_DIR)/taxonomy.tre \
	$(STEP_6_DIR)/taxonomy.tre

# default is "all"
all: $(ARTIFACTS)
	

STEP_5_INPUT_FILENAMES=$(shell cat $(STEP_1_DIR)/tree-ranking.txt)
STEP_5_PHYLO_PATHS=$(addprefix $(STEP_4_DIR)/, $(STEP_5_INPUT_FILENAMES))

$(STEP_4_DIR)/step-5-phylo-args.txt : $(STEP_1_DIR)/taxonomy.tre $(STEP_1_DIR)/tree-ranking.txt
	rm -f $(STEP_4_DIR)/step-5-phylo-args.txt
	for f in $(STEP_5_PHYLO_PATHS) ; do echo $$f >> $(STEP_4_DIR)/step-5-phylo-args.txt ; done

$(STEP_5_DIR)/taxonomy.tre : $(STEP_5_PHYLO_PATHS) $(STEP_1_DIR)/taxonomy.tre $(STEP_4_DIR)/step-5-phylo-args.txt
	otc-nonterminals-to-exemplars -e$(STEP_5_DIR) $(STEP_1_DIR)/taxonomy.tre -f$(STEP_4_DIR)/step-5-phylo-args.txt

$(STEP_6_DIR)/taxonomy.tre : $(STEP_5_DIR)/taxonomy.tre
	cp $(STEP_5_DIR)/taxonomy.tre $(STEP_6_DIR)/taxonomy.tre

clean:
	rm -f $(ARTIFACTS)



