SHELL := /bin/bash

default: A3.sqn

working_submission_files.tgz: additional_oc_sequences_for_submission.txt additional_oc_sequences_for_submission.fsa additional_oc_sequences_for_submission.tbl additional_oc_sequences_for_submission.A2.sqn additional_oc_sequences_for_submission.A3.sqn
	cp $^ working
	tar -czvf $@ working

additional_oc_sequences_for_submission.fsa additional_oc_sequences_for_submission.gff additional_oc_sequences_for_submission.tbl:
	python map_tsv_to_fasta_header.py new_files/GenBank_Sequences_A2_A3_CLEAN_fixed_strain.txt new_files/TD_n1257_GenBank.fasta  new_files/TD_n1257_GenBank.gff

#additional_oc_sequences_for_submission.tbl: additional_oc_sequences_for_submission.gff
#	gff_to_feature_table <(grep -v 'Name=$$' $^) > $@

A3.sqn: A3.fsa A3.tbl
	./table2asn -indir . -a s -V vbt -t A3_template.sbt
	cat *.stats

additional_oc_sequences_for_submission.A3.sqn: additional_oc_sequences_for_submission.fsa additional_oc_sequences_for_submission.tbl
	table2asn -indir . -a s -V vb -t A3_template.sbt
	mv tillie.stats tillie.A3.stats
	mv additional_oc_sequences_for_submission.gbf additional_oc_sequences_for_submission.A3.gbf
	mv additional_oc_sequences_for_submission.sqn additional_oc_sequences_for_submission.A3.sqn
	mv additional_oc_sequences_for_submission.val additional_oc_sequences_for_submission.A3.val 
	cat tillie.A3.stats
