SHELL := /bin/bash

default: additional_oc_sequences_for_submission.fsa additional_oc_sequences_for_submission.tbl additional_oc_sequences_for_submission.A2.sqn additional_oc_sequences_for_submission.A3.sqn

additional_oc_sequences_for_submission.fsa additional_oc_sequences_for_submission.gff:
	python map_tsv_to_fasta_header.py new_files/GenBank_Sequences_A2_A3_CLEAN_fixed_strain.txt new_files/TD_n1257_GenBank.fasta  new_files/TD_n1257_GenBank.gff

additional_oc_sequences_for_submission.tbl: additional_oc_sequences_for_submission.gff
	gff_to_feature_table <(grep -v 'Name=$$' $^) > $@

additional_oc_sequences_for_submission.A2.sqn: additional_oc_sequences_for_submission.fsa additional_oc_sequences_for_submission.tbl
	table2asn -indir . -a s -V v -t A2_template.sbt
	mv additional_oc_sequences_for_submission.sqn additional_oc_sequences_for_submission.A2.sqn

additional_oc_sequences_for_submission.A3.sqn: additional_oc_sequences_for_submission.fsa additional_oc_sequences_for_submission.tbl
	table2asn -indir . -a s -V v -t A3_template.sbt
	mv additional_oc_sequences_for_submission.sqn additional_oc_sequences_for_submission.A3.sqn
