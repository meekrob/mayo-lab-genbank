default: additional_oc_sequences_for_submission.fsa additional_oc_sequences_for_submission.tbl #additional_oc_sequences_for_submission.sqn


additional_oc_sequences_for_submission.fsa: TD_n1257_fixed_strain.fsa
	ln -v $^ $@

TD_n1257_fixed_strain.fsa: 
	python map_tsv_to_fasta_header.py new_files/GenBank_Sequences_A2_A3_CLEAN_fixed_strain.txt new_files/TD_n1257_GenBank.fasta > $@

additional_oc_sequences_for_submission.tbl:
	gff_to_feature_table <(grep -v 'Editing History' new_files/TD_n1257_GenBank.gff) > $@

additional_oc_sequences_for_submission.sqn: additional_oc_sequences_for_submission.fsa additional_oc_sequences_for_submission.tbl
	table2asn -indir . -a s -V v -t A2_template.sbt
