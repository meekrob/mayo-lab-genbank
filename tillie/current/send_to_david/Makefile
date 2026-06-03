TIMESTAMP := $(shell date "+%Y-%m-%d-%H-%M-%S")
PREFIX := tillie_submission_files-$(TIMESTAMP)
default:
	$(MAKE) -C A2_submission_files 
	$(MAKE) -C A3_submission_files 
	$(MAKE) -C FABADRU_SA_submission_files
	mkdir $(PREFIX)
	cp -R A2_submission_files A3_submission_files FABADRU_SA_submission_files $(PREFIX)
	tar -cvzf $(PREFIX).tgz --exclude "*.sbt" --exclude table2asn $(PREFIX)
	rm -Rf $(PREFIX)
    
