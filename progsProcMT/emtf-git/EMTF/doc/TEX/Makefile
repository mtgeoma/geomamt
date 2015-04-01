TARGETS:	clean install EMTF_doc matlab_doc Z_files MMT_doc

EMTF_doc:	EMTF_doc.tex
		latex EMTF_doc
		dvips EMTF_doc	

MMT_doc:	MMT_doc.tex
		latex MMT_doc
		dvips MMT_doc	

matlab_doc:	matlab_doc.tex
		latex matlab_doc
		dvips matlab_doc	

Z_files:	Z_files.tex
		latex Z_files
		dvips Z_files	

clean:
	rm -f *.log *.aux *.bak *dvi

install:
	mv *.ps ../PS
