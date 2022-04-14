PANDOC := pandoc
PANDOC_FLAGS := -s -f markdown -t html





all : helloworld helloworld_par vectoradd_par vectoradd_parvec vectoradd_vec instructions.html


instructions.html : INSTRUCTIONS.md
	$(PANDOC) $(PANDOC_FLAGS) -o instructions.html INSTRUCTIONS.md



helloworld : 
	make -C helloworld


helloworld_par :
	make -C helloworld_par



vectoradd_par : 
	make -f Makefile.docs -C vectoradd_par


vectoradd_parvec :
	make -f Makefile.docs -C vectoradd_parvec


vectoradd_vec :
	make -f Makefile.docs -C vectoradd_vec





.PHONY : all clean helloworld helloworld_par vectoradd_par vectoradd_parvec vectoradd_vec



clean :
	make clean -C helloworld
	make clean -C helloworld_par
	make clean -C vectoradd_par
	make clean -C vectoradd_parvec
	make clean -C vectoradd_vec
	rm -rf ./instructions.html
