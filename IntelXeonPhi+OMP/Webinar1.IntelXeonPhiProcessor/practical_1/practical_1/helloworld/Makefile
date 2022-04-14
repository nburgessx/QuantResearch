PANDOC := pandoc
PANDOC_FLAGS := -s -f markdown -t html

instructions.html : INSTRUCTIONS.md
	$(PANDOC) $(PANDOC_FLAGS) -o instructions.html INSTRUCTIONS.md






.PHONY : clean 



clean :
	rm -rf ./*.html
