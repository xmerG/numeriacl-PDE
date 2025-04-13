all:
	$(MAKE) -C src

cleanoutput:
	rm -f output/*.json output/*.csv

report:
	$(MAKE) -C doc  

clean:
	$(MAKE) -C src clean
	$(MAKE) -C doc clean

.PHONY: all cleanoutput draw report clean
