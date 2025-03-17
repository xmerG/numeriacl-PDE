all:
	$(MAKE) -C src

cleanoutput:
	rm -f output/*.json

draw:
	python3 draw/plotA.py 

report:
	$(MAKE) -C doc  

clean:
	$(MAKE) -C src clean
	$(MAKE) -C doc clean
	rm -f output/*.json  

.PHONY: all cleanoutput draw report clean
