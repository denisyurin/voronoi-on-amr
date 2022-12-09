BINNAME	:= drect-log-grid
SRCNAME	:= rect-log-grid
FLAGS := -w -O3 -pthread

define COMPILE
  g++ -I. $(SRCNAME).cpp -o $(BINNAME)  $(FLAGS)
endef

RUN := ./$(BINNAME)

standard: 
	clear
	$(COMPILE)
	$(RUN)

bg: 
	clear
	$(COMPILE)
	$(RUN) 1> log.out.txt 2> log.err.txt &
