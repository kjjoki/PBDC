# Makefile for the proximal bundle method PBDC

FF = f95 
FFLAGS = 
OPEN = -fopenmp
# -pg, -fopenmp 

all: tpbdc

tpbdc: constants.o bundle1.o bundle2.o functions.o norm_min.o pbdc.o tpbdc.o fun.o plqdf1.o pvmm.o pqsubs.o mqsubs.o
	$(FF) -o tpbdc $(FFLAGS) $(OPEN) constants.o bundle1.o bundle2.o functions.o norm_min.o pbdc.o tpbdc.o fun.o plqdf1.o pvmm.o pqsubs.o mqsubs.o 

constants.mod: constants.o constants.f95
	$(FF) -c $(FFLAGS) $(OPEN) constants.f95
	
constants.o: constants.f95
	$(FF) -c $(FFLAGS) $(OPEN) constants.f95

bundle1.mod: constants.mod bundle1.o bundle1.f95 
	$(FF) -c $(FFLAGS) $(OPEN) bundle1.f95
	
bundle1.o: constants.mod bundle1.f95
	$(FF) -c $(FFLAGS) $(OPEN) bundle1.f95 

bundle2.mod: constants.mod bundle2.o bundle2.f95
	$(FF) -c $(FFLAGS) $(OPEN) bundle2.f95 
	
bundle2.o: constants.mod bundle2.f95
	$(FF) -c $(FFLAGS) $(OPEN) bundle2.f95 

functions.mod: constants.mod functions.o functions.f95
	$(FF) -c $(FFLAGS) $(OPEN) functions.f95 
	
functions.o: constants.mod functions.f95
	$(FF) -c $(FFLAGS) $(OPEN) functions.f95 

norm_min.mod: constants.mod bundle1.mod bundle2.mod functions.mod norm_min.o norm_min.f95
	$(FF) -c $(FFLAGS) $(OPEN) norm_min.f95 
	
norm_min.o: constants.mod bundle1.mod bundle2.mod functions.mod norm_min.f95
	$(FF) -c $(FFLAGS) $(OPEN) norm_min.f95 

pbdc.mod: constants.mod bundle1.mod bundle2.mod functions.mod norm_min.mod pbdc.o pbdc.f95 
	$(FF) -c $(FFLAGS) $(OPEN) pbdc.f95	 
	
pbdc.o: constants.mod bundle1.mod bundle2.mod functions.mod norm_min.mod pbdc.f95
	$(FF) -c $(FFLAGS) $(OPEN) pbdc.f95 
	
tpbdc.o: constants.mod bundle1.mod bundle2.mod functions.mod norm_min.mod pbdc.mod tpbdc.f95
	$(FF) -c $(FFLAGS) $(OPEN) tpbdc.f95 

fun.o: constants.mod norm_min.mod fun.f95
	$(FF) -c $(FFLAGS) $(OPEN) fun.f95 

plqdf1.o: plqdf1.f
	$(FF) -c $(FFLAGS) $(OPEN) plqdf1.f 

pqsubs.o: pqsubs.f
	$(FF) -c $(FFLAGS) $(OPEN) pqsubs.f 

mqsubs.o: mqsubs.f
	$(FF) -c $(FFLAGS) $(OPEN) mqsubs.f 

pvmm.o: pvmm.f
	$(FF) -c $(FFLAGS) $(OPEN) pvmm.f 	

clean:	
	rm tpbdc constants.mod constants.o bundle1.mod bundle1.o bundle2.mod bundle2.o functions.mod functions.o norm_min.mod norm_min.o pbdc.mod pbdc.o tpbdc.o fun.o plqdf1.o pvmm.o pqsubs.o mqsubs.o 
	echo Clean done