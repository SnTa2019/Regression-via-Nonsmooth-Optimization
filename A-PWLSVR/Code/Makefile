# Makefile for the A-PWLSVR

FF = gfortran
FFLAGS = -fbounds-check
#FFLAGS = -fbounds-check -Wall
OPEN = 
#OPEN =  -fopenmp 
#RM = del   #in windows
RM = rm    #in unix/linux/mac

all: spr

spr: constants.o initspr.o bundle1.o bundle2.o functions.o dbdc.o spr.o plqdf1.o 
	$(FF) -o spr $(FFLAGS) $(OPEN) constants.o initspr.o bundle1.o bundle2.o functions.o dbdc.o spr.o plqdf1.o 

constants.mod: constants.o constants.f03
	$(FF) -c $(FFLAGS) $(OPEN) constants.f03
	
constants.o: constants.f03
	$(FF) -c $(FFLAGS) $(OPEN) constants.f03

initspr.o: constants.mod initspr.f03
	$(FF) -c $(FFLAGS) $(OPEN) initspr.f03

bundle1.mod: constants.mod bundle1.o bundle1.f03 
	$(FF) -c $(FFLAGS) $(OPEN) bundle1.f03
	
bundle1.o: constants.mod bundle1.f03
	$(FF) -c $(FFLAGS) $(OPEN) bundle1.f03 

bundle2.mod: constants.mod bundle2.o bundle2.f03
	$(FF) -c $(FFLAGS) $(OPEN) bundle2.f03 
	
bundle2.o: constants.mod bundle2.f03
	$(FF) -c $(FFLAGS) $(OPEN) bundle2.f03 

functions.mod: constants.mod functions.o functions.f03
	$(FF) -c $(FFLAGS) $(OPEN) functions.f03 
	
functions.o: constants.mod functions.f03
	$(FF) -c $(FFLAGS) $(OPEN) functions.f03 

dbdc.mod: constants.mod bundle1.mod bundle2.mod functions.mod dbdc.o dbdc.f03 
	$(FF) -c $(FFLAGS) $(OPEN) dbdc.f03	 
	
dbdc.o: constants.mod bundle1.mod bundle2.mod functions.mod dbdc.f03
	$(FF) -c $(FFLAGS) $(OPEN) dbdc.f03 
	
spr.o: constants.mod bundle1.mod bundle2.mod functions.mod dbdc.mod spr.f03
	$(FF) -c $(FFLAGS) $(OPEN) spr.f03 

plqdf1.o: plqdf1.f
	$(FF) -c $(FFLAGS) $(OPEN) plqdf1.f 

clean:	
	$(RM) spr *.mod *.o  
	echo Clean done