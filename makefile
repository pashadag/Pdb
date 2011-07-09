SYSTEM     = x86-64_RHEL3.0_3.2
LIBFORMAT  = static_pic


COPT =  -g -Wall -Wno-sign-compare -fPIC -fexceptions 

# Create a list of the object files using macro substitutions
OBJECTS = $(SOURCES:.cpp=.o)

# redefine inbuilt macros and define new ones
CC    = g++
CFLAGS = $(COPT) 


#.SUFFIXES:       # remove inbuilt definitions
#.SUFFIXES: .cpp .o # define default compilation procedure
#.cpp.o:            # .o files depend on .c files
#	$(CC) $(CFLAGS) $*.c # start with a tab character!!

# Default target
all: pabNL buildsa qcheck2 filter_genome ab pab generateReads singleGenerateReads pabNLpar pabFull short abshort
#	@echo "compilation done"

validate: validate.o
	$(CC) validate.o -o validate -lpthread 
validate.o : validate.cpp
	$(CC) -c $(COPT) validate.cpp -o validate.o
pabFull: pabFull.o
	$(CC) pabFull.o -o pabFull -lpthread 
pabFull.o : pabFull.cpp
	$(CC) -c $(COPT) pabFull.cpp -o pabFull.o
pabNLpar: pabNLpar.o
	$(CC) pabNLpar.o -o pabNLpar -lpthread 
pabNLpar.o : pabNLpar.cpp
	$(CC) -c $(COPT) pabNLpar.cpp -o pabNLpar.o
pabNL: pabNL.o
	$(CC) pabNL.o -o pabNL
pabNL.o : pabNL.cpp
	$(CC) -c $(COPT) pabNL.cpp -o pabNL.o
buildsa: buildsa.o
	$(CC) buildsa.o -o buildsa
buildsa.o : buildsa.cpp
	$(CC) -c $(COPT) buildsa.cpp -o buildsa.o
filter_genome: filter_genome.o
	$(CC) filter_genome.o -o filter_genome
filter_genome.o : filter_genome.cpp
	$(CC) -c $(COPT) filter_genome.cpp -o filter_genome.o
qcheck2: qcheck2.o
	$(CC) qcheck2.o -o qcheck2
qcheck2.o : qcheck2.cpp
	$(CC) -c $(COPT) qcheck2.cpp -o qcheck2.o
abshort: abshort.o
	$(CC) abshort.o -o abshort
abshort.o : abshort.cpp
	$(CC) -c $(COPT) abshort.cpp -o abshort.o
ab: ab.o
	$(CC) ab.o -o ab
ab.o : ab.cpp
	$(CC) -c $(COPT) ab.cpp -o ab.o
pab: pab.o
	$(CC) pab.o -o pab
pab.o : pab.cpp
	$(CC) -c $(COPT) pab.cpp -o pab.o
generateReads: generateReads.o
	$(CC) generateReads.o -o generateReads
generateReads.o : generateReads.cpp
	$(CC) -c $(COPT) generateReads.cpp -o generateReads.o
singleGenerateReads: singleGenerateReads.o
	$(CC) singleGenerateReads.o -o singleGenerateReads
singleGenerateReads.o : singleGenerateReads.cpp
	$(CC) -c $(COPT) singleGenerateReads.cpp -o singleGenerateReads.o
short: short.o
	$(CC) short.o -o short 
short.o : short.cpp
	$(CC) -c $(COPT) short.cpp -o short.o





# Target deleting unwanted files
clean:
	rm -f *.o *~ core mppcore
