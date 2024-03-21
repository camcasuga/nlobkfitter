AMPLITUDEINCLUDE=$(NLOBKFITTER_LIB_PATH)/amplitudelib
AMPLITUDEPATH=$(NLOBKFITTER_LIB_PATH)/amplitudelib/build/lib
CUBAPATH=$(NLOBKFITTER_LIB_PATH)/Cuba-4.2.2/lib/
MINUITINC=$(NLOBKFITTER_LIB_PATH)/Minuit2/include/
MINUITLIBDIR=$(NLOBKFITTER_LIB_PATH)/Minuit2/lib/
CXXFLAGS = `gsl-config --cflags` -std=c++11 -O3 -fopenmp -pedantic -I$(AMPLITUDEINCLUDE) -I src/ -I$(MINUITINC)
LDFLAGS = `gsl-config --libs` -L$(MINUITLIBDIR) -lMinuit2 -L$(CUBAPATH) -lcuba -lm

include filelist.m

all: fitex swarmscan lightqscan

%.o: %.cpp $(DEPENDENCIES) # Works but recompiles all objects if a single header changes. Safe.
	g++ -c -o $@ $< $(CXXFLAGS)

#$(OBJECTS): %.o: %.cpp %.hpp # Recompiles only the corresponding object after a header change. Uncertain if this is enough, possibly unsafe.
#	g++ -c -o $@ $< $(CXXFLAGS)

fit: $(OBJECTS) src/main.o
	g++ $(OBJECTS) src/main.o amplitudelib_v2/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o fit $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

fitex: $(OBJECTS) src/main_switches.o 
	g++ $(OBJECTS) src/main_switches.o $(AMPLITUDEPATH)/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o fitex $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

lightqfitex: $(OBJECTS) src/lightq-fit.o 
	g++ $(OBJECTS) src/lightq-fit.o $(AMPLITUDEPATH)/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o lightqfitex $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)


swarmscan: $(OBJECTS) src/swarmscan.o src/nlodissigmar.hpp
	g++ $(OBJECTS) src/swarmscan.o $(AMPLITUDEPATH)/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o swarmscan $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

massivescan: $(OBJECTS) src/massive_swarmscan.o src/nlodissigmar.hpp
	g++ $(OBJECTS) src/massive_swarmscan.o $(AMPLITUDEPATH)/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o massivescan $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

lightqscan: $(OBJECTS) src/lightq_scan.o
	g++ $(OBJECTS) src/lightq_scan.o $(AMPLITUDEPATH)/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o lightqscan $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

nuclearstrufun: $(OBJECTS) src/unint_nucl_structfun.o 
	g++ $(OBJECTS) src/unint_nucl_structfun.o $(AMPLITUDEPATH)/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o nuclearstrufun $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

plottool: $(OBJECTS) src/tool_plot.o 
	g++ $(OBJECTS) src/tool_plot.o $(AMPLITUDEPATH)/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o plottool $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

fltool: $(OBJECTS) src/tool_fl.o
	g++ $(OBJECTS) src/tool_fl.o $(AMPLITUDEPATH)/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o fltool $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

lofit: $(OBJECTS) src/lo-fit-plot.o 
	g++ $(OBJECTS) src/lo-fit-plot.o $(AMPLITUDEPATH)/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o lofit $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

tool: $(OBJECTS) src/toolF2.o 
	g++ $(OBJECTS) src/toolF2.o $(AMPLITUDEPATH)/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o tool $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

subunsubtest: $(OBJECTS) src/sub_unsub_test.o 
	g++ $(OBJECTS) src/sub_unsub_test.o $(AMPLITUDEPATH)/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o subunsub $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

clean:
	rm -f $(OBJECTS) $(AMPLITUDELIBO) src/main.o src/toolF2.o src/main_switches.o src/tool_plot.o src/tool_fl.o src/swarmscan.o src/massive_swarmscan.o src/lightq_scan.o src/lightq-fit.o src/sub_unsub_test.o src/lo-fit-plot.o
	rm -f fit fitex tool plottool fltool swarmscan massivescan lightqscan subunsub lightqfitex lofit
