prefix = /home/baonguye/guinea-pig/install
exec_prefix = ${prefix}
bindir = ${exec_prefix}/bin
libdir = ${exec_prefix}/lib

CC = gcc -std=gnu11
CXX = g++
CXXFLAGS = 
CPPFLAGS =  -DUSE_NON_CONST
FFTW_CFLAGS = -DUSE_FFT_LOCAL
ROOT_CFLAGS = 
FFTW_LDFLAGS = 
ROOT_LDFLAGS = 
CFLAGS = $(CPPFLAGS) $(FFTW_CFLAGS) $(ROOT_CFLAGS) -O3 -ffast-math -fno-keep-inline-functions -Wall -Wundef -fPIC -m64 -DCOMPUTER64b
LDFLAGS =   

MODULES=besselCPP.o beamParametersCPP.o beamCPP.o backgroundCPP.o pairsCPP.o minijetsCPP.o \
	mathematicalTools.o lumiCPP.o jetParameterCPP.o guineapigCPP.o grvCPP.o gridCPP.o fourierCPP.o \
	fieldCPP.o bhabhaSamplesCPP.o timerCPP.o switchesCPP.o splineCPP.o rndmCPP.o resultsCPP.o principal.o \
	physicalTools.o particlesCPP.o particleBeamCPP.o readData.o memory.o tridentCPP.o option_args.o

all: guinea guinearoot

guinea: $(MODULES) 
	$(CXX) $(CFLAGS) $(MODULES) $(LDFLAGS) $(FFTW_LDFLAGS) -o guinea

guinearoot: $(MODULES) 
	$(CXX) $(CFLAGS) $(MODULES) $(LDFLAGS) $(FFTW_LDFLAGS) $(ROOT_LDFLAGS) -o guinearoot

.cc.o: 
	$(CXX) $(CFLAGS) $(CXXFLAGS) $< -c -o $@

.c.o:
	$(CC) $(CFLAGS) $< -c -o $@


dep:
	g++ -M -MM -E -MD *.cc


besselCPP.o: besselCPP.cc besselCPP.h splineCPP.h

beamParametersCPP.o: beamParametersCPP.cc beamParametersCPP.h \
  abstractParticle.h typeDefs.h mathematicalEntities.h \
  mathematicalTools.h besselCPP.h splineCPP.h rndmCPP.h define.h \
  physconst.h abstractIOclass.h parametersCPP.h readData.h memory.h

beamCPP.o: beamCPP.cc beamCPP.h beamParametersCPP.h abstractParticle.h \
  typeDefs.h mathematicalEntities.h mathematicalTools.h besselCPP.h \
  splineCPP.h rndmCPP.h define.h physconst.h abstractIOclass.h \
  parametersCPP.h readData.h memory.h particleBeamCPP.h \
  particlesCPP.h physicalTools.h mathconst.h fileInputOutput.h \
  fileInputOutputAscii.h IfileInputOutput.h

backgroundCPP.o: backgroundCPP.cc backgroundCPP.h typeDefs.h \
  fileInputOutput.h abstractParticle.h mathematicalEntities.h \
  mathematicalTools.h besselCPP.h splineCPP.h rndmCPP.h define.h \
  physconst.h abstractIOclass.h fileInputOutputAscii.h IfileInputOutput.h \
  meshCPP.h pairsCPP.h resultsCPP.h mathconst.h beamCPP.h \
  beamParametersCPP.h parametersCPP.h readData.h memory.h \
  particleBeamCPP.h particlesCPP.h physicalTools.h switchesCPP.h \
  jetParameterCPP.h minijetsCPP.h grvCPP.h

pairsCPP.o: pairsCPP.cc pairsCPP.h resultsCPP.h define.h physconst.h \
  mathconst.h beamCPP.h beamParametersCPP.h abstractParticle.h typeDefs.h \
  mathematicalEntities.h mathematicalTools.h besselCPP.h splineCPP.h \
  rndmCPP.h abstractIOclass.h parametersCPP.h readData.h memory.h \
  particleBeamCPP.h particlesCPP.h physicalTools.h fileInputOutput.h \
  fileInputOutputAscii.h IfileInputOutput.h switchesCPP.h meshCPP.h

minijetsCPP.o: minijetsCPP.cc minijetsCPP.h typeDefs.h fileInputOutput.h \
  abstractParticle.h mathematicalEntities.h mathematicalTools.h \
  besselCPP.h splineCPP.h rndmCPP.h define.h physconst.h \
  abstractIOclass.h fileInputOutputAscii.h IfileInputOutput.h \
  resultsCPP.h mathconst.h beamCPP.h beamParametersCPP.h parametersCPP.h \
  readData.h memory.h particleBeamCPP.h particlesCPP.h \
  physicalTools.h switchesCPP.h jetParameterCPP.h grvCPP.h pairsCPP.h \
  meshCPP.h

mathematicalTools.o: mathematicalTools.cc mathematicalTools.h besselCPP.h \
  splineCPP.h typeDefs.h rndmCPP.h define.h physconst.h

lumiCPP.o: lumiCPP.cc lumiCPP.h rndmCPP.h define.h meshCPP.h \
  abstractParticle.h typeDefs.h mathematicalEntities.h \
  mathematicalTools.h besselCPP.h splineCPP.h physconst.h \
  abstractIOclass.h fileInputOutput.h fileInputOutputAscii.h \
  IfileInputOutput.h

jetParameterCPP.o: jetParameterCPP.cc jetParameterCPP.h typeDefs.h \
  define.h physconst.h

guineapigCPP.o: guineapigCPP.cc guineapigCPP.h rndmCPP.h define.h \
  beamParametersCPP.h abstractParticle.h typeDefs.h \
  mathematicalEntities.h mathematicalTools.h besselCPP.h splineCPP.h \
  physconst.h abstractIOclass.h parametersCPP.h readData.h memory.h \
  gridCPP.h meshCPP.h beamCPP.h particleBeamCPP.h \
  particlesCPP.h physicalTools.h mathconst.h fileInputOutput.h \
  fileInputOutputAscii.h IfileInputOutput.h extraPhotonCPP.h \
  switchesCPP.h resultsCPP.h pairsCPP.h lumiCPP.h jetParameterCPP.h \
  backgroundCPP.h minijetsCPP.h grvCPP.h bhabhaSamplesCPP.h fourierCPP.h \
  fieldCPP.h timerCPP.h

grvCPP.o: grvCPP.cc grvCPP.h

gridCPP.o: gridCPP.cc gridCPP.h typeDefs.h meshCPP.h rndmCPP.h \
  define.h beamCPP.h beamParametersCPP.h abstractParticle.h \
  mathematicalEntities.h mathematicalTools.h besselCPP.h splineCPP.h \
  physconst.h abstractIOclass.h parametersCPP.h readData.h memory.h \
  particleBeamCPP.h particlesCPP.h physicalTools.h mathconst.h \
  fileInputOutput.h fileInputOutputAscii.h IfileInputOutput.h \
  extraPhotonCPP.h switchesCPP.h resultsCPP.h pairsCPP.h lumiCPP.h \
  jetParameterCPP.h backgroundCPP.h minijetsCPP.h grvCPP.h \
  bhabhaSamplesCPP.h fourierCPP.h fieldCPP.h

fourierCPP.o: fourierCPP.cc fourierCPP.h

fieldCPP.o: fieldCPP.cc fieldCPP.h typeDefs.h define.h mathconst.h \
  fourierCPP.h

bhabhaSamplesCPP.o: bhabhaSamplesCPP.cc bhabhaSamplesCPP.h particlesCPP.h \
  define.h typeDefs.h abstractParticle.h mathematicalEntities.h \
  mathematicalTools.h besselCPP.h splineCPP.h rndmCPP.h physconst.h \
  abstractIOclass.h physicalTools.h pairsCPP.h resultsCPP.h mathconst.h \
  beamCPP.h beamParametersCPP.h parametersCPP.h readData.h memory.h \
  particleBeamCPP.h fileInputOutput.h fileInputOutputAscii.h \
  IfileInputOutput.h switchesCPP.h meshCPP.h

timerCPP.o: timerCPP.cc timerCPP.h define.h

switchesCPP.o: switchesCPP.cc switchesCPP.h typeDefs.h define.h \
  abstractParticle.h mathematicalEntities.h mathematicalTools.h \
  besselCPP.h splineCPP.h rndmCPP.h physconst.h abstractIOclass.h \
  parametersCPP.h readData.h memory.h

splineCPP.o: splineCPP.cc splineCPP.h fileInputOutput.h \
  abstractParticle.h typeDefs.h mathematicalEntities.h \
  mathematicalTools.h besselCPP.h rndmCPP.h define.h physconst.h \
  abstractIOclass.h fileInputOutputAscii.h IfileInputOutput.h

rndmCPP.o: rndmCPP.cc rndmCPP.h define.h mathconst.h

resultsCPP.o: resultsCPP.cc resultsCPP.h define.h physconst.h \
  mathconst.h beamCPP.h beamParametersCPP.h abstractParticle.h typeDefs.h \
  mathematicalEntities.h mathematicalTools.h besselCPP.h splineCPP.h \
  rndmCPP.h abstractIOclass.h parametersCPP.h readData.h memory.h \
  particleBeamCPP.h particlesCPP.h physicalTools.h fileInputOutput.h \
  fileInputOutputAscii.h IfileInputOutput.h switchesCPP.h

principal.o: principal.cc guineapigCPP.h rndmCPP.h define.h \
  beamParametersCPP.h abstractParticle.h typeDefs.h \
  mathematicalEntities.h mathematicalTools.h besselCPP.h splineCPP.h \
  physconst.h abstractIOclass.h parametersCPP.h readData.h memory.h \
  gridCPP.h meshCPP.h beamCPP.h particleBeamCPP.h \
  particlesCPP.h physicalTools.h mathconst.h fileInputOutput.h \
  fileInputOutputAscii.h IfileInputOutput.h extraPhotonCPP.h \
  switchesCPP.h resultsCPP.h pairsCPP.h lumiCPP.h jetParameterCPP.h \
  backgroundCPP.h minijetsCPP.h grvCPP.h bhabhaSamplesCPP.h fourierCPP.h \
  fieldCPP.h timerCPP.h

physicalTools.o: physicalTools.cc physicalTools.h besselCPP.h splineCPP.h \
  typeDefs.h rndmCPP.h define.h physconst.h mathematicalEntities.h \
  mathematicalTools.h abstractIOclass.h

particlesCPP.o: particlesCPP.cc particlesCPP.h define.h typeDefs.h \
  abstractParticle.h mathematicalEntities.h mathematicalTools.h \
  besselCPP.h splineCPP.h rndmCPP.h physconst.h abstractIOclass.h \
  physicalTools.h tridentCPP.h

particleBeamCPP.o: particleBeamCPP.cc particleBeamCPP.h particlesCPP.h \
  define.h typeDefs.h abstractParticle.h mathematicalEntities.h \
  mathematicalTools.h besselCPP.h splineCPP.h rndmCPP.h physconst.h \
  abstractIOclass.h physicalTools.h mathconst.h fileInputOutput.h \
  fileInputOutputAscii.h IfileInputOutput.h

readData.o: readData.cc readData.h

memory.o: memory.cc memory.h

tridentCPP.o: tridentCPP.cc tridentCPP.h physicalTools.h rndmCPP.h \
  gridCPP.h \

clean :

	rm -f $(MODULES) guinea guinearoot
