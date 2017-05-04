# ===============================
# Makefile for Medida
# ===============================

ROOTGLIBS     = $(shell root-config --glibs) -lcln -lginac -lMathCore -lMathMore
#ROOTGLIBS     = -L/usr/lib/root -lGui -lCore -lCint -lPostscript -lRIO -lNet -lHist -lGraf -lGpad -lTree -lRint -lMatrix -lMathMore -lMathCore -lThread -pthread -lm -ldl -rdynamic -lcln -lginac
ROOTCFLAGS    = $(shell root-config --cflags) 

CXX           = g++
CXXFLAGS      = -g -Wall -fPIC -std=c++11
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)

# ===============================
teste: main.cpp model.h BGL.h MCMC.h widthcalc.h multivector.h Formulas.h
# -------------------------------
	$(CXX) $(CXXFLAGS) -c main.cpp -o teste.o 
	$(CXX) teste.o -o teste $(ROOTGLIBS) 

# ===============================
update: update.cpp
# -------------------------------
	$(CXX) $(CXXFLAGS) -c update.cpp -o update.o 
	$(CXX) update.o -o update $(ROOTGLIBS) 

# ===============================
draw: draw.cpp
# -------------------------------
	$(CXX) $(CXXFLAGS) -c draw.cpp -o draw.o 
	$(CXX) draw.o -o draw $(ROOTGLIBS) 

# ===============================
scriptplots: scriptplots.cpp model.h BGL.h MCMC.h widthcalc.h multivector.h Formulas.h
# -------------------------------
	$(CXX) $(CXXFLAGS) -c scriptplots.cpp -o scriptplots.o 
	$(CXX) scriptplots.o -o scriptplots $(ROOTGLIBS) 

# ===============================
clean:
	rm teste.o update.o draw.o
# ===============================

