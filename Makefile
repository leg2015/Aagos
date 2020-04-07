# Project-specific settings
PROJECT := Aagos
PROJECT_TEST := AagosTests
EMP_DIR := ../Empirical/source

# Flags to use regardless of compiler
CFLAGS_all := -Wall -Wno-unused-function -std=c++17 -I$(EMP_DIR)/

# Native compiler information
CXX_nat := g++-9
CFLAGS_nat := -O3 -DNDEBUG $(CFLAGS_all)
CFLAGS_nat_debug := -g -pedantic -DEMP_TRACK_MEM  -Wnon-virtual-dtor -Wcast-align -Woverloaded-virtual -Wconversion -Weffc++ $(CFLAGS_all)
CFLAGS_nat_profile := -O3 -DNDEBUG $(CFLAGS_all) -pg

# Emscripten compiler information
CXX_web := emcc
OFLAGS_web_all := -s TOTAL_MEMORY=67108864 --js-library $(EMP_DIR)/web/library_emp.js -s EXPORTED_FUNCTIONS="['_main', '_empCppCallback']" -s DISABLE_EXCEPTION_CATCHING=1 -s NO_EXIT_RUNTIME=1 #--embed-file configs
OFLAGS_web := -Oz -DNDEBUG
OFLAGS_web_debug := -g4 -Oz -pedantic -Wno-dollar-in-identifier-extension

CFLAGS_web := $(CFLAGS_all) $(OFLAGS_web) $(OFLAGS_web_all)
CFLAGS_web_debug := $(CFLAGS_all) $(OFLAGS_web_debug) $(OFLAGS_web_all)


default: $(PROJECT)
native: $(PROJECT)
web: $(PROJECT).js
all: $(PROJECT) $(PROJECT).js

debug:	CFLAGS_nat := $(CFLAGS_nat_debug)
debug:	$(PROJECT)



debug-web:	CFLAGS_web := $(CFLAGS_web_debug)
debug-web:	$(PROJECT).js

web-debug:	debug-web

$(PROJECT):	source/native/$(PROJECT).cc
	$(CXX_nat) $(CFLAGS_nat) source/native/$(PROJECT).cc -o $(PROJECT)
	@echo To build the web version use: make web
	@echo To build the test version use: make $(PROJECT_TEST)
	@echo To build the profile version use: make profile

profile:	CFLAGS_nat_profile := $(CFLAGS_nat_profile)
profile:    source/native/$(PROJECT).cc
	$(CXX_nat) $(CFLAGS_nat_profile) source/native/$(PROJECT).cc -o $(PROJECT)
	@echo To build the web version use: make web
	@echo To build the test version use: make $(PROJECT_TEST)




$(PROJECT_TEST): source/native/$(PROJECT_TEST).cc
	$(CXX_nat) $(CFLAGS_nat) source/native/$(PROJECT_TEST).cc -o $(PROJECT_TEST)

debugTest:	CFLAGS_nat := $(CFLAGS_nat_debug)
debugTest: source/native/$(PROJECT_TEST).cc
	$(CXX_nat) $(CFLAGS_nat) source/native/$(PROJECT_TEST).cc -o $(PROJECT_TEST)



$(PROJECT).js: source/web/$(PROJECT)-web.cc
	$(CXX_web) $(CFLAGS_web) source/web/$(PROJECT)-web.cc -o web/$(PROJECT).js

clean:
	rm -f $(PROJECT) $(PROJECT_TEST) web/$(PROJECT).js web/*.js.map web/*.js.map *~ source/*.o

# Debugging information
print-%: ; @echo '$(subst ','\'',$*=$($*))'
