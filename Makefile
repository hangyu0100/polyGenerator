make_dir := $(shell mkdir -p ./obj)  

CPPFLAGS = -I./include -I/usr/include/flint -D _TBB
CXXFLAGS = -Wall -Wextra -pedantic -O2 -std=c++11 -g 
LDLIBS = -lpplp -lglpk -lm -lapron -lgmpxx -lpolkaMPQ -lmpfr -lgmp -lflint -ltbb

SOURCES = $(wildcard src/*.cc)
OBJS = $(patsubst %.cc, obj/%.o, $(notdir ${SOURCES}))

all: genPoly

./obj/%.o: ./src/%.cc ./include/*.h 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

genPoly: ${OBJS}
	$(CXX) $(LDFLAGS) $+ $(LDLIBS) -o genPoly

.PHONY:clean 
clean:
	rm -rf obj/*.o genPoly

