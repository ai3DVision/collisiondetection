include paths.mk
CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))
AR_FLAGS := rcs
CC_FLAGS := -g -Wall -MMD -O3 -fno-omit-frame-pointer -Iinclude -msse4.2 -DNDEBUG -std=c++11

bin/libcollisions.a: $(OBJ_FILES)
	ar $(AR_FLAGS) $@ $^

obj/%.o: src/%.cpp
	g++ $(CC_FLAGS) -I$(EIGEN_DIR) -c -o $@ $<

clean:
	rm obj/*
	rm bin/*
