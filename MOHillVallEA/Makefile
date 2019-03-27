CC = g++ -O3 
CFLAGS = -std=c++11 -MMD

HVEA_DIR := ./core
HVEA_SRC_FILES := $(wildcard $(HVEA_DIR)/*.cpp)
HVEA_OBJ_FILES := $(patsubst $(HVEA_DIR)/%.cpp,$(HVEA_DIR)/%.o,$(HVEA_SRC_FILES))
HVEA_DEP_FILES := $(patsubst $(HVEA_DIR)/%.cpp,$(HVEA_DIR)/%.d,$(HVEA_SRC_FILES))

all: mo-hillvallea.app

mo-hillvallea.app: mo-hillvallea.o $(HVEA_OBJ_FILES)
	$(CC) $(CFLAGS) -o $@ mo-hillvallea.o $(HVEA_OBJ_FILES)
	
%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(CEC_OBJ_FILES) $(CEC_DEP_FILES) $(HVEA_OBJ_FILES) $(HVEA_DEP_FILES) *.d *.o *.app
