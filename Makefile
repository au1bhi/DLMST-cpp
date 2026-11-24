CXX      := g++
CXXFLAGS := -std=c++20 -O2 -Wall -Wextra -Wpedantic -Iinclude
ASANFLAGS := -std=c++20 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -Wall -Wextra -Wpedantic -Iinclude

BUILD_DIR := build

.PHONY: all test asan bench clean

all: $(BUILD_DIR)/dlmst_bench $(BUILD_DIR)/test_dlmst

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/dlmst_bench: src/main.cpp include/dlmst/*.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) src/main.cpp -o $@

$(BUILD_DIR)/test_dlmst: tests/test_dlmst.cpp include/dlmst/*.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) tests/test_dlmst.cpp -o $@

$(BUILD_DIR)/test_dlmst_asan: tests/test_dlmst.cpp include/dlmst/*.hpp | $(BUILD_DIR)
	$(CXX) $(ASANFLAGS) tests/test_dlmst.cpp -o $@

test: $(BUILD_DIR)/test_dlmst
	./$(BUILD_DIR)/test_dlmst

asan: $(BUILD_DIR)/test_dlmst_asan
	./$(BUILD_DIR)/test_dlmst_asan

bench: $(BUILD_DIR)/dlmst_bench
	./$(BUILD_DIR)/dlmst_bench

clean:
	rm -rf $(BUILD_DIR)
