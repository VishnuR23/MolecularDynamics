.PHONY: build test reproduce clean

build:
	cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j

test: build
	ctest --test-dir build --output-on-failure

reproduce: build
	bash experiments/run_all.sh

clean:
	rm -rf build
