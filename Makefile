.PHONY: build test reproduce venv clean

VENV := .venv
VENV_PYTHON := $(VENV)/bin/python3

build:
	cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j

test: build
	ctest --test-dir build --output-on-failure

# The plotting virtualenv. `make reproduce` needs it (run_all.sh renders a
# figure after every experiment) and .venv/ is gitignored, so on a fresh
# clone nothing creates it -- `make reproduce` used to fail immediately.
# It is a real prerequisite of the target, so it is built like one.
#
# uv if it is on PATH, plain venv + pip otherwise: `make reproduce` is the
# literal command in the spec's definition of done, and it should not
# require a tool the spec never asked for.
$(VENV_PYTHON):
	@if command -v uv >/dev/null 2>&1; then \
	  echo "bootstrapping $(VENV) with uv"; \
	  uv venv $(VENV) && uv pip install --python $(VENV_PYTHON) -r experiments/requirements.txt; \
	else \
	  echo "uv not found; bootstrapping $(VENV) with python3 -m venv"; \
	  python3 -m venv $(VENV) && $(VENV_PYTHON) -m pip install --quiet --upgrade pip && \
	  $(VENV_PYTHON) -m pip install --quiet -r experiments/requirements.txt; \
	fi

venv: $(VENV_PYTHON)

reproduce: build venv
	bash experiments/run_all.sh

clean:
	rm -rf build
