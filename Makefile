.PHONY: build test verify reproduce venv clean

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

# The fast path. Everything here runs in well under a minute and needs no
# Python: the unit suite, then the two experiments that gate the engine's
# correctness -- the NIST reference table and the energy-conservation
# order. This is what CI runs and what a reviewer should run first.
#
# `reproduce` below is the slow path: it re-runs every physics experiment
# from scratch, which takes minutes, and is only needed to regenerate the
# committed results rather than to check that they hold.
verify: test
	bash experiments/exp01_nist_table.sh
	bash experiments/exp02_energy_conservation.sh
	@echo
	@echo "verify: unit suite + NIST reference table + conservation order all pass."
	@echo "        run 'make reproduce' to regenerate every committed result (minutes)."

reproduce: build venv
	bash experiments/run_all.sh

clean:
	rm -rf build
