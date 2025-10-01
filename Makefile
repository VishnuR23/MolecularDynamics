.PHONY: cpp test train format clean

cpp:
	python setup.py build_ext --inplace

test:
	pytest -q

train:
	moldynnet train --config configs/train_toy.yaml

format:
	black src tests

clean:
	rm -rf build dist *.egg-info *.so *.pyd
