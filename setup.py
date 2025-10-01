from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
import sys
import pybind11

ext_modules = [
    Extension(
        "moldyn_core",
        sources=[
            "cpp/bindings.cpp",
            "cpp/backend.cpp",
        ],
        include_dirs=[pybind11.get_include(), "cpp"],
        language="c++",
        extra_compile_args=["-O3", "-std=c++17"]
    )
]

class BuildExt(build_ext):
    def build_extensions(self):
        ct = self.compiler.compiler_type
        if ct == "msvc":
            for e in self.extensions:
                e.extra_compile_args = ["/O2", "/std:c++17"]
        super().build_extensions()

setup(
    name="moldynnet",
    version="0.1.0",
    packages=find_packages("src"),
    package_dir={"": "src"},
    install_requires=[],
    ext_modules=ext_modules,
    cmdclass={"build_ext": BuildExt},
    entry_points={
        "console_scripts": [
            "moldynnet = moldynnet.cli:app",
        ]
    },
    include_package_data=True,
)
