import subprocess, sys, os, json, pathlib

def test_cli_help():
    out = subprocess.check_output([sys.executable, '-m', 'moldynnet.cli', '--help'])
    assert b'MolDynamicsNet CLI' in out
