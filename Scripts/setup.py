import subprocess
import sys

def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", "--user", package])


libraries = ['termcolor', 'networkx', 'numpy', 'scipy', 'matplotlib', 'requests', 'tk']


if __name__ == "__main__":
    for library in libraries:
        install(library)