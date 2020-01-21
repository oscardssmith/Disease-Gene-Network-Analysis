import subprocess
import sys

def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])


libraries = ['networkx', 'numpy', 'scipy', 'matplotlib', 'requests']


if __name__ == "__main__":
    for library in libraries:
        install(library)