import pip

def install(package):
    if hasattr(pip, 'main'):
        pip.main(['install', package])
    else:
        pip._internal.main(['install', package])


libraries = ['networkx', 'numpy', 'scipy', 'matplotlib', 'requests']


if __name__ == "__main__":
    for library in libraries:
        install(library)