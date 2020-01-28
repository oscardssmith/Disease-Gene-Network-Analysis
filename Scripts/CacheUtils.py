import os
import os.path as path
import pickle
import tempfile

def compute_if_not_cached(f, fileName = None, *args):
    """
     Tool for caching large calculation.
     When run, it checks if a temp file containing the result exists.
         This file is either filename,
         or the name of function if fileName not specified
     If so, it unpickles and returns that file.
     If not, it runs f(*args), pickles the result to the file and returns it.
    """
    if fileName is None:
        fileName = f.__name__
    path = path.join(tempfile.gettempdir(), "CompBioCompsCache", fileName + ".pickle")
    if os.path.isfile(path):
        print("pickled matrix file exists, loading matrix from file")
        try:
            with open(path, 'rb') as handle:
                return pickle.load(handle)
        except pickle.UnpicklingError: #If previous pickling was broken or corrupted.
            os.remove(path)
            return compute_if_not_cached(f, fileName, *args)
    else:
        result = f(*args)
        with open(path, 'wb') as handle:
            pickle.dump(result, handle)
        return result

