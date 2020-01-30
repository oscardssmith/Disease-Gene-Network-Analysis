import os
import os.path as path
import pickle
import tempfile

def compute_if_not_cached(f, *args, fileName=None):
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
    file_path = path.join(tempfile.gettempdir(), "CompBioCompsCache", fileName + ".pickle")
    if path.isfile(file_path):
        print("pickled matrix file exists, loading matrix from file")
        try:
            with open(file_path, 'rb') as handle:
                return pickle.load(handle)
        except pickle.UnpicklingError: #If previous pickling was broken or corrupted.
            os.remove(file_path)
            return compute_if_not_cached(f, fileName, *args)
    else:
        result = f(*args)
        with open(file_path, 'wb') as handle:
            pickle.dump(result, handle)
        return result

