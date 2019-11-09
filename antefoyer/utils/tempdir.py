
import contextlib
import tempfile
import shutil

@contextlib.contextmanager
def temporary_directory():
    import shutil
    import tempfile
    tmp_dir = tempfile.mkdtemp()
    try:
        yield tmp_dir
    finally:
        shutil.rmtree(tmp_dir)

@contextlib.contextmanager
def temporary_cd(dir_path):
    import os
    prev_dir = os.getcwd()
    os.chdir(os.path.abspath(dir_path))
    try:
        yield
    finally:
        os.chdir(prev_dir)


