import pytest

@pytest.fixture(scope="function")
def tmp_folder():
    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as tmpdirname:
        yield tmpdirname


