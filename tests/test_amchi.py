"""amchi tests."""

import amchi


def test_stub() -> None:
    """Stub test to ensure the test suite runs."""
    print(amchi.__version__)  # noqa: T201


def test__greet() -> None:
    """Test the greet function."""
    assert amchi.greet("World") == "Hello, World!"


def test__greet_jim() -> None:
    """Test the greet_jim function."""
    assert amchi.greet_jim() == "Hello, Jim!"
