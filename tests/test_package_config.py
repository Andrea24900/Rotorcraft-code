"""Tests for package configuration and metadata."""

import sys
from pathlib import Path
import tomli  # For Python < 3.11, use tomli; for 3.11+ use tomllib


def test_pyproject_toml_exists():
    """Test that pyproject.toml exists in the project root."""
    project_root = Path(__file__).parent.parent
    pyproject_path = project_root / "pyproject.toml"
    assert pyproject_path.exists(), "pyproject.toml should exist"


def test_pyproject_toml_valid():
    """Test that pyproject.toml is valid TOML and has required fields."""
    project_root = Path(__file__).parent.parent
    pyproject_path = project_root / "pyproject.toml"

    # Read and parse TOML
    try:
        if sys.version_info >= (3, 11):
            import tomllib
            with open(pyproject_path, "rb") as f:
                config = tomllib.load(f)
        else:
            import tomli
            with open(pyproject_path, "rb") as f:
                config = tomli.load(f)
    except Exception as e:
        raise AssertionError(f"Failed to parse pyproject.toml: {e}")

    # Test build-system section
    assert "build-system" in config, "build-system section is required"
    assert "requires" in config["build-system"], "build-system.requires is required"
    assert "build-backend" in config["build-system"], "build-backend is required"
    assert config["build-system"]["build-backend"] == "setuptools.build_meta"

    # Test project section
    assert "project" in config, "project section is required"
    project = config["project"]

    # Required fields
    assert "name" in project, "project.name is required"
    assert project["name"] == "rotor_dynamics", f"Expected 'rotor_dynamics', got '{project['name']}'"

    assert "version" in project, "project.version is required"
    assert project["version"] == "0.1.0", f"Expected '0.1.0', got '{project['version']}'"

    assert "description" in project, "project.description is required"
    assert len(project["description"]) > 0, "description should not be empty"

    assert "requires-python" in project, "project.requires-python is required"
    assert project["requires-python"] == ">=3.8"


def test_dependencies_defined():
    """Test that dependencies are properly defined."""
    project_root = Path(__file__).parent.parent
    pyproject_path = project_root / "pyproject.toml"

    if sys.version_info >= (3, 11):
        import tomllib
        with open(pyproject_path, "rb") as f:
            config = tomllib.load(f)
    else:
        import tomli
        with open(pyproject_path, "rb") as f:
            config = tomli.load(f)

    dependencies = config["project"]["dependencies"]

    # Check required dependencies
    assert any("numpy" in dep for dep in dependencies), "numpy should be in dependencies"
    assert any("scipy" in dep for dep in dependencies), "scipy should be in dependencies"
    assert any("matplotlib" in dep for dep in dependencies), "matplotlib should be in dependencies"

    # Check version constraints
    numpy_deps = [dep for dep in dependencies if "numpy" in dep]
    assert any(">=1.24" in dep for dep in numpy_deps), "numpy version should be >=1.24.0"


def test_optional_dependencies():
    """Test that optional dependencies are defined."""
    project_root = Path(__file__).parent.parent
    pyproject_path = project_root / "pyproject.toml"

    if sys.version_info >= (3, 11):
        import tomllib
        with open(pyproject_path, "rb") as f:
            config = tomllib.load(f)
    else:
        import tomli
        with open(pyproject_path, "rb") as f:
            config = tomli.load(f)

    optional_deps = config["project"].get("optional-dependencies", {})

    # Check dev dependencies
    assert "dev" in optional_deps, "dev optional dependencies should be defined"
    dev_deps = optional_deps["dev"]
    assert any("pytest" in dep for dep in dev_deps), "pytest should be in dev dependencies"


def test_authors_defined():
    """Test that authors are properly defined."""
    project_root = Path(__file__).parent.parent
    pyproject_path = project_root / "pyproject.toml"

    if sys.version_info >= (3, 11):
        import tomllib
        with open(pyproject_path, "rb") as f:
            config = tomllib.load(f)
    else:
        import tomli
        with open(pyproject_path, "rb") as f:
            config = tomli.load(f)

    authors = config["project"]["authors"]
    assert len(authors) > 0, "At least one author should be defined"

    first_author = authors[0]
    assert "name" in first_author, "Author should have a name"
    assert "email" in first_author, "Author should have an email"


def test_classifiers_defined():
    """Test that classifiers are properly defined."""
    project_root = Path(__file__).parent.parent
    pyproject_path = project_root / "pyproject.toml"

    if sys.version_info >= (3, 11):
        import tomllib
        with open(pyproject_path, "rb") as f:
            config = tomllib.load(f)
    else:
        import tomli
        with open(pyproject_path, "rb") as f:
            config = tomli.load(f)

    classifiers = config["project"]["classifiers"]
    assert len(classifiers) > 0, "Classifiers should be defined"

    # Check for common classifier categories
    has_dev_status = any("Development Status" in c for c in classifiers)
    has_intended_audience = any("Intended Audience" in c for c in classifiers)
    has_language = any("Programming Language :: Python" in c for c in classifiers)

    assert has_dev_status, "Should have Development Status classifier"
    assert has_intended_audience, "Should have Intended Audience classifier"
    assert has_language, "Should have Programming Language classifier"


def test_package_discovery_configured():
    """Test that package discovery is properly configured."""
    project_root = Path(__file__).parent.parent
    pyproject_path = project_root / "pyproject.toml"

    if sys.version_info >= (3, 11):
        import tomllib
        with open(pyproject_path, "rb") as f:
            config = tomllib.load(f)
    else:
        import tomli
        with open(pyproject_path, "rb") as f:
            config = tomli.load(f)

    # Check setuptools configuration
    assert "tool" in config, "tool section should exist"
    assert "setuptools" in config["tool"], "tool.setuptools should exist"
    assert "packages" in config["tool"]["setuptools"] or "packages.find" in config["tool"]["setuptools"], \
        "Package discovery should be configured"


def test_license_defined():
    """Test that license is defined."""
    project_root = Path(__file__).parent.parent
    pyproject_path = project_root / "pyproject.toml"

    if sys.version_info >= (3, 11):
        import tomllib
        with open(pyproject_path, "rb") as f:
            config = tomllib.load(f)
    else:
        import tomli
        with open(pyproject_path, "rb") as f:
            config = tomli.load(f)

    assert "license" in config["project"], "License should be defined"
    license_value = config["project"]["license"]

    # Can be either a string (modern) or a dict (old style)
    if isinstance(license_value, str):
        assert len(license_value) > 0, "License string should not be empty"
    elif isinstance(license_value, dict):
        assert "text" in license_value, "License dict should have 'text' field"


def test_readme_defined():
    """Test that README is defined and exists."""
    project_root = Path(__file__).parent.parent
    pyproject_path = project_root / "pyproject.toml"

    if sys.version_info >= (3, 11):
        import tomllib
        with open(pyproject_path, "rb") as f:
            config = tomllib.load(f)
    else:
        import tomli
        with open(pyproject_path, "rb") as f:
            config = tomli.load(f)

    readme = config["project"].get("readme")
    if readme:
        readme_path = project_root / readme
        assert readme_path.exists(), f"README file {readme} should exist"


if __name__ == "__main__":
    # Run tests manually
    import traceback

    tests = [
        test_pyproject_toml_exists,
        test_pyproject_toml_valid,
        test_dependencies_defined,
        test_optional_dependencies,
        test_authors_defined,
        test_classifiers_defined,
        test_package_discovery_configured,
        test_license_defined,
        test_readme_defined,
    ]

    passed = 0
    failed = 0

    for test in tests:
        try:
            test()
            print(f"PASS {test.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"FAIL {test.__name__}: {e}")
            failed += 1
        except Exception as e:
            print(f"FAIL {test.__name__}: {type(e).__name__}: {e}")
            traceback.print_exc()
            failed += 1

    print(f"\n{passed} passed, {failed} failed")
    sys.exit(0 if failed == 0 else 1)
