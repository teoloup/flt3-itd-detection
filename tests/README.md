# Tests directory
This directory contains test files for the FLT3 ITD Detection Pipeline.

## Running Tests

```bash
# Run all tests
python -m pytest tests/

# Run with coverage
python -m pytest --cov=. tests/

# Run specific test file
python -m pytest tests/test_flt3_itd.py

# Run with verbose output
python -m pytest -v tests/
```

## Test Structure

- `test_flt3_itd.py` - Main functionality tests
- Add more test files as needed for specific modules

## Test Data

Test data files should be placed in a `test_data/` subdirectory (not tracked in git).
Use small, anonymized BAM files for testing.
