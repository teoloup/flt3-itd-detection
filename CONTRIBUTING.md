# Contributing to FLT3 ITD Detection Pipeline

Thank you for your interest in contributing! This document provides guidelines for contributing to the project.

## Development Setup

1. **Fork and Clone**
   ```bash
   git clone https://github.com/teoloup/flt3-itd-detection.git
   cd flt3-itd-detection
   ```

2. **Create Virtual Environment**
   ```bash
   python -m venv venv
   source venv/bin/activate  # Linux/Mac
   venv\Scripts\activate     # Windows
   ```

3. **Install Development Dependencies**
   ```bash
   pip install -r requirements.txt
   pip install -e .
   ```

## Code Style

- **Formatting**: Use `black` for code formatting
- **Linting**: Use `flake8` for linting
- **Type Hints**: Use type hints where possible
- **Docstrings**: Follow Google-style docstrings

```bash
# Format code
black .

# Check linting
flake8 .

# Type checking
mypy .
```

## Testing

- Write tests for new features
- Maintain test coverage above 80%
- Run tests before submitting PRs

```bash
# Run tests
pytest

# Run with coverage
pytest --cov=.
```

## Pull Request Process

1. **Create Feature Branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make Changes**
   - Write clean, documented code
   - Add tests for new functionality
   - Update documentation if needed

3. **Test Your Changes**
   ```bash
   pytest
   black .
   flake8 .
   ```

4. **Commit with Clear Messages**
   ```bash
   git commit -m "feat: add new ITD validation method"
   ```

5. **Push and Create PR**
   ```bash
   git push origin feature/your-feature-name
   ```

## Commit Message Format

Use conventional commits:
- `feat:` - New features
- `fix:` - Bug fixes
- `docs:` - Documentation changes
- `test:` - Test additions/changes
- `refactor:` - Code refactoring
- `perf:` - Performance improvements

## Issue Reporting

When reporting issues, please include:
- Python version
- Operating system
- Input file characteristics
- Complete error messages
- Steps to reproduce

## Questions?

Feel free to open an issue for questions or reach out to the maintainers.
