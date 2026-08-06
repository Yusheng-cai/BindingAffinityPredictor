# Tests

Add focused unit tests for schemas, parsing, unit conversion, score direction, filtering, aggregation, and provenance capture. Small adapter smoke tests should use fixtures rather than downloading model weights during the test suite.

Run the dependency-free initial suite from the repository root:

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -v
```
