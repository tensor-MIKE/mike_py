# MIKE: Module Isogeny Key Exchange

A SageMath implementation of the isogeny based NIKE: MIKE

## Parameter sets

This implementation currently has the three parameter sets targetting NIST levels I, III, V for the "FastMIKE" family.

## Testing

To install dependencies and run MIKE tests run the command

```
sage -python -m pip install -r requirements.txt
sage -python -m pytest
```

If `sage -python` is not recognised (e.g. Sage installed via conda or pip), then
you should be able to run the tests directly within python:

```
python -m pip install -r requirements.txt
python -m pytest
```

## Implementation Notes

The purpose of this implementation is as a clear description of the paper and in particular has not been written to squeeze out performance from slow parts of SageMath. All finite field arithmetic and (essentially all) elliptic curve arithmetic is performed using SageMath functions rather than faster, custom x-only arithmetic.
