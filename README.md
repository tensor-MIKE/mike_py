# MIKE: Module Isogeny Key Exchange

A SageMath implementation of the isogeny based NIKE: MIKE

## Parameter sets

This implementation currently has the three parameter sets targetting NIST levels I, III, V for the "FastMIKE" family.

## Testing

Run tests with the command

```
sage -python -m pytest
```

## Implementation Notes

The purpose of this implementation is as a clear description of the paper and in particular has not been written to squeeze out performance from slow parts of SageMath. All finite field arithmetic and (essentially all) elliptic curve arithmetic is performed using SageMath functions rather than faster, custom x-only arithmetic.
