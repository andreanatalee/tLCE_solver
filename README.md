# Please README


## Requirements

1. **SageMath** installed

## List of input arguments required

```bash
usage: tLCE_solver.py [-h] -n CODE_LENGTH -k CODE_DIMENSION -q PRIME -t SAMPLES

Parses command.

options:
  -h, --help            show this help message and exit
  -n CODE_LENGTH, --code_length CODE_LENGTH
                        length of code
  -k CODE_DIMENSION, --code_dimension CODE_DIMENSION
                        code dimension
  -t SAMPLES, --samples SAMPLES
                        number of samples
  -q PRIME, --prime PRIME
                        Field characteristic
```

## Script for the solver (Algorithm 3)

```bash
sage -python tLCE_solver.py -n CODE_LENGTH -k CODE_DIMENSION -q PRIME -t SAMPLES
```

### Concerning the experimental validation of the Equation 9

```bash
% sage -python test_probability.py -h
usage: test_probability.py [-h] -n CODE_LENGTH -k CODE_DIMENSION -q PRIME -t SAMPLES -m ITERATIONS

options:
  -h, --help  show this help message and exit
  -n CODE_LENGTH        Code length n
  -k CODE_DIMENSIONS    Code dimension k
  -q PRIME              Modulo q
  -t SAMPLES            Number of samples
  -m ITERATIONS         Number of experiment
```


## License

Apache License Version 2.0, January 2004
