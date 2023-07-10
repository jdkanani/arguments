### Arguments

Simple knowledge arguments

### Installation

To install all dependencies (this project uses poetry for package management)

```bash
$ make install
```

### Usage

To run any module

```bash
$ poetry run python -m utils.utils
```

To test kate (KZG) commitment

```bash
$ poetry run python -m arguments.pcs.kzg
```

To create simple multiset argument proof and verify it

```bash
$ poetry run python -m arguments.multiset
```

To create simple permutation argument proof and verify it

```bash
$ poetry run python -m arguments.permutation
```

### Credits and references

This repository heavily uses library and methods from https://github.com/ETHorHIL/Plonk_Py
