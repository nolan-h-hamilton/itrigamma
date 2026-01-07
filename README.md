# itrigamma

A self-contained implementation of the inverse-trigamma function used in variational inference, empirical bayes, etc.

* Written in C with fast bindings for Python (`python/citrigamma.pyx`) and R (`r/RCall_itrigamma.c`).

**Install**

`python -m pip install itrigamma`

**Test**

```
python -c "import time, numpy as np
from itrigamma import itrigamma
x = np.random.default_rng(0).random(10_000_000)
t = time.perf_counter()
y = itrigamma(x)
dt = time.perf_counter() - t
print(f'{(x.size/dt) / 1_000_000:.1f} million evals per sec.')
"
```

--> `8.2 million evals per sec.`
