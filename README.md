# itrigamma

A self-contained implementation of the inverse-trigamma function, $\psi'^{-1}$, for use in variational inference, empirical bayes, and other settings involving variance estimation wrt log-gamma random variables.

Written in C with fast bindings for

* Python (`python/citrigamma.pyx`)
* R (`r/RCall_itrigamma.c`)


**Install (Python)**

`pip install itrigamma`

**Test (Python)**

Copy+paste into shell:

```bash
python <<'TEST_PY'
import time, numpy as np, math
from itrigamma import itrigamma

x = np.random.default_rng(0).random(10_000_000)

t = time.perf_counter()
y = itrigamma(x)
dt = time.perf_counter() - t

print(f'Approx. {(x.size/dt) / 1_000_000:.1f} million evals per sec.')

KNOWN = [
    (1.0e-8, 100000000.499999999167),
    (0.05, 20.49583523915460591),
    (0.10, 10.491681821078422372),
    (0.20, 5.4834517798530476662),
    (0.50, 2.4599529483523074137),
    (1.00, 1.4262551202150789904),
    (2.00, 0.87666407746426017692),
    (5.00, 0.49616873470410694489),
    (10.0, 0.33508104437803564513),
    (20.0, 0.23077632915242632912),
    (50.0, 0.14337981235286147770),
]

for i, (x0, y0) in enumerate(KNOWN):
    assert math.isclose(itrigamma(x0), y0, abs_tol=1e-12)
    print(f'  --> check {i}:\t{(x0,y0)}\tGOOD')

TEST_PY
```

*Output* (Exact results vary wrt machine precision.)

```bash
TEST_PY
Approx. 8.2 million evals per sec.
  --> check 0:	(1e-08, 100000000.5)	GOOD
  --> check 1:	(0.05, 20.495835239154605)	GOOD
  --> check 2:	(0.1, 10.491681821078423)	GOOD
  --> check 3:	(0.2, 5.4834517798530475)	GOOD
  --> check 4:	(0.5, 2.4599529483523073)	GOOD
  --> check 5:	(1.0, 1.426255120215079)	GOOD
  --> check 6:	(2.0, 0.8766640774642602)	GOOD
  --> check 7:	(5.0, 0.49616873470410694)	GOOD
  --> check 8:	(10.0, 0.33508104437803565)	GOOD
  --> check 9:	(20.0, 0.23077632915242632)	GOOD
  --> check 10:	(50.0, 0.14337981235286149)	GOOD
```


**Install (R)**

First clone this repository

```bash
git clone https://github.com/nolan-h-hamilton/itrigamma.git
cd itrigamma
```

Then build with base R's [r/utils/SHLIB](https://rdrr.io/r/utils/SHLIB.html)
```bash
PKG_CPPFLAGS="-I$(pwd)/src" R CMD SHLIB r/RCall_itrigamma.c src/itrigamma.c -o r_itrigamma.so
```

Within R, you can then load the shared library and use base R's [.Call](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Interface-functions-_002eCall-and-_002eExternal) for a simple wrapper, e.g.,

```R

# load
dyn.load("r_itrigamma.so")
# easiest: include the following as a wrapper around .Call
itrigamma = function(x) .Call("r_itrigamma", as.double(x))

# e.g.,
x_ = seq(1.0e-8, 100, length.out = 50)
y = itrigamma(x_)
plot(x_,y, type='l')
```
