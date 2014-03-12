module Language.Passage.Measure where

import Language.Passage.AST
import Language.Passage.Term(logGamma, tcase)

logit :: Floating a => a -> a
logit p = log(p/(1-p))

logBeta :: Expr -> Expr -> Expr
logBeta x y = logGamma x + logGamma y - logGamma (x + y)

logFact :: Expr -> Expr
logFact n = logGamma (n + 1)

logComb :: Expr -> Expr -> Expr
logComb n k = logFact n - logFact k - logFact (n - k)

-- | A normal Measure with mean 0 and precision 1
stdNormal :: Measure
stdNormal = Measure
  { measureName = "N(0,1)"
  , measureParams = []
  , measureSupport = Real
  , measureLogDensity = \x -> -0.5 * x**2
  }

-- | A normal Measure, with a mean and precision
normal :: Expr -> Expr -> Measure
normal m t = Measure
  { measureName    = "N"
  , measureParams  = [m, t]
  , measureSupport = Real
  , measureLogDensity      =  \x ->   log t        / 2
                         - t * (x ** 2) / 2
                         + t * x * m
                         - t * (m ** 2) / 2
  }

--- | A standard uniform Measure with parameters 0 and 1
standardUniform :: Measure
standardUniform = Measure
  { measureName    = "SU"
  , measureParams  = [0, 1]
  , measureSupport = Interval 0 1
  , measureLogDensity      = \_ -> 0
  }

--- | A uniform Measure with lower and upper bounds
uniform :: Expr -> Expr -> Measure
uniform lo hi = Measure
  { measureName    = "U"
  , measureParams  = [lo, hi]
  , measureSupport = Interval lo hi
    -- NB: Uniform Measure, independent of the variable (hence constant function)
  , measureLogDensity      = \_ -> - (log (hi - lo))
  }

discreteUniform :: Expr -> Measure
discreteUniform n = Measure
  { measureName    = "DisreteUniform"
  , measureParams  = [0, n]
  , measureSupport = Discrete (Just n)
    -- NB: Uniform Measure, independent of the variable (hence constant function)
  , measureLogDensity      = \_ -> - (log (n + 1))
  }

geometric :: Expr -> Measure
geometric p = Measure
  { measureName    = "Geometric"
  , measureParams  = [p]
  , measureSupport = Discrete Nothing
  , measureLogDensity      = \x -> x * log (1 - p) + log p
  }

-- | A categorical Measure with given support size and probabilities
-- | Probabilities are assumed to add to one (not checked here)
categorical :: Expr -> [Expr] -> Measure
categorical n ps = Measure
  { measureName    = "Categorical"
  , measureParams  = n:ps
  , measureSupport = Discrete (Just (n - 1))
  , measureLogDensity      = \x -> log (tcase x ps) 
  }

-- | A Bernoulli Measure with a mean
bernoulli :: Expr -> Measure
bernoulli p = Measure
  { measureName    = "B"
  , measureParams  = [p]
  , measureSupport = Discrete (Just 1)
  , measureLogDensity      = \x -> log (1 - p) + logit p * x
  }

-- | A binomial Measure with given number of samples and probability of success
-- | Number of samples is assumed to be fixed
binomial :: Expr -> Expr -> Measure
binomial n p = Measure
  { measureName    = "Binomial"
  , measureParams  = [n, p]
  , measureSupport = Discrete (Just n)
  , measureLogDensity      = \x -> logComb n x + x * logit p + n * log (1 - p)
  }

negBinomial :: Expr -> Expr -> Measure
negBinomial r p = Measure
  { measureName    = "NegativeBinomial"
  , measureParams  = [r, p]
  , measureSupport = PosReal
  , measureLogDensity      = \x -> logComb (x+r-1) x + r * log (1 - p) + x * log p
  }

poisson :: Expr -> Measure
poisson lambda = Measure
  { measureName    = "Poisson"
  , measureParams  = [lambda]
  , measureSupport = Discrete Nothing
  , measureLogDensity      = \x -> x * log lambda - logFact x - lambda
  }

-- | A beta Measure with the given prior sample sizes.
beta :: Expr -> Expr -> Measure
beta a b =
  Measure
    { measureName    = "Beta"
    , measureParams  = [a, b]
    , measureSupport = Interval 0 1
    , measureLogDensity      = \x -> (a - 1) * log x + (b - 1) * log (1 - x) - logBeta a b
    }

-- | A gamma Measure with the given prior sample sizes.
dgamma :: Expr -> Expr -> Measure
dgamma a b =
  Measure
    { measureName    = "Gamma"
    , measureParams  = [a, b]
    , measureSupport = PosReal
    , measureLogDensity      = \x -> a * log b - logGamma a + (a - 1) * log x - b * x
    }

-- | An improper uniform Measure; has no impact on likelihood
improperUniform :: Measure
improperUniform =
  Measure
    { measureName    = "ImproperUniform"
    , measureParams  = []
    , measureSupport = Real
    , measureLogDensity      = const 0
    }
    
-- | An improper scale
improperScale :: Measure
improperScale =
  Measure
    { measureName    = "ImproperScale"
    , measureParams  = []
    , measureSupport = PosReal
    , measureLogDensity      = \x -> -log x
    }
