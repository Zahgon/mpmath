// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package libmp

import (
	"math"
	"math/big"
)

// Float helpers standing in for the C library calls the original makes
// through Python's math module and float operators. glibc returns the
// correctly rounded double for these; Go's math package can be several ulps
// away, and the original feeds such estimates into integer code. The values
// are computed at extra precision with the history-dependent caches set
// aside, so calling them leaves no trace on later results.

const pymathPrec = 160

var allMemos = []*memo{ln2Memo, ln10Memo, piMemo, eMemo, phiMemo, lnSqrt2piMemo, sqrtpiMemo}

func isolated(f func() MPF) float64 {
	type saved struct {
		prec int
		val  *big.Int
	}
	keep := make([]saved, len(allMemos))
	for i, m := range allMemos {
		keep[i] = saved{m.prec, m.val}
		m.prec, m.val = -1, nil
	}
	cacheMu.Lock()
	oldLogInt, oldLogT, oldAtanT, oldCosSin := logIntCache, logTaylorCache, atanTaylorCache, cosSinCache
	logIntCache = map[int]fixedPair{}
	logTaylorCache = map[pairKey]fixedPair{}
	atanTaylorCache = map[pairKey]fixedPair{}
	cosSinCache = map[string]fixedPair{}
	cacheMu.Unlock()
	defer func() {
		for i, m := range allMemos {
			m.prec, m.val = keep[i].prec, keep[i].val
		}
		cacheMu.Lock()
		logIntCache, logTaylorCache, atanTaylorCache, cosSinCache = oldLogInt, oldLogT, oldAtanT, oldCosSin
		cacheMu.Unlock()
	}()
	return ToFloat(f(), false, RoundNearest)
}

func finite(xs ...float64) bool {
	for _, x := range xs {
		if math.IsInf(x, 0) || math.IsNaN(x) {
			return false
		}
	}
	return true
}

func ff(x float64) MPF { return FromFloat(x, 0, RoundFast) }

// PyPow is C pow(x, y) for the finite, positive x the original passes.
func PyPow(x, y float64) float64 {
	if !finite(x, y) || x <= 0 || y == 0 || x == 1 {
		return math.Pow(x, y)
	}
	return isolated(func() MPF { return MpfPow(ff(x), ff(y), pymathPrec, RoundNearest) })
}

// PyLog is C log(x).
func PyLog(x float64) float64 {
	if !finite(x) || x <= 0 || x == 1 {
		return math.Log(x)
	}
	return isolated(func() MPF { return MpfLn(ff(x), pymathPrec, RoundNearest) })
}

// PyAtan is C atan(x). It cannot go through MpfAtan, whose cache is filled
// by atanNewton, which calls this.
func PyAtan(x float64) float64 {
	if !finite(x) || x == 0 {
		return math.Atan(x)
	}
	return isolated(func() MPF {
		wp := pymathPrec + 40
		one := lsh(bigOne, wp)
		t := ToFixed(ff(math.Abs(x)), wp)
		recip := t.Cmp(one) > 0
		if recip {
			t = floorDiv(lsh(one, wp), t)
		}
		const halvings = 8
		for i := 0; i < halvings; i++ {
			t = floorDiv(lsh(t, wp), add(one, Isqrt(add(lsh(one, wp), mul(t, t)))))
		}
		t2 := rsh(mul(t, t), wp)
		s, term := I(0), t
		for k := int64(0); term.Sign() != 0; k++ {
			q := floorDiv(term, I(2*k+1))
			if k&1 == 1 {
				s = sub(s, q)
			} else {
				s = add(s, q)
			}
			term = rsh(mul(term, t2), wp)
		}
		s = lsh(s, halvings)
		if recip {
			s = sub(rsh(PiFixed(wp), 1), s)
		}
		if x < 0 {
			s = neg(s)
		}
		return FromManExp(s, -wp, pymathPrec, RoundNearest)
	})
}

// PyAtan2 is C atan2(y, x); signed zeros and infinities go to math.Atan2,
// whose results there are exact.
func PyAtan2(y, x float64) float64 {
	if !finite(x, y) || x == 0 || y == 0 {
		return math.Atan2(y, x)
	}
	return isolated(func() MPF { return MpfAtan2(ff(y), ff(x), pymathPrec, RoundNearest) })
}

// PyCos is C cos(x).
func PyCos(x float64) float64 {
	if !finite(x) || x == 0 {
		return math.Cos(x)
	}
	return isolated(func() MPF { return MpfCos(ff(x), pymathPrec, RoundNearest) })
}

// PySin is C sin(x).
func PySin(x float64) float64 {
	if !finite(x) || x == 0 {
		return math.Sin(x)
	}
	return isolated(func() MPF { return MpfSin(ff(x), pymathPrec, RoundNearest) })
}

// PyHypot is C hypot(x, y).
func PyHypot(x, y float64) float64 {
	if !finite(x, y) || x == 0 || y == 0 {
		return math.Hypot(x, y)
	}
	return isolated(func() MPF { return MpfHypot(ff(x), ff(y), pymathPrec, RoundNearest) })
}
