// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package libmp

import (
	"math"
	"math/big"
)

// Port of mpmath/libmp/libmpf.py: raw binary floating-point values and their
// correctly rounded arithmetic.

// MPF is a raw mpf value (-1)**Sign * Man * 2**Exp with BC = Man.BitLen().
// The special values carry a zero mantissa and a nonzero exponent, exactly as
// the original encodes them. Man is never modified once the value is built.
type MPF struct {
	Sign int
	Man  *big.Int
	Exp  int
	BC   int
}

// Rnd is a rounding mode: 'n' nearest, 'f' floor, 'c' ceiling, 'u' up
// (away from zero) or 'd' down (toward zero).
type Rnd byte

// The rounding modes.
const (
	RoundNearest Rnd = 'n'
	RoundFloor   Rnd = 'f'
	RoundCeiling Rnd = 'c'
	RoundUp      Rnd = 'u'
	RoundDown    Rnd = 'd'
	RoundFast        = RoundDown
)

// The constant raw values of the original.
var (
	Fzero = MPF{0, I(0), 0, 0}
	Fone  = MPF{0, I(1), 0, 1}
	Fnone = MPF{1, I(1), 0, 1}
	Ftwo  = MPF{0, I(1), 1, 1}
	Ften  = MPF{0, I(5), 1, 3}
	Fhalf = MPF{0, I(1), -1, 1}
	Fnan  = MPF{0, I(0), -123, -1}
	Finf  = MPF{0, I(0), -456, -2}
	Fninf = MPF{1, I(0), -789, -3}
)

const blog2_10 = 3.3219280948873626

// Eq is tuple equality of two raw values (nan equals nan here, as for tuples).
func (s MPF) Eq(t MPF) bool {
	return s.Sign == t.Sign && s.Exp == t.Exp && s.BC == t.BC && s.Man.Cmp(t.Man) == 0
}

// IsSpecial reports a zero mantissa with a nonzero exponent (inf or nan).
func (s MPF) IsSpecial() bool { return s.Man.Sign() == 0 && s.Exp != 0 }

// IsZero reports the value zero.
func (s MPF) IsZero() bool { return s.Man.Sign() == 0 && s.Exp == 0 }

// IsNan reports nan.
func (s MPF) IsNan() bool { return s.Eq(Fnan) }

// IsInf reports +inf or -inf.
func (s MPF) IsInf() bool { return s.Eq(Finf) || s.Eq(Fninf) }

// PrecToDps is prec_to_dps.
func PrecToDps(n int) int {
	return maxInt(1, int(math.RoundToEven(float64(n)/blog2_10)-1))
}

// DpsToPrec is dps_to_prec.
func DpsToPrec(n int) int {
	return maxInt(1, int(math.RoundToEven(float64(n+1)*blog2_10)))
}

// ReprDps is repr_dps.
func ReprDps(n int) int {
	dps := PrecToDps(n)
	if dps == 15 {
		return 17
	}
	return dps + 3
}

// hmask(n) is (1 << (n-1)) - 1.
func hmask(n int) *big.Int {
	return sub(lsh(bigOne, n-1), bigOne)
}

// RoundInt rounds the integer x / 2**n in the given direction.
func RoundInt(x *big.Int, n int, rnd Rnd) *big.Int {
	switch rnd {
	case RoundNearest:
		if x.Sign() >= 0 {
			t := rsh(x, n-1)
			if isOdd(t) && (t.Bit(1) == 1 || new(big.Int).And(x, hmask(n)).Sign() != 0) {
				return addi(rsh(t, 1), 1)
			}
			return rsh(t, 1)
		}
		return neg(RoundInt(neg(x), n, rnd))
	case RoundFloor:
		return rsh(x, n)
	case RoundCeiling:
		return neg(rsh(neg(x), n))
	case RoundDown:
		if x.Sign() >= 0 {
			return rsh(x, n)
		}
		return neg(rsh(neg(x), n))
	case RoundUp:
		if x.Sign() >= 0 {
			return neg(rsh(neg(x), n))
		}
		return rsh(x, n)
	}
	return nil
}

// shiftsDown is shifts_down[rnd][sign].
func shiftsDown(rnd Rnd, sign int) bool {
	switch rnd {
	case RoundFloor:
		return sign == 0
	case RoundCeiling:
		return sign == 1
	case RoundDown:
		return true
	case RoundUp:
		return false
	}
	return false
}

// Normalize builds a canonical raw value from a finite sign, nonnegative
// mantissa, exponent and exact bit count, rounding to prec bits.
func Normalize(sign int, man *big.Int, exp, bc, prec int, rnd Rnd) MPF {
	if bc != man.BitLen() || man.Sign() < 0 {
		Raise("AssertionError", "")
	}
	if man.Sign() == 0 {
		return Fzero
	}
	n := bc - prec
	if n > 0 {
		if rnd == RoundNearest {
			t := rsh(man, n-1)
			if isOdd(t) && (t.Bit(1) == 1 || new(big.Int).And(man, hmask(n)).Sign() != 0) {
				man = addi(rsh(t, 1), 1)
			} else {
				man = rsh(t, 1)
			}
		} else if shiftsDown(rnd, sign) {
			man = rsh(man, n)
		} else {
			man = neg(rsh(neg(man), n))
		}
		exp += n
		bc = prec
	}
	if !isOdd(man) {
		t := int(man.TrailingZeroBits())
		man = rsh(man, t)
		exp += t
		bc -= t
	}
	if man.Cmp(bigOne) == 0 {
		bc = 1
	}
	return MPF{sign, man, exp, bc}
}

// FromManExp creates a raw value from a signed mantissa and exponent,
// rounding to prec bits (prec 0 stores it exactly).
func FromManExp(man *big.Int, exp, prec int, rnd Rnd) MPF {
	sign := 0
	if man.Sign() < 0 {
		sign = 1
		man = neg(man)
	}
	bc := man.BitLen()
	if prec == 0 {
		if man.Sign() == 0 {
			return Fzero
		}
		if !isOdd(man) {
			t := int(man.TrailingZeroBits())
			man = rsh(man, t)
			exp += t
			bc -= t
		}
		return MPF{sign, man, exp, bc}
	}
	return Normalize(sign, man, exp, bc, prec, rnd)
}

// FromInt creates a raw value from an integer.
func FromInt(n *big.Int, prec int, rnd Rnd) MPF {
	return FromManExp(n, 0, prec, rnd)
}

// FromInt64 creates an exact raw value from a machine integer.
func FromInt64(n int64) MPF { return FromManExp(big.NewInt(n), 0, 0, RoundFast) }

// ToManExp returns the mantissa (signed if asked) and exponent.
func ToManExp(s MPF, signed bool) (*big.Int, int) {
	if s.IsSpecial() {
		Raise("ValueError", "mantissa and exponent are defined for finite numbers only")
	}
	man := s.Man
	if signed && s.Sign == 1 {
		man = neg(man)
	}
	return man, s.Exp
}

// ToInt converts to an integer, rounding in the given direction.
func ToInt(s MPF, rnd Rnd) *big.Int {
	if s.IsSpecial() {
		if s.Eq(Fnan) {
			Raise("ValueError", "cannot convert nan to int")
		}
		Raise("OverflowError", "cannot convert infinity to int")
	}
	if s.Exp >= 0 {
		if s.Sign == 1 {
			return lsh(neg(s.Man), s.Exp)
		}
		return lsh(s.Man, s.Exp)
	}
	if rnd == RoundFast {
		if s.Sign == 1 {
			return neg(rsh(s.Man, -s.Exp))
		}
		return rsh(s.Man, -s.Exp)
	}
	if s.Sign == 1 {
		return RoundInt(neg(s.Man), -s.Exp, rnd)
	}
	return RoundInt(s.Man, -s.Exp, rnd)
}

func mpfRoundInt(s MPF, rnd Rnd) MPF {
	if s.IsSpecial() {
		return s
	}
	if s.Exp >= 0 {
		return s
	}
	mag := s.Exp + s.BC
	if mag < 1 {
		switch rnd {
		case RoundCeiling:
			if s.Sign == 1 {
				return Fzero
			}
			return Fone
		case RoundFloor:
			if s.Sign == 1 {
				return Fnone
			}
			return Fzero
		case RoundNearest:
			if mag < 0 || s.Man.Cmp(bigOne) == 0 {
				return Fzero
			} else if s.Sign == 1 {
				return Fnone
			}
			return Fone
		default:
			Raise("NotImplementedError", "")
		}
	}
	return MpfPos(s, minInt(s.BC, mag), rnd)
}

// MpfFloor is floor(s), rounded to prec when prec is nonzero.
func MpfFloor(s MPF, prec int, rnd Rnd) MPF {
	v := mpfRoundInt(s, RoundFloor)
	if prec != 0 {
		v = MpfPos(v, prec, rnd)
	}
	return v
}

// MpfCeil is ceil(s).
func MpfCeil(s MPF, prec int, rnd Rnd) MPF {
	v := mpfRoundInt(s, RoundCeiling)
	if prec != 0 {
		v = MpfPos(v, prec, rnd)
	}
	return v
}

// MpfNint is the nearest integer to s (ties to even).
func MpfNint(s MPF, prec int, rnd Rnd) MPF {
	v := mpfRoundInt(s, RoundNearest)
	if prec != 0 {
		v = MpfPos(v, prec, rnd)
	}
	return v
}

// MpfFrac is s - floor(s).
func MpfFrac(s MPF, prec int, rnd Rnd) MPF {
	return MpfSub(s, MpfFloor(s, 0, RoundFast), prec, rnd)
}

// FromFloat creates a raw value from a float64.
func FromFloat(x float64, prec int, rnd Rnd) MPF {
	if x != x {
		return Fnan
	}
	if math.IsInf(x, 1) {
		return Finf
	}
	if math.IsInf(x, -1) {
		return Fninf
	}
	m, e := math.Frexp(x)
	return FromManExp(big.NewInt(int64(m*(1<<53))), e-53, prec, rnd)
}

// ToFloat converts to a float64 (overflow gives inf, underflow 0).
func ToFloat(s MPF, strict bool, rnd Rnd) float64 {
	if s.Man.Sign() == 0 {
		switch {
		case s.Eq(Fzero):
			return 0
		case s.Eq(Finf):
			return math.Inf(1)
		case s.Eq(Fninf):
			return math.Inf(-1)
		}
		return math.NaN()
	}
	sign, man, exp, bc := s.Sign, s.Man, s.Exp, s.BC
	if bc > 53 {
		r := Normalize(sign, man, exp, bc, 53, rnd)
		sign, man, exp, bc = r.Sign, r.Man, r.Exp, r.BC
	}
	m := float64(man.Int64())
	if sign == 1 {
		m = -m
	}
	var v float64
	if exp > 2000 {
		v = math.Inf(1)
	} else if exp < -2200 {
		v = 0
	} else {
		v = math.Ldexp(m, exp)
	}
	if math.IsInf(v, 0) {
		if strict {
			Raise("OverflowError", "math range error")
		}
		if exp+bc > 0 {
			if sign == 1 {
				return math.Inf(-1)
			}
			return math.Inf(1)
		}
		return 0
	}
	if v == 0 && sign == 1 {
		return math.Copysign(0, -1)
	}
	return v
}

// FromRational is p/q rounded to prec bits.
func FromRational(p, q *big.Int, prec int, rnd Rnd) MPF {
	return MpfDiv(FromInt(p, 0, RoundFast), FromInt(q, 0, RoundFast), prec, rnd)
}

// ToRational returns integers p, q with s = p/q exactly.
func ToRational(s MPF) (*big.Int, *big.Int) {
	if s.Eq(Fnan) {
		Raise("ValueError", "cannot convert nan to a rational number")
	}
	if s.Eq(Finf) || s.Eq(Fninf) {
		Raise("OverflowError", "cannot convert infinity to a rational number")
	}
	man := s.Man
	if s.Sign == 1 {
		man = neg(man)
	}
	if s.Exp >= 0 {
		return lsh(man, s.Exp), I(1)
	}
	return man, lsh(bigOne, -s.Exp)
}

// ToFixed converts to a fixed-point integer with prec fractional bits.
func ToFixed(s MPF, prec int) *big.Int {
	offset := s.Exp + prec
	man := s.Man
	if s.Sign == 1 {
		man = neg(man)
	}
	return lsh(man, offset)
}

// MpfEq tests equality (false when either side is nan).
func MpfEq(s, t MPF) bool {
	if s.Man.Sign() == 0 || t.Man.Sign() == 0 {
		if s.Eq(Fnan) || t.Eq(Fnan) {
			return false
		}
	}
	return s.Eq(t)
}

// MpfSign is -1, 0 or 1 (nan gives 0).
func MpfSign(s MPF) int {
	if s.Man.Sign() == 0 {
		if s.Eq(Finf) {
			return 1
		}
		if s.Eq(Fninf) {
			return -1
		}
		return 0
	}
	if s.Sign == 1 {
		return -1
	}
	return 1
}

// MpfCmp compares s and t: -1, 0 or 1 (with the original's nan convention).
func MpfCmp(s, t MPF) int {
	if s.Man.Sign() == 0 || t.Man.Sign() == 0 {
		if s.Eq(Fzero) {
			return -MpfSign(t)
		}
		if t.Eq(Fzero) {
			return MpfSign(s)
		}
		if s.Eq(t) {
			return 0
		}
		if t.Eq(Fnan) {
			return 1
		}
		if s.Eq(Finf) {
			return 1
		}
		if t.Eq(Fninf) {
			return 1
		}
		return -1
	}
	if s.Sign != t.Sign {
		if s.Sign == 0 {
			return 1
		}
		return -1
	}
	if s.Exp == t.Exp {
		c := s.Man.Cmp(t.Man)
		if c == 0 {
			return 0
		}
		if c > 0 {
			if s.Sign == 1 {
				return -1
			}
			return 1
		}
		if s.Sign == 1 {
			return 1
		}
		return -1
	}
	a := s.BC + s.Exp
	b := t.BC + t.Exp
	if s.Sign == 1 {
		if a < b {
			return 1
		}
		if a > b {
			return -1
		}
	} else {
		if a < b {
			return -1
		}
		if a > b {
			return 1
		}
	}
	delta := MpfSub(s, t, 5, RoundFloor)
	if delta.Sign == 1 {
		return -1
	}
	return 1
}

// MpfLt is s < t.
func MpfLt(s, t MPF) bool {
	if s.Eq(Fnan) || t.Eq(Fnan) {
		return false
	}
	return MpfCmp(s, t) < 0
}

// MpfLe is s <= t.
func MpfLe(s, t MPF) bool {
	if s.Eq(Fnan) || t.Eq(Fnan) {
		return false
	}
	return MpfCmp(s, t) <= 0
}

// MpfGt is s > t.
func MpfGt(s, t MPF) bool {
	if s.Eq(Fnan) || t.Eq(Fnan) {
		return false
	}
	return MpfCmp(s, t) > 0
}

// MpfGe is s >= t.
func MpfGe(s, t MPF) bool {
	if s.Eq(Fnan) || t.Eq(Fnan) {
		return false
	}
	return MpfCmp(s, t) >= 0
}

// MpfPos rounds s to prec bits (prec 0 leaves it unchanged).
func MpfPos(s MPF, prec int, rnd Rnd) MPF {
	if prec != 0 {
		if s.IsSpecial() {
			return s
		}
		return Normalize(s.Sign, s.Man, s.Exp, s.BC, prec, rnd)
	}
	return s
}

// MpfNeg is -s.
func MpfNeg(s MPF, prec int, rnd Rnd) MPF {
	if s.Man.Sign() == 0 {
		if s.Exp != 0 {
			if s.Eq(Finf) {
				return Fninf
			}
			if s.Eq(Fninf) {
				return Finf
			}
		}
		return s
	}
	if prec == 0 {
		return MPF{1 - s.Sign, s.Man, s.Exp, s.BC}
	}
	return Normalize(1-s.Sign, s.Man, s.Exp, s.BC, prec, rnd)
}

// MpfAbs is abs(s).
func MpfAbs(s MPF, prec int, rnd Rnd) MPF {
	if s.IsSpecial() {
		if s.Eq(Fninf) {
			return Finf
		}
		return s
	}
	if prec == 0 {
		if s.Sign == 1 {
			return MPF{0, s.Man, s.Exp, s.BC}
		}
		return s
	}
	return Normalize(0, s.Man, s.Exp, s.BC, prec, rnd)
}

func orPrec(prec, bc int) int {
	if prec != 0 {
		return prec
	}
	return bc
}

func mpfAdd(s, t MPF, prec int, rnd Rnd, subtract int) MPF {
	ssign, sman, sexp, sbc := s.Sign, s.Man, s.Exp, s.BC
	tsign, tman, texp, tbc := t.Sign, t.Man, t.Exp, t.BC
	tsign ^= subtract
	if sman.Sign() != 0 && tman.Sign() != 0 {
		offset := sexp - texp
		if offset != 0 {
			if offset > 0 {
				if offset > 100 && prec != 0 {
					delta := sbc + sexp - tbc - texp
					if delta > prec+4 {
						offset = prec + 4
						sman = lsh(sman, offset)
						if tsign == ssign {
							sman = addi(sman, 1)
						} else {
							sman = addi(sman, -1)
						}
						return Normalize(ssign, sman, sexp-offset, sman.BitLen(), prec, rnd)
					}
				}
				var man *big.Int
				if ssign == tsign {
					man = add(tman, lsh(sman, offset))
				} else {
					if ssign == 1 {
						man = sub(tman, lsh(sman, offset))
					} else {
						man = sub(lsh(sman, offset), tman)
					}
					if man.Sign() >= 0 {
						ssign = 0
					} else {
						man = neg(man)
						ssign = 1
					}
				}
				bc := man.BitLen()
				return Normalize(ssign, man, texp, bc, orPrec(prec, bc), rnd)
			}
			if offset < -100 && prec != 0 {
				delta := tbc + texp - sbc - sexp
				if delta > prec+4 {
					offset = prec + 4
					tman = lsh(tman, offset)
					if ssign == tsign {
						tman = addi(tman, 1)
					} else {
						tman = addi(tman, -1)
					}
					return Normalize(tsign, tman, texp-offset, tman.BitLen(), prec, rnd)
				}
			}
			var man *big.Int
			if ssign == tsign {
				man = add(sman, lsh(tman, -offset))
			} else {
				if tsign == 1 {
					man = sub(sman, lsh(tman, -offset))
				} else {
					man = sub(lsh(tman, -offset), sman)
				}
				if man.Sign() >= 0 {
					ssign = 0
				} else {
					man = neg(man)
					ssign = 1
				}
			}
			bc := man.BitLen()
			return Normalize(ssign, man, sexp, bc, orPrec(prec, bc), rnd)
		}
		var man *big.Int
		if ssign == tsign {
			man = add(tman, sman)
		} else {
			if ssign == 1 {
				man = sub(tman, sman)
			} else {
				man = sub(sman, tman)
			}
			if man.Sign() >= 0 {
				ssign = 0
			} else {
				man = neg(man)
				ssign = 1
			}
		}
		bc := man.BitLen()
		return Normalize(ssign, man, texp, bc, orPrec(prec, bc), rnd)
	}
	if subtract != 0 {
		t = MpfNeg(t, 0, RoundFast)
	}
	if sman.Sign() == 0 {
		if sexp != 0 {
			if s.Eq(t) || t.Man.Sign() != 0 || t.Exp == 0 {
				return s
			}
			return Fnan
		}
		if t.Man.Sign() != 0 {
			return Normalize(t.Sign, t.Man, t.Exp, t.BC, orPrec(prec, t.BC), rnd)
		}
		return t
	}
	if t.Exp != 0 {
		return t
	}
	if sman.Sign() != 0 {
		return Normalize(ssign, sman, sexp, sbc, orPrec(prec, sbc), rnd)
	}
	return s
}

// MpfAdd is s + t rounded to prec bits (prec 0 is exact).
func MpfAdd(s, t MPF, prec int, rnd Rnd) MPF { return mpfAdd(s, t, prec, rnd, 0) }

// MpfSub is s - t.
func MpfSub(s, t MPF, prec int, rnd Rnd) MPF { return mpfAdd(s, t, prec, rnd, 1) }

// MpfSum sums xs with a single final rounding.
func MpfSum(xs []MPF, prec int, rnd Rnd, absolute bool) MPF {
	man := I(0)
	exp := 0
	maxExtraPrec := prec * 2
	if maxExtraPrec == 0 {
		maxExtraPrec = 1000000
	}
	var special *MPF
	for _, x := range xs {
		xman, xexp, xbc := x.Man, x.Exp, x.BC
		if xman.Sign() != 0 {
			if x.Sign == 1 && !absolute {
				xman = neg(xman)
			}
			delta := xexp - exp
			if xexp >= exp {
				if delta > maxExtraPrec && (man.Sign() == 0 || delta-man.BitLen() > maxExtraPrec) {
					man = xman
					exp = xexp
				} else {
					man = add(man, lsh(xman, delta))
				}
			} else {
				delta = -delta
				if delta-xbc > maxExtraPrec {
					if man.Sign() == 0 {
						man, exp = xman, xexp
					}
				} else {
					man = add(lsh(man, delta), xman)
					exp = xexp
				}
			}
		} else if xexp != 0 {
			if absolute {
				x = MpfAbs(x, 0, RoundFast)
			}
			base := Fzero
			if special != nil {
				base = *special
			}
			v := MpfAdd(base, x, 1, RoundFast)
			special = &v
		}
	}
	if special != nil {
		return *special
	}
	return FromManExp(man, exp, prec, rnd)
}

// MpfMul is s * t.
func MpfMul(s, t MPF, prec int, rnd Rnd) MPF {
	sign := s.Sign ^ t.Sign
	man := mul(s.Man, t.Man)
	if man.Sign() != 0 {
		bc := man.BitLen()
		if prec != 0 {
			return Normalize(sign, man, s.Exp+t.Exp, bc, prec, rnd)
		}
		return MPF{sign, man, s.Exp + t.Exp, bc}
	}
	sSpecial := s.IsSpecial()
	tSpecial := t.IsSpecial()
	if !sSpecial && !tSpecial {
		return Fzero
	}
	if s.Eq(Fnan) || t.Eq(Fnan) {
		return Fnan
	}
	if tSpecial {
		s, t = t, s
	}
	if t.Eq(Fzero) {
		return Fnan
	}
	if MpfSign(s)*MpfSign(t) == 1 {
		return Finf
	}
	return Fninf
}

// MpfMulInt multiplies by an integer.
func MpfMulInt(s MPF, n *big.Int, prec int, rnd Rnd) MPF {
	if s.Man.Sign() == 0 {
		return MpfMul(s, FromInt(n, 0, RoundFast), prec, rnd)
	}
	if n.Sign() == 0 {
		return Fzero
	}
	sign := s.Sign
	if n.Sign() < 0 {
		sign ^= 1
		n = neg(n)
	}
	man := mul(s.Man, n)
	bc := s.BC + n.BitLen() - 1
	bc += int(rsh(man, bc).Int64())
	return Normalize(sign, man, s.Exp, bc, prec, rnd)
}

// MpfShift multiplies by 2**n exactly.
func MpfShift(s MPF, n int) MPF {
	if s.Man.Sign() == 0 {
		return s
	}
	return MPF{s.Sign, s.Man, s.Exp + n, s.BC}
}

// MpfFrexp returns (y, n) with x = y*2**n and abs(y) in [0.5, 1).
func MpfFrexp(x MPF) (MPF, int) {
	if x.Man.Sign() == 0 {
		if x.Eq(Fzero) {
			return Fzero, 0
		}
		Raise("ValueError", "")
	}
	return MpfShift(x, -x.BC-x.Exp), x.BC + x.Exp
}

// MpfDiv is s / t.
func MpfDiv(s, t MPF, prec int, rnd Rnd) MPF {
	sman, sexp, sbc := s.Man, s.Exp, s.BC
	tman, texp, tbc := t.Man, t.Exp, t.BC
	if sman.Sign() == 0 || tman.Sign() == 0 {
		if s.Eq(Fzero) {
			if t.Eq(Fzero) {
				Raise("ZeroDivisionError", "")
			}
			if t.Eq(Fnan) {
				return Fnan
			}
			return Fzero
		}
		if t.Eq(Fzero) {
			Raise("ZeroDivisionError", "")
		}
		sSpecial := s.IsSpecial()
		tSpecial := t.IsSpecial()
		if sSpecial && tSpecial {
			return Fnan
		}
		if s.Eq(Fnan) || t.Eq(Fnan) {
			return Fnan
		}
		if !tSpecial {
			if t.Eq(Fzero) {
				return Fnan
			}
			if MpfSign(s)*MpfSign(t) == 1 {
				return Finf
			}
			return Fninf
		}
		return Fzero
	}
	sign := s.Sign ^ t.Sign
	if tman.Cmp(bigOne) == 0 {
		return Normalize(sign, sman, sexp-texp, sbc, prec, rnd)
	}
	var extra int
	if prec == 0 {
		extra = maxInt(sbc, tbc) - sbc + tbc + 5
	} else {
		extra = prec - sbc + tbc + 5
	}
	if extra < 5 {
		extra = 5
	}
	quot, rem := divmod(lsh(sman, extra), tman)
	if rem.Sign() != 0 {
		quot = addi(lsh(quot, 1), 1)
		extra++
	}
	bc := quot.BitLen()
	return Normalize(sign, quot, sexp-texp-extra, bc, orPrec(prec, bc), rnd)
}

// MpfRdivInt is n / t for an integer n.
func MpfRdivInt(n *big.Int, t MPF, prec int, rnd Rnd) MPF {
	if n.Sign() == 0 || t.Man.Sign() == 0 {
		return MpfDiv(FromInt(n, 0, RoundFast), t, prec, rnd)
	}
	sign := t.Sign
	if n.Sign() < 0 {
		sign ^= 1
		n = neg(n)
	}
	extra := prec + t.BC + 5
	quot, rem := divmod(lsh(n, extra), t.Man)
	if rem.Sign() != 0 {
		quot = addi(lsh(quot, 1), 1)
		extra++
	}
	return Normalize(sign, quot, -t.Exp-extra, quot.BitLen(), prec, rnd)
}

// MpfMod is s % t with Python's sign convention.
func MpfMod(s, t MPF, prec int, rnd Rnd) MPF {
	if s.IsSpecial() || t.IsSpecial() {
		if t.Eq(Finf) || t.Eq(Fninf) {
			return s
		}
		return Fnan
	}
	if s.Sign == t.Sign && t.Exp > s.Exp+s.BC {
		return s
	}
	if t.Man.Cmp(bigOne) == 0 && s.Exp > t.Exp+t.BC {
		return Fzero
	}
	base := minInt(s.Exp, t.Exp)
	sman := s.Man
	if s.Sign == 1 {
		sman = neg(sman)
	}
	tman := t.Man
	if t.Sign == 1 {
		tman = neg(tman)
	}
	man := pyMod(lsh(sman, s.Exp-base), lsh(tman, t.Exp-base))
	sign := 0
	if man.Sign() < 0 {
		man = neg(man)
		sign = 1
	}
	return Normalize(sign, man, base, man.BitLen(), prec, rnd)
}

// ReciprocalRnd maps a rounding mode to the one for a reciprocal.
func ReciprocalRnd(rnd Rnd) Rnd {
	switch rnd {
	case RoundDown:
		return RoundUp
	case RoundUp:
		return RoundDown
	case RoundFloor:
		return RoundCeiling
	case RoundCeiling:
		return RoundFloor
	}
	return RoundNearest
}

// NegativeRnd maps a rounding mode to the one for a negated value.
func NegativeRnd(rnd Rnd) Rnd {
	switch rnd {
	case RoundDown:
		return RoundDown
	case RoundUp:
		return RoundUp
	case RoundFloor:
		return RoundCeiling
	case RoundCeiling:
		return RoundFloor
	}
	return RoundNearest
}

// MpfPowInt is s**n for an integer n.
func MpfPowInt(s MPF, n *big.Int, prec int, rnd Rnd) MPF {
	if s.IsSpecial() {
		if s.Eq(Finf) {
			if n.Sign() > 0 {
				return s
			}
			if n.Sign() == 0 {
				return Fone
			}
			return Fzero
		}
		if s.Eq(Fninf) {
			if n.Sign() > 0 {
				if isOdd(n) {
					return Fninf
				}
				return Finf
			}
			if n.Sign() == 0 {
				return Fone
			}
			return Fzero
		}
		if n.Sign() == 0 {
			return Fone
		}
		return Fnan
	}
	if n.Sign() == 0 {
		return Fone
	}
	if n.Cmp(bigOne) == 0 {
		return MpfPos(s, prec, rnd)
	}
	if n.Cmp(bigTwo) == 0 {
		man, exp, bc := s.Man, s.Exp, s.BC
		if man.Sign() == 0 {
			return Fzero
		}
		man = mul(man, man)
		if man.Cmp(bigOne) == 0 {
			return MPF{0, I(1), exp + exp, 1}
		}
		bc = bc + bc - 2
		bc += rsh(man, bc).BitLen()
		return Normalize(0, man, exp+exp, bc, prec, rnd)
	}
	if n.Cmp(I(-1)) == 0 {
		return MpfDiv(Fone, s, prec, rnd)
	}
	if n.Sign() < 0 {
		inverse := MpfPowInt(s, neg(n), prec+5, ReciprocalRnd(rnd))
		return MpfDiv(Fone, inverse, prec, rnd)
	}
	resultSign := s.Sign & int(n.Bit(0))
	man, exp, bc := s.Man, s.Exp, s.BC
	if man.Cmp(bigOne) == 0 {
		return MPF{resultSign, I(1), exp * int(n.Int64()), 1}
	}
	if n.IsInt64() && n.Int64() < 1000 && int64(bc)*n.Int64() < 1000 {
		ni := int(n.Int64())
		man = powi(man, ni)
		return Normalize(resultSign, man, exp*ni, man.BitLen(), prec, rnd)
	}
	roundsDown := rnd == RoundNearest || shiftsDown(rnd, resultSign)
	workprec := prec + 4*n.BitLen() + 4
	pm, pe, pbc := I(1), 0, 1
	nn := new(big.Int).Set(n)
	for {
		if nn.Bit(0) == 1 {
			pm = mul(pm, man)
			pe = pe + exp
			pbc += bc - 2
			pbc = pbc + rsh(pm, pbc).BitLen()
			if pbc > workprec {
				if roundsDown {
					pm = rsh(pm, pbc-workprec)
				} else {
					pm = neg(rsh(neg(pm), pbc-workprec))
				}
				pe += pbc - workprec
				pbc = workprec
			}
			nn.Sub(nn, bigOne)
			if nn.Sign() == 0 {
				break
			}
		}
		man = mul(man, man)
		exp = exp + exp
		bc = bc + bc - 2
		bc = bc + rsh(man, bc).BitLen()
		if bc > workprec {
			if roundsDown {
				man = rsh(man, bc-workprec)
			} else {
				man = neg(rsh(neg(man), bc-workprec))
			}
			exp += bc - workprec
			bc = workprec
		}
		nn.Rsh(nn, 1)
	}
	return Normalize(resultSign, pm, pe, pbc, prec, rnd)
}

// MpfPowInt64 is MpfPowInt with a machine integer exponent.
func MpfPowInt64(s MPF, n int64, prec int, rnd Rnd) MPF {
	return MpfPowInt(s, big.NewInt(n), prec, rnd)
}

// MpfPerturb is x + eps with directed rounding (eps of the given sign).
func MpfPerturb(x MPF, epsSign int, prec int, rnd Rnd) MPF {
	if rnd == RoundNearest {
		return MpfPos(x, prec, rnd)
	}
	eps := MPF{epsSign, I(1), x.Exp + x.BC - prec - 1, 1}
	var away bool
	if x.Sign == 1 {
		away = (rnd == RoundDown || rnd == RoundCeiling) != (epsSign == 1)
	} else {
		away = (rnd == RoundUp || rnd == RoundCeiling) != (epsSign == 1)
	}
	if away {
		return MpfAdd(x, eps, prec, rnd)
	}
	return MpfPos(x, prec, rnd)
}

// MpfSqrt is the correctly rounded square root of a nonnegative value.
func MpfSqrt(s MPF, prec int, rnd Rnd) MPF {
	sign, man, exp, bc := s.Sign, s.Man, s.Exp, s.BC
	if sign == 1 {
		Raise("ComplexResult", "square root of a negative number")
	}
	if man.Sign() == 0 {
		return s
	}
	if exp&1 != 0 {
		exp--
		man = lsh(man, 1)
		bc++
	} else if man.Cmp(bigOne) == 0 {
		return Normalize(sign, man, floorDivInt(exp, 2), bc, prec, rnd)
	}
	shift := maxInt(4, 2*prec-bc+4)
	shift += shift & 1
	if rnd == RoundFloor || rnd == RoundDown {
		man = Isqrt(lsh(man, shift))
	} else {
		var rem *big.Int
		man, rem = Sqrtrem(lsh(man, shift))
		if rem.Sign() != 0 {
			man = addi(lsh(man, 1), 1)
			shift += 2
		}
	}
	return FromManExp(man, floorDivInt(exp-shift, 2), prec, rnd)
}

// MpfHypot is sqrt(x**2 + y**2).
func MpfHypot(x, y MPF, prec int, rnd Rnd) MPF {
	if y.Eq(Fzero) {
		return MpfAbs(x, prec, rnd)
	}
	if x.Eq(Fzero) {
		return MpfAbs(y, prec, rnd)
	}
	hypot2 := MpfAdd(MpfMul(x, x, 0, RoundFast), MpfMul(y, y, 0, RoundFast), prec+4, RoundFast)
	return MpfSqrt(hypot2, prec, rnd)
}
