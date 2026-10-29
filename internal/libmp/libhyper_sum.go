// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package libmp

import "math/big"

// HypCoeff is one parameter of a generalized hypergeometric series. Kind is
// 'Z' (integer), 'Q' (rational P/Q), 'R' (real) or 'C' (complex).
type HypCoeff struct {
	Kind byte
	Int  *big.Int
	P, Q *big.Int
	R    MPF
	C    MPC
}

func fixedOf(s MPF, wp int) *big.Int {
	m := s.Man
	if s.Sign != 0 {
		m = neg(m)
	}
	off := s.Exp + wp
	if off >= 0 {
		return lsh(m, off)
	}
	return rsh(m, -off)
}

func prodI(xs []*big.Int) *big.Int {
	r := I(1)
	for _, x := range xs {
		r = mul(r, x)
	}
	return r
}

// Hypsum sums pFq(coeffs; z) in fixed point exactly as the summators built by
// the original's make_hyp_summator do. z is real when zc is nil. It returns
// the real or complex sum, whether the result is complex, and its magnitude.
// maxterms <= 0 selects the default of wp*100.
func Hypsum(p, q int, coeffs []HypCoeff, zr MPF, zc *MPC, prec, wp, epsshift int,
	magnitudeCheck map[int]int, maxterms int) (MPF, MPC, bool, int) {
	haveComplexParam := false
	for _, c := range coeffs {
		if c.Kind == 'C' {
			haveComplexParam = true
		}
	}
	haveComplexArg := zc != nil
	haveComplex := haveComplexParam || haveComplexArg

	MAX := maxterms
	if MAX <= 0 {
		MAX = wp * 100
	}
	if epsshift < 0 || wp < 0 {
		Raise("ValueError", "negative shift count")
	}
	HIGH := lsh(bigOne, epsshift)
	LOW := neg(HIGH)
	one := lsh(bigOne, wp)
	SRE, PRE := one, one
	SIM, PIM := I(0), I(0)

	var ZRE, ZIM *big.Int
	if haveComplexArg {
		ZRE = fixedOf(zc.Re, wp)
		ZIM = fixedOf(zc.Im, wp)
	} else {
		ZRE = fixedOf(zr, wp)
	}

	var aint, bint []*big.Int
	var aratP, aratQ, bratP, bratQ []*big.Int
	var areal, breal []*big.Int
	var acre, acim, bcre, bcim []*big.Int
	for i, c := range coeffs {
		isB := i >= p
		switch c.Kind {
		case 'Z':
			v := new(big.Int).Set(c.Int)
			if isB {
				bint = append(bint, v)
			} else {
				aint = append(aint, v)
			}
		case 'Q':
			if isB {
				bratP = append(bratP, new(big.Int).Set(c.P))
				bratQ = append(bratQ, c.Q)
			} else {
				aratP = append(aratP, new(big.Int).Set(c.P))
				aratQ = append(aratQ, c.Q)
			}
		case 'R':
			v := fixedOf(c.R, wp)
			if isB {
				breal = append(breal, v)
			} else {
				areal = append(areal, v)
			}
		case 'C':
			re, im := fixedOf(c.C.Re, wp), fixedOf(c.C.Im, wp)
			if isB {
				bcre, bcim = append(bcre, re), append(bcim, im)
			} else {
				acre, acim = append(acre, re), append(acim, im)
			}
		default:
			Raise("ValueError", "")
		}
	}
	nCancel := minInt(len(areal), len(breal))
	haveMultiplier := len(aint)+len(aratP)+len(bratP) > 0

	for n := 1; n < 100000000; n++ {
		if _, ok := magnitudeCheck[n]; ok {
			pMag := PRE.BitLen()
			if haveComplex {
				pMag = maxInt(pMag, PIM.BitLen())
			}
			magnitudeCheck[n] = wp - pMag
		}
		var mulv *big.Int
		if haveMultiplier {
			mulv = mul(mul(prodI(aint), prodI(aratP)), prodI(bratQ))
		}
		div := muli(mul(mul(prodI(bint), prodI(bratP)), prodI(aratQ)), int64(n))
		if div.Sign() == 0 {
			if haveMultiplier && mulv.Sign() == 0 {
				break
			}
			Raise("ZeroDivisionError", "")
		}

		applyReal := func(P *big.Int) *big.Int {
			for k := 0; k < nCancel; k++ {
				P = floorDiv(mul(P, areal[k]), breal[k])
			}
			for _, a := range areal[nCancel:] {
				P = rsh(mul(P, a), wp)
			}
			for _, b := range breal[nCancel:] {
				P = floorDiv(lsh(P, wp), b)
			}
			return P
		}

		if haveComplex {
			PRE = applyReal(PRE)
			PIM = applyReal(PIM)
			if haveComplexArg {
				re := sub(mul(PRE, ZRE), mul(PIM, ZIM))
				im := add(mul(PIM, ZRE), mul(PRE, ZIM))
				if haveMultiplier {
					re, im = mul(mulv, re), mul(mulv, im)
				}
				PRE = rsh(floorDiv(re, div), wp)
				PIM = rsh(floorDiv(im, div), wp)
			} else {
				if haveMultiplier {
					PRE = floorDiv(rsh(mul(mul(mulv, PRE), ZRE), wp), div)
					PIM = floorDiv(rsh(mul(mul(mulv, PIM), ZRE), wp), div)
				} else {
					PRE = floorDiv(rsh(mul(PRE, ZRE), wp), div)
					PIM = floorDiv(rsh(mul(PIM, ZRE), wp), div)
				}
			}
			for i := range acre {
				PRE, PIM = sub(mul(PRE, acre[i]), mul(PIM, acim[i])), add(mul(PIM, acre[i]), mul(PRE, acim[i]))
				PRE = rsh(PRE, wp)
				PIM = rsh(PIM, wp)
			}
			for i := range bcre {
				mag := add(mul(bcre[i], bcre[i]), mul(bcim[i], bcim[i]))
				re := add(mul(PRE, bcre[i]), mul(PIM, bcim[i]))
				im := sub(mul(PIM, bcre[i]), mul(PRE, bcim[i]))
				PRE = floorDiv(lsh(re, wp), mag)
				PIM = floorDiv(lsh(im, wp), mag)
			}
		} else {
			PRE = applyReal(PRE)
			if haveMultiplier {
				PRE = floorDiv(rsh(mul(mul(PRE, mulv), ZRE), wp), div)
			} else {
				PRE = floorDiv(rsh(mul(PRE, ZRE), wp), div)
			}
		}

		if haveComplex {
			SRE = add(SRE, PRE)
			SIM = add(SIM, PIM)
			if HIGH.Cmp(PRE) > 0 && PRE.Cmp(LOW) > 0 && HIGH.Cmp(PIM) > 0 && PIM.Cmp(LOW) > 0 {
				break
			}
		} else {
			SRE = add(SRE, PRE)
			if HIGH.Cmp(PRE) > 0 && PRE.Cmp(LOW) > 0 {
				break
			}
		}
		if n > MAX {
			Raise("NoConvergence", "Hypergeometric series converges too slowly. Try increasing maxterms.")
		}
		for i := range aint {
			aint[i] = addi(aint[i], 1)
		}
		for i := range bint {
			bint[i] = addi(bint[i], 1)
		}
		for i := range aratP {
			aratP[i] = add(aratP[i], aratQ[i])
		}
		for i := range bratP {
			bratP[i] = add(bratP[i], bratQ[i])
		}
		for i := range areal {
			areal[i] = add(areal[i], one)
		}
		for i := range breal {
			breal[i] = add(breal[i], one)
		}
		for i := range acre {
			acre[i] = add(acre[i], one)
		}
		for i := range bcre {
			bcre[i] = add(bcre[i], one)
		}
	}

	if haveComplex {
		a := FromManExp(SRE, -wp, prec, RoundNearest)
		b := FromManExp(SIM, -wp, prec, RoundNearest)
		var magn int
		switch {
		case SRE.Sign() != 0 && SIM.Sign() != 0:
			magn = maxInt(a.Exp+a.BC, b.Exp+b.BC)
		case SRE.Sign() != 0:
			magn = a.Exp + a.BC
		case SIM.Sign() != 0:
			magn = b.Exp + b.BC
		default:
			magn = -wp + 1
		}
		return MPF{}, MPC{a, b}, true, magn
	}
	a := FromManExp(SRE, -wp, prec, RoundNearest)
	magn := -wp + 1
	if SRE.Sign() != 0 {
		magn = a.Exp + a.BC
	}
	return a, MPC{}, false, magn
}
