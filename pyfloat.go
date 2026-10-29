// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package mpmath

import (
	"fmt"
	"math"
	"strconv"
	"strings"
)

// pyFloatFormat is CPython's float 'r' formatting; dot0 adds ".0" to
// integral fixed-point output, as repr(float) does.
func pyFloatFormat(f float64, dot0 bool) string {
	switch {
	case math.IsNaN(f):
		return "nan"
	case math.IsInf(f, 1):
		return "inf"
	case math.IsInf(f, -1):
		return "-inf"
	}
	e := strconv.FormatFloat(f, 'e', -1, 64)
	sign := ""
	if e[0] == '-' {
		sign, e = "-", e[1:]
	}
	mant, expStr, _ := strings.Cut(e, "e")
	digits := strings.Replace(mant, ".", "", 1)
	x, _ := strconv.Atoi(expStr)
	if f == 0 {
		x = 0
	}
	if x < -4 || x >= 16 {
		s := digits[:1]
		if len(digits) > 1 {
			s += "." + digits[1:]
		}
		return sign + s + fmt.Sprintf("e%+03d", x)
	}
	var s string
	switch {
	case x < 0:
		s = "0." + strings.Repeat("0", -x-1) + digits
	case len(digits) > x+1:
		s = digits[:x+1] + "." + digits[x+1:]
	default:
		s = digits + strings.Repeat("0", x+1-len(digits))
		if dot0 {
			s += ".0"
		}
	}
	return sign + s
}

// pyFloatRepr is repr(float).
func pyFloatRepr(f float64) string { return pyFloatFormat(f, true) }

// pyComplexRepr is repr(complex).
func pyComplexRepr(z complex128) string {
	re, im := real(z), imag(z)
	is := pyFloatFormat(im, false)
	if re == 0 && !math.Signbit(re) {
		return is + "j"
	}
	if !strings.HasPrefix(is, "-") {
		is = "+" + is
	}
	return "(" + pyFloatFormat(re, false) + is + "j)"
}
