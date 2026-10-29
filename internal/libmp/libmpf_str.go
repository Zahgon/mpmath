// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package libmp

import (
	"fmt"
	"math"
	"math/big"
	"strconv"
	"strings"
	"unicode/utf8"
)

// Radix conversion, to_str, from_str and format-spec formatting from
// mpmath/libmp/libmpf.py.

func powerOfTwoBase(base int) (int, bool) {
	e2 := int(math.Log2(float64(base)))
	return e2, 1<<e2 == base
}

// ToDigitsExp returns (sign, digits, exponent) of s with dps digits in the
// given base, the digits truncated toward zero. dps may be fractional.
func ToDigitsExp(s MPF, dps float64, base int) (string, string, int) {
	sign := ""
	if s.Sign == 1 {
		sign = "-"
		s = MpfNeg(s, 0, RoundFast)
	}
	man, exp, bc := s.Man, s.Exp, s.BC
	if man.Sign() == 0 {
		return "", strings.Repeat("0", int(dps)), 0
	}
	var blog2 float64
	if base == 10 {
		blog2 = blog2_10
	} else if e2, ok := powerOfTwoBase(base); ok {
		blog2 = float64(e2)
	} else {
		Raise("NotImplementedError", "")
	}
	bitprec := int(dps*blog2) + 10
	expFrom1 := exp + bc
	exponent := 0
	if base == 10 && absInt(expFrom1) > 3500 {
		expprec := intBitLen(exp) + 5
		tmp := FromInt64(int64(exp))
		tmp = MpfMul(tmp, MpfLn2(expprec, RoundFast), 0, RoundFast)
		tmp = MpfDiv(tmp, MpfLn10(expprec, RoundFast), expprec, RoundFast)
		b := ToInt(tmp, RoundFast)
		s = MpfDiv(s, MpfPowInt(Ften, b, bitprec, RoundFast), bitprec, RoundFast)
		man, exp, bc = s.Man, s.Exp, s.BC
		exponent = int(b.Int64())
	}
	fixprec := maxInt(bitprec-exp-bc, 0)
	fixdps := int(float64(fixprec)/blog2 + 0.5)
	sf := ToFixed(s, fixprec)
	sb := BinToRadix(sf, fixprec, base, fixdps)
	digits := Numeral(sb, base)
	exponent += len(digits) - fixdps - 1
	return sign, digits, exponent
}

func digitVal(c byte) int {
	switch {
	case c >= '0' && c <= '9':
		return int(c - '0')
	case c >= 'a' && c <= 'z':
		return int(c-'a') + 10
	}
	return int(c-'A') + 10
}

// RoundDigits rounds a digit string to dps digits; it returns the digits and
// how far the point moved.
func RoundDigits(sign bool, digits string, dps, base int, rnd Rnd, fixed bool) (string, int) {
	if len(digits) <= dps {
		Raise("AssertionError", "")
	}
	if rnd == RoundCeiling {
		if sign {
			rnd = RoundDown
		} else {
			rnd = RoundUp
		}
	} else if rnd == RoundFloor {
		if sign {
			rnd = RoundUp
		} else {
			rnd = RoundDown
		}
	}
	exponent := 0
	var rndDigs string
	switch rnd {
	case RoundDown:
		return digits[:dps], 0
	case RoundNearest:
		rndDigs = stddigits[base/2+base%2 : base]
	default:
		rndDigs = stddigits[1:base]
	}
	tieDown, tieUp := false, false
	if rnd == RoundNearest {
		if digits[dps] == rndDigs[0] {
			tieDown = true
			for i := dps + 1; i < len(digits); i++ {
				if digits[i] != '0' {
					tieDown = false
					break
				}
			}
			prev := dps - 1
			if prev < 0 {
				prev += len(digits)
			}
			odd := ""
			for i := 1; i < base; i += 2 {
				odd += string(stddigits[i])
			}
			if strings.IndexByte(odd, digits[prev]) >= 0 {
				tieDown = false
			}
		}
	} else if rnd == RoundUp {
		if digits[dps] == '0' {
			for i := dps + 1; i < len(digits); i++ {
				if digits[i] != '0' {
					tieUp = true
					break
				}
			}
		}
	}
	if tieDown {
		digits = digits[:dps] + string(stddigits[digitVal(digits[dps])-1])
	} else if tieUp {
		digits = digits[:dps] + "1"
	}
	if strings.IndexByte(rndDigs, digits[dps]) >= 0 {
		digits = digits[:dps]
		i := dps - 1
		dig := stddigits[base-1]
		for i >= 0 && digits[i] == dig {
			i--
		}
		if i >= 0 {
			digits = digits[:i] + string(stddigits[digitVal(digits[i])+1]) + strings.Repeat("0", dps-i-1)
		} else {
			n := dps - 1
			if fixed {
				n = dps
			}
			digits = "1" + strings.Repeat("0", maxInt(n, 0))
			exponent++
		}
	} else {
		digits = digits[:dps]
	}
	return digits, exponent
}

// StrOpts are to_str's keyword arguments; nil bounds take the defaults.
type StrOpts struct {
	StripZeros       bool
	MinFixed         *float64
	MaxFixed         *float64
	ShowZeroExponent bool
	Base             int
	BinaryExp        bool
	Rnd              Rnd
}

// DefaultStrOpts are to_str's defaults.
func DefaultStrOpts() StrOpts {
	return StrOpts{StripZeros: true, Base: 10, Rnd: RoundNearest}
}

// ToStr is to_str(s, dps, ...).
func ToStr(s MPF, dps int, o StrOpts) string {
	base := o.Base
	sep := "e"
	if base > 10 {
		sep = "@"
	}
	if o.BinaryExp {
		sep = "p"
		if base != 2 && base != 16 {
			Raise("ValueError", "binary_exp option could be used for base 2 and 16")
		}
	}
	switch o.Rnd {
	case RoundNearest, RoundFloor, RoundCeiling, RoundUp, RoundDown:
	default:
		Raise("ValueError", "rnd should be one of n, f, c, u, d.")
	}
	prefix := ""
	switch base {
	case 2:
		prefix = "0b"
	case 8:
		prefix = "0o"
	case 16:
		prefix = "0x"
	}
	if s.Man.Sign() == 0 {
		switch {
		case s.Eq(Fzero):
			t := "0.0"
			if dps == 0 {
				t = ".0"
			}
			if o.ShowZeroExponent {
				t += sep + "+0"
			}
			return prefix + t
		case s.Eq(Finf):
			return "inf"
		case s.Eq(Fninf):
			return "-inf"
		case s.Eq(Fnan):
			return "nan"
		}
		Raise("ValueError", "")
	}
	minFixed := float64(minInt(-floorDivInt(dps, 3), -5))
	if o.MinFixed != nil {
		minFixed = *o.MinFixed
	}
	maxFixed := float64(dps)
	if o.MaxFixed != nil {
		maxFixed = *o.MaxFixed
	}
	sign, digits, exponent := ToDigitsExp(s, float64(dps+10), base)
	rndDigs := stddigits[base/2+base%2 : base]
	if dps == 0 {
		if strings.IndexByte(rndDigs, digits[0]) >= 0 {
			exponent++
		}
		digits = ".0"
	} else {
		if o.BinaryExp && base == 16 {
			exponent *= 4
			if d0 := digitVal(digits[0]); d0 > 1 {
				shift := int(math.Floor(math.Log2(float64(d0))))
				exponent += shift
				n, _ := new(big.Int).SetString(digits, 16)
				digits = rsh(n, shift).Text(16)
			}
		}
		var expAdd int
		digits, expAdd = RoundDigits(s.Sign == 1, digits, dps, base, o.Rnd, false)
		exponent += expAdd
		var split int
		if !o.BinaryExp && minFixed < float64(exponent) && float64(exponent) < maxFixed {
			if exponent < 0 {
				digits = strings.Repeat("0", -exponent) + digits
				split = 1
			} else {
				split = exponent + 1
				if split > dps {
					digits += strings.Repeat("0", split-dps)
				}
			}
			exponent = 0
		} else {
			split = 1
		}
		digits = digits[:split] + "." + digits[split:]
		if o.StripZeros {
			digits = strings.TrimRight(digits, "0")
			if digits[len(digits)-1] == '.' {
				digits += "0"
			}
		}
	}
	sign += prefix
	if exponent == 0 && dps != 0 && !o.ShowZeroExponent {
		return sign + digits
	}
	return sign + digits + sep + fmt.Sprintf("%+d", exponent)
}

// StrToManExp is str_to_man_exp.
func StrToManExp(x string, base int) (*big.Int, int) {
	x = strings.ReplaceAll(strings.TrimRight(PyLower(x), "l"), "_", "")
	sep := "e"
	if base > 10 {
		sep = "@"
	}
	e2, pow2 := powerOfTwoBase(base)
	if pow2 && (e2 == 1 || e2 == 4) && strings.Contains(x, "p") {
		sep = "p"
	}
	parts := strings.Split(x, sep)
	exp := 0
	switch len(parts) {
	case 1:
	case 2:
		x = parts[0]
		exp = int(PyInt(parts[1], 10).Int64())
	default:
		Raise("ValueError", "couldn't convert a str to mpf")
	}
	parts = strings.Split(x, ".")
	if len(parts) == 2 {
		a, b := parts[0], strings.TrimRight(parts[1], "0")
		n := utf8.RuneCountInString(b)
		if sep != "p" {
			exp -= n
		} else {
			exp -= n * e2
		}
		if a == "" {
			a = "0"
		}
		x = a + b
	}
	return PyInt(x, base), exp
}

var specialStr = map[string]MPF{
	"inf": Finf, "+inf": Finf, "-inf": Fninf, "nan": Fnan,
	"oo": Finf, "+oo": Finf, "-oo": Fninf,
}

// FromStr is from_str(x, prec, rnd, base); base 0 detects a 0b/0o/0x prefix.
func FromStr(x string, prec int, rnd Rnd, base int) MPF {
	x = PyStrip(PyLower(x))
	if v, ok := specialStr[x]; ok {
		return v
	}
	if base == 0 {
		switch {
		case strings.HasPrefix(x, "0b") || strings.HasPrefix(x, "-0b"):
			base = 2
		case strings.HasPrefix(x, "0x") || strings.HasPrefix(x, "-0x"):
			base = 16
		case strings.HasPrefix(x, "0o") || strings.HasPrefix(x, "-0o"):
			base = 8
		default:
			base = 10
		}
	}
	if strings.Contains(x, "/") {
		parts := strings.Split(x, "/")
		if len(parts) > 2 {
			Raise("ValueError", "too many values to unpack (expected 2)")
		}
		p, q := strings.TrimRight(parts[0], "l"), strings.TrimRight(parts[1], "l")
		pi := PyInt(p, base)
		qi := PyInt(q, base)
		return FromRational(pi, qi, prec, rnd)
	}
	man, exp := StrToManExp(x, base)
	if base == 10 {
		if absInt(exp) > 400 {
			s := FromInt(man, prec+10, RoundFast)
			return MpfMul(s, MpfPowInt64(Ften, int64(exp), prec+10, RoundFast), prec, rnd)
		}
		if exp >= 0 {
			return FromInt(mul(man, pow10(exp)), prec, rnd)
		}
		return FromRational(man, pow10(-exp), prec, rnd)
	}
	if e2, ok := powerOfTwoBase(base); ok {
		if !strings.Contains(x, "p") {
			return FromManExp(man, exp*e2, prec, rnd)
		}
		return FromManExp(man, exp, prec, rnd)
	}
	Raise("NotImplementedError", "")
	return Fzero
}

// MpfHash is hash() of an mpf (CPython's numeric hash; not defined for nan).
func MpfHash(s MPF) int64 {
	if s.Man.Sign() == 0 {
		if s.Eq(Finf) {
			return 314159
		}
		if s.Eq(Fninf) {
			return -314159
		}
		if s.Eq(Fzero) {
			return 0
		}
	}
	modulus := sub(lsh(bigOne, 61), bigOne)
	h := pyMod(s.Man, modulus)
	sexp := s.Exp
	if sexp >= 0 {
		sexp = sexp % 61
	} else {
		sexp = 61 - 1 - pyModInt(-1-sexp, 61)
	}
	h = pyMod(lsh(h, sexp), modulus)
	if s.Sign == 1 {
		h = neg(h)
	}
	if h.Cmp(I(-1)) == 0 {
		h = I(-2)
	}
	return h.Int64()
}

// FormatSpec is a parsed format specification.
type FormatSpec struct {
	Fill        string
	Align       byte
	Sign        byte
	NoNeg0      bool
	Alternate   bool
	ThousandSep string
	FracSep     string
	Width       int
	Precision   int
	Type        string
	Rounding    Rnd
	HasRounding bool
}

func isAlign(r rune) bool { return r == '<' || r == '>' || r == '=' || r == '^' }

// ReadFormatSpec parses a format specification exactly as the original's
// regular expression does.
func ReadFormatSpec(spec string) FormatSpec {
	runes := []rune(spec)
	type start struct {
		fill  string
		align rune
		pos   int
	}
	var starts []start
	if len(runes) >= 2 && isAlign(runes[1]) {
		starts = append(starts, start{string(runes[0]), runes[1], 2})
	}
	if len(runes) >= 1 && isAlign(runes[0]) {
		starts = append(starts, start{"", runes[0], 1})
	}
	starts = append(starts, start{"", 0, 0})
	for _, st := range starts {
		if f, ok := parseSpecRest(runes, st.pos); ok {
			f.Fill, f.Align = " ", '>'
			if st.fill != "" {
				f.Fill = st.fill
			}
			if st.align != 0 {
				f.Align = byte(st.align)
			}
			if f.zeropad {
				if st.align == 0 {
					f.Align = '='
				}
				if st.fill == "" {
					f.Fill = "0"
				}
			}
			t := strings.ToLower(f.Type)
			if f.Precision < 0 && t != "" && t != "a" && t != "b" {
				f.Precision = 6
			}
			return f.FormatSpec
		}
	}
	Raise("ValueError", "Invalid format specifier '%s'", spec)
	return FormatSpec{}
}

type parsedSpec struct {
	FormatSpec
	zeropad bool
}

func parseSpecRest(r []rune, i int) (parsedSpec, bool) {
	f := parsedSpec{FormatSpec: FormatSpec{Sign: '-', Width: -1, Precision: -1}}
	at := func(k int) rune {
		if k < len(r) {
			return r[k]
		}
		return 0
	}
	if c := at(i); c == '-' || c == '+' || c == ' ' {
		f.Sign = byte(c)
		i++
	}
	if at(i) == 'z' {
		f.NoNeg0 = true
		i++
	}
	if at(i) == '#' {
		f.Alternate = true
		i++
	}
	if at(i) == '0' {
		k := i + 1
		for at(k) == '0' {
			k++
		}
		if c := at(k); c >= '1' && c <= '9' {
			f.zeropad = true
			i++
		}
	}
	w := i
	for c := at(i); c >= '0' && c <= '9'; c = at(i) {
		i++
	}
	if i > w {
		f.Width, _ = strconv.Atoi(string(r[w:i]))
	}
	if c := at(i); c == ',' || c == '_' {
		f.ThousandSep = string(c)
		i++
	}
	if at(i) == '.' {
		if c := at(i + 1); c == ',' || c == '_' || (c >= '0' && c <= '9') {
			i++
			p := i
			for c := at(i); c >= '0' && c <= '9'; c = at(i) {
				i++
			}
			if i > p {
				f.Precision, _ = strconv.Atoi(string(r[p:i]))
			}
			if c := at(i); c == ',' || c == '_' {
				f.FracSep = string(c)
				i++
			}
		}
	}
	switch at(i) {
	case 'U':
		f.Rounding, f.HasRounding = RoundCeiling, true
		i++
	case 'D':
		f.Rounding, f.HasRounding = RoundFloor, true
		i++
	case 'Y':
		f.Rounding, f.HasRounding = RoundUp, true
		i++
	case 'Z':
		f.Rounding, f.HasRounding = RoundDown, true
		i++
	case 'N':
		f.Rounding, f.HasRounding = RoundNearest, true
		i++
	}
	if c := at(i); i < len(r) && strings.ContainsRune("aAbeEfFgG%", c) {
		f.Type = string(c)
		i++
	}
	return f, i == len(r)
}

func calcPadding(nchars, width int, align byte) (int, int) {
	ntotal := maxInt(nchars, width)
	switch align {
	case '>', '=':
		return ntotal - nchars, 0
	case '^':
		l := (ntotal - nchars) / 2
		return l, ntotal - nchars - l
	}
	return 0, ntotal - nchars
}

func formatFixed(s MPF, dps int, rnd Rnd) (string, string) {
	_, _, exponent := ToDigitsExp(s, 1, 10)
	_, digits, exponent := ToDigitsExp(s, float64(maxInt(dps+exponent+4, int(float64(s.BC)/blog2_10))), 10)
	origDps := dps
	dps += exponent + 1
	if dps < 0 {
		return "0", strings.Repeat("0", origDps)
	}
	digits, expAdd := RoundDigits(s.Sign == 1, digits, dps, 10, rnd, true)
	exponent += expAdd
	var split int
	if exponent < 0 {
		digits = strings.Repeat("0", -exponent) + digits
		split = 1
	} else {
		split = exponent + 1
	}
	if split > len(digits) {
		split = len(digits)
	}
	intPart := digits[:split]
	if origDps == 0 {
		return intPart, ""
	}
	return intPart, digits[split:]
}

func formatScientific(s MPF, dps int, rnd Rnd) (string, string, string) {
	dps++
	_, digits, exponent := ToDigitsExp(s, float64(maxInt(dps+10, int(float64(s.BC)/blog2_10)+10)), 10)
	digits, expAdd := RoundDigits(s.Sign == 1, digits, dps, 10, rnd, false)
	exponent += expAdd
	return digits[:1], digits[1:], fmt.Sprintf("e%+03d", exponent)
}

func formatHexadecimal(s MPF, dps int, rnd Rnd) (string, string, string) {
	prec := s.Man.BitLen()
	if dps >= 0 {
		prec = 4*dps + 1
	}
	var exponent int
	var fracDigits, digits string
	if s.Man.Sign() != 0 {
		s = MpfPos(s, prec, rnd)
		exponent = s.Exp + s.BC - 1
		man := new(big.Int).Or(s.Man, lsh(bigOne, s.BC+2))
		man = lsh(man, 1+4*floorDivInt(s.BC+3, 4)-s.BC)
		fracDigits = man.Text(16)[1:]
		digits = "1"
	} else {
		digits = "0"
	}
	if dps >= 0 {
		if len(fracDigits) > dps {
			fracDigits = fracDigits[:dps]
		}
		fracDigits += strings.Repeat("0", dps-len(fracDigits))
	} else {
		fracDigits = strings.TrimRight(fracDigits, "0")
	}
	return digits, fracDigits, fmt.Sprintf("p%+d", exponent)
}

func formatBinary(s MPF, dps int, rnd Rnd) (string, string, string) {
	prec := s.Man.BitLen()
	if dps >= 0 {
		prec = dps + 1
	}
	s = MpfPos(s, prec, rnd)
	digits := s.Man.Text(2)
	if n := dps + 1 - len(digits); n > 0 {
		digits += strings.Repeat("0", n)
	}
	exponent := s.Exp
	if s.Man.Sign() != 0 {
		exponent += s.Man.BitLen() - 1
	}
	return digits[:1], digits[1:], fmt.Sprintf("p%+d", exponent)
}

func fillSep(digits, sep, prev string, nmod int) string {
	var parts []string
	for pos := nmod; pos < len(digits); pos += 3 {
		end := pos + 3
		if end > len(digits) {
			end = len(digits)
		}
		parts = append(parts, digits[pos:end])
	}
	return prev + strings.Join(parts, sep)
}

func formatDigits(num MPF, f FormatSpec, prec int, rnd Rnd, prettyReprDps bool) (string, string) {
	capitalize := strings.Contains("AFGE", f.Type) && f.Type != ""
	fmtType := strings.ToLower(f.Type)
	percent := false
	if fmtType == "%" {
		percent = true
		fmtType = "f"
		num = MpfMul(num, FromInt64(100), prec, RoundNearest)
	}
	dps := f.Precision
	intPart, exponent := "", ""
	stripLastZero, stripZeros := false, false
	if f.HasRounding {
		rnd = f.Rounding
	}
	if fmtType == "" || fmtType == "g" {
		if !f.Alternate {
			stripZeros = true
			if fmtType == "g" {
				stripLastZero = true
			}
		}
		if dps < 0 {
			if prettyReprDps {
				dps = ReprDps(prec)
			} else {
				dps = PrecToDps(prec)
			}
		}
		if dps == 0 {
			dps = 1
		}
		_, tdigits, exp := ToDigitsExp(num, math.Max(53/blog2_10, float64(dps)), 10)
		if num.Man.Sign() != 0 {
			_, expAdd := RoundDigits(true, tdigits, dps, 10, rnd, false)
			exp += expAdd
		}
		fix0 := 1
		if fmtType != "" {
			fix0 = 0
		}
		if -4 <= exp && exp < dps-fix0 {
			dps = maxInt(0, dps-exp-1)
		} else {
			fmtType = "e"
			dps = maxInt(0, dps-1)
		}
	}
	var fracPart string
	switch {
	case num.Eq(Finf) || num.Eq(Fninf) || num.Eq(Fnan):
		fracPart = "inf"
		if num.Eq(Fnan) {
			fracPart = "nan"
		}
		if capitalize {
			fracPart = strings.ToUpper(fracPart)
		}
	case fmtType == "e":
		intPart, fracPart, exponent = formatScientific(num, dps, rnd)
		if stripZeros {
			fracPart = strings.TrimRight(fracPart, "0")
		}
		if fracPart != "" || f.Alternate {
			fracPart = "." + fracPart
		}
		if capitalize {
			exponent = strings.ReplaceAll(exponent, "e", "E")
		}
	case fmtType == "a":
		intPart, fracPart, exponent = formatHexadecimal(num, dps, rnd)
		if capitalize {
			intPart = "0X" + intPart
			fracPart = strings.ToUpper(fracPart)
			exponent = strings.ReplaceAll(exponent, "p", "P")
		} else {
			intPart = "0x" + intPart
		}
		if fracPart != "" || f.Alternate {
			fracPart = "." + fracPart
		}
	case fmtType == "b":
		intPart, fracPart, exponent = formatBinary(num, dps, rnd)
		if fracPart != "" || f.Alternate {
			fracPart = "." + fracPart
		}
	default:
		intPart, fracPart = formatFixed(num, dps, rnd)
		if stripZeros {
			fracPart = strings.TrimRight(fracPart, "0")
		}
		if fracPart == "" && fmtType == "" {
			fracPart = "0"
		}
		if fracPart != "" || f.Alternate || (dps != 0 && !stripLastZero) {
			fracPart = "." + fracPart
		}
	}
	if f.FracSep != "" && fracPart != "" {
		fracPart = fillSep(fracPart, f.FracSep, fracPart[:1], 1)
	}
	digits := fracPart + exponent
	sign := ""
	if num.Sign == 1 {
		sign = "-"
	}
	if sign != "-" && f.Sign != '-' {
		sign = string(f.Sign)
	}
	if fmtType == "f" && f.NoNeg0 {
		allZero := true
		for _, c := range digits {
			if !strings.ContainsRune("0._,", c) {
				allZero = false
				break
			}
		}
		if intPart == "0" && allZero {
			if f.Sign == '-' {
				sign = ""
			} else {
				sign = string(f.Sign)
			}
		}
	}
	if percent {
		digits += "%"
	}
	sep := f.ThousandSep
	width := f.Width
	minLeading := width - len(digits) - len(sign)
	if intPart != "" && fmtType != "a" && fmtType != "b" && f.Fill == "0" && f.Align == '=' &&
		minLeading > len(intPart) {
		n := minLeading
		if sep != "" {
			n = 3*minLeading/4 + 1
		}
		if len(intPart) < n {
			intPart = strings.Repeat("0", n-len(intPart)) + intPart
		}
	}
	split := len(intPart)
	if sep != "" && split > 3 {
		nmod := split % 3
		prev := ""
		if nmod != 0 {
			prev = intPart[:nmod] + sep
		}
		intPart = fillSep(intPart, sep, prev, nmod)
	}
	return sign, intPart + digits
}

// FormatMpf is format(mpf, spec).
func FormatMpf(num MPF, spec string, prec int, rnd Rnd, prettyReprDps bool) string {
	f := ReadFormatSpec(spec)
	sign, digits := formatDigits(num, f, prec, rnd, prettyReprDps)
	nchars := len(digits) + len(sign)
	lpad, rpad := calcPadding(nchars, f.Width, f.Align)
	if f.Align == '=' {
		return sign + strings.Repeat(f.Fill, lpad) + digits + strings.Repeat(f.Fill, rpad)
	}
	return strings.Repeat(f.Fill, lpad) + sign + digits + strings.Repeat(f.Fill, rpad)
}

// FormatMpc is format(mpc, spec).
func FormatMpc(re, im MPF, spec string, prec int, rnd Rnd, prettyReprDps bool) string {
	f := ReadFormatSpec(spec)
	if f.Fill == "0" {
		Raise("ValueError", "Zero padding is not allowed in complex format specifier.")
	}
	if f.Align == '=' {
		Raise("ValueError", "'=' alignment flag is not allowed in complex format specifier.")
	}
	if f.Type == "%" {
		Raise("ValueError", "%s", "'%' formatting type is not allowed in complex format specifier.")
	}
	fmtType := strings.ToLower(f.Type)
	if fmtType == "" {
		f.Type = "g"
	}
	signRe, digitsRe := formatDigits(re, f, prec, rnd, prettyReprDps)
	fmtSign := f.Sign
	f.Sign = '+'
	signIm, digitsIm := formatDigits(im, f, prec, rnd, prettyReprDps)
	digitsIm += "j"
	if fmtType == "" {
		if re.Eq(Fzero) {
			signRe, digitsRe = "", ""
			if signIm == "+" {
				if fmtSign == ' ' || fmtSign == '+' {
					signIm = string(fmtSign)
				} else {
					signIm = ""
				}
			}
		} else {
			signRe = "(" + signRe
			digitsIm += ")"
		}
	}
	nchars := len(signRe) + len(digitsRe) + len(signIm) + len(digitsIm)
	lpad, rpad := calcPadding(nchars, f.Width, f.Align)
	return strings.Repeat(f.Fill, lpad) + signRe + digitsRe + signIm + digitsIm + strings.Repeat(f.Fill, rpad)
}
