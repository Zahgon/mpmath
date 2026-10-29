// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package libmp

import (
	"fmt"
	"math/big"
	"strings"
	"unicode"
	"unicode/utf8"
)

// CPython 3.12 string behaviour that the original inherits: str.isspace,
// str.strip, repr(str) and int(text, base).

// PyIsSpace is str.isspace for one character.
func PyIsSpace(r rune) bool {
	switch {
	case r >= 0x09 && r <= 0x0d, r >= 0x1c && r <= 0x20, r == 0x85, r == 0xa0, r == 0x1680,
		r >= 0x2000 && r <= 0x200a, r == 0x2028, r == 0x2029, r == 0x202f, r == 0x205f, r == 0x3000:
		return true
	}
	return false
}

// PyStrip is str.strip() with no argument.
func PyStrip(s string) string {
	return strings.TrimFunc(s, PyIsSpace)
}

// PyLower is str.lower().
func PyLower(s string) string {
	var b strings.Builder
	for _, r := range s {
		if r == 0x130 {
			b.WriteString("i\u0307")
			continue
		}
		b.WriteRune(unicode.ToLower(r))
	}
	return b.String()
}

// PyRepr is repr() of a str.
func PyRepr(s string) string {
	quote := byte('\'')
	if strings.ContainsRune(s, '\'') && !strings.ContainsRune(s, '"') {
		quote = '"'
	}
	var b strings.Builder
	b.WriteByte(quote)
	for _, r := range s {
		switch {
		case r == rune(quote) || r == '\\':
			b.WriteByte('\\')
			b.WriteRune(r)
		case r == '\t':
			b.WriteString(`\t`)
		case r == '\n':
			b.WriteString(`\n`)
		case r == '\r':
			b.WriteString(`\r`)
		case r < 0x20 || r == 0x7f:
			fmt.Fprintf(&b, `\x%02x`, r)
		case r < 0x7f:
			b.WriteRune(r)
		case pyPrintable(r):
			b.WriteRune(r)
		case r <= 0xff:
			fmt.Fprintf(&b, `\x%02x`, r)
		case r <= 0xffff:
			fmt.Fprintf(&b, `\u%04x`, r)
		default:
			fmt.Fprintf(&b, `\U%08x`, r)
		}
	}
	b.WriteByte(quote)
	return b.String()
}

func pyPrintable(r rune) bool {
	if r == ' ' {
		return true
	}
	if unicode.In(r, unicode.Cc, unicode.Cf, unicode.Cs, unicode.Co, unicode.Zl, unicode.Zp, unicode.Zs) {
		return false
	}
	return unicode.In(r, unicode.L, unicode.M, unicode.N, unicode.P, unicode.S)
}

// decimalValue returns the digit value of a Unicode decimal digit.
func decimalValue(r rune) (int, bool) {
	if r >= '0' && r <= '9' {
		return int(r - '0'), true
	}
	if r < 0x80 || !unicode.Is(unicode.Nd, r) {
		return 0, false
	}
	start := r
	for start > 0 && unicode.Is(unicode.Nd, start-1) && r-start < 9 {
		start--
	}
	return int(r-start) % 10, true
}

// truncRepr is the %.200R of CPython's int() error message.
func truncRepr(s string) string {
	r := PyRepr(s)
	if utf8.RuneCountInString(r) > 200 {
		runes := []rune(r)
		r = string(runes[:200])
	}
	return r
}

// PyInt is int(text, base) for 2 <= base <= 36, raising ValueError as CPython does.
func PyInt(text string, base int) *big.Int {
	v, ok := pyIntParse(text, base)
	if !ok {
		Raise("ValueError", "invalid literal for int() with base %d: %s", base, truncRepr(text))
	}
	return v
}

func pyIntParse(text string, base int) (*big.Int, bool) {
	// CPython first maps Unicode decimal digits to ASCII and whitespace to spaces.
	var mapped []rune
	for _, r := range text {
		if PyIsSpace(r) {
			mapped = append(mapped, ' ')
		} else if d, ok := decimalValue(r); ok {
			mapped = append(mapped, rune('0'+d))
		} else if r >= 0x80 {
			return nil, false
		} else {
			mapped = append(mapped, r)
		}
	}
	s := strings.TrimSpace(string(mapped))
	negative := false
	if len(s) > 0 && (s[0] == '+' || s[0] == '-') {
		negative = s[0] == '-'
		s = s[1:]
	}
	if len(s) >= 2 && s[0] == '0' {
		p := s[1] | 0x20
		if (base == 16 && p == 'x') || (base == 8 && p == 'o') || (base == 2 && p == 'b') {
			s = s[2:]
			if strings.HasPrefix(s, "_") {
				s = s[1:]
			}
		}
	}
	if s == "" {
		return nil, false
	}
	var digits strings.Builder
	prevUnderscore := false
	for i := 0; i < len(s); i++ {
		c := s[i]
		if c == '_' {
			if i == 0 || prevUnderscore || i == len(s)-1 {
				return nil, false
			}
			prevUnderscore = true
			continue
		}
		prevUnderscore = false
		var d int
		switch {
		case c >= '0' && c <= '9':
			d = int(c - '0')
		case c >= 'a' && c <= 'z':
			d = int(c-'a') + 10
		case c >= 'A' && c <= 'Z':
			d = int(c-'A') + 10
		default:
			return nil, false
		}
		if d >= base {
			return nil, false
		}
		digits.WriteByte(c)
	}
	v, ok := new(big.Int).SetString(digits.String(), base)
	if !ok {
		return nil, false
	}
	if negative {
		v.Neg(v)
	}
	return v, true
}
