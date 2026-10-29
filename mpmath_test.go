// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package mpmath

import "testing"

func text(t *testing.T, n Number, err error) string {
	t.Helper()
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	return n.String()
}

func TestValuesMatchTheOriginal(t *testing.T) {
	c := NewContext()
	c.SetDps(50)
	pi, err := c.Constant("pi", PrecOptions{})
	if got, want := text(t, pi, err), "3.1415926535897932384626433832795028841971693993751"; got != want {
		t.Errorf("str(pi) at dps 50 = %q, want %q", got, want)
	}
	c.SetPrec(53)
	cases := []struct {
		name string
		f    func() (Number, error)
		want string
	}{
		{"gamma(0.5)", func() (Number, error) { return c.Gamma("0.5") }, "1.77245385090552"},
		{"zeta(2)", func() (Number, error) { return c.Zeta(2) }, "1.64493406684823"},
		{"sqrt(-4)", func() (Number, error) { return c.Sqrt(-4) }, "(0.0 + 2.0j)"},
		{"erf(1)", func() (Number, error) { return c.Erf(1) }, "0.842700792949715"},
		{"besselj(2, 3)", func() (Number, error) { return c.Besselj(2, 3) }, "0.486091260585891"},
	}
	for _, tc := range cases {
		n, err := tc.f()
		if got := text(t, n, err); got != tc.want {
			t.Errorf("%s = %q, want %q", tc.name, got, tc.want)
		}
	}
	two, err := c.Sqrt(2)
	if err != nil {
		t.Fatal(err)
	}
	if got, want := two.Repr(), "mpf('1.4142135623730951')"; got != want {
		t.Errorf("repr(sqrt(2)) = %q, want %q", got, want)
	}
	e, _ := c.Constant("e", PrecOptions{})
	if got, _ := c.Nstr(e, 30, NstrOptions{}); got != "2.71828182845904509079559829843" {
		t.Errorf("nstr(e, 30) = %q, want %q", got, "2.71828182845904509079559829843")
	}
	x, _ := c.Mpf("1234.5678")
	if got, _ := c.Format(x, ",.2f"); got != "1,234.57" {
		t.Errorf("format(x, ',.2f') = %q, want %q", got, "1,234.57")
	}
}

func TestRoundingModes(t *testing.T) {
	c := NewContext()
	for _, tc := range []struct {
		mode string
		man  string
		exp  int
	}{{"f", "7205759403792793", -56}, {"c", "3602879701896397", -55}} {
		if err := c.SetRounding(tc.mode); err != nil {
			t.Fatal(err)
		}
		v, err := c.Mpf("0.1")
		if err != nil {
			t.Fatal(err)
		}
		_, man, exp, _ := v.Raw()
		if man.String() != tc.man || exp != tc.exp {
			t.Errorf("mpf('0.1') rounding %s = (%s, %d), want (%s, %d)", tc.mode, man, exp, tc.man, tc.exp)
		}
	}
	if err := c.SetRounding("x"); err == nil || err.Error() != "ValueError: invalid rounding mode" {
		t.Errorf("SetRounding(\"x\") error = %v, want ValueError: invalid rounding mode", err)
	}
}

func TestErrorsMatchTheOriginal(t *testing.T) {
	c := NewContext()
	for _, tc := range []struct {
		name string
		f    func() (Number, error)
		want string
	}{
		{"gamma(-2)", func() (Number, error) { return c.Gamma(-2) }, "ValueError: gamma function pole"},
		{"mpf('abc')", func() (Number, error) { return c.Mpf("abc") }, "ValueError: invalid literal for int() with base 10: 'abc'"},
		{"mpf(1)/0", func() (Number, error) { return c.Div(c.one(), 0) }, "ZeroDivisionError: "},
		{"zeta(1)", func() (Number, error) { return c.Zeta(1) }, "ValueError: zeta(1) pole"},
	} {
		_, err := tc.f()
		if err == nil || err.Error() != tc.want {
			t.Errorf("%s error = %v, want %q", tc.name, err, tc.want)
		}
	}
}
