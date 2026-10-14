//! Portable maths for output paths. These go through the pure-Rust `libm` crate, so results are
//! bit-identical on every platform. Use them instead of `f64::ln`, `f64::exp` and friends
//! anywhere the result can reach generated data.

/// Natural logarithm.
pub fn ln(x: f64) -> f64 {
    libm::log(x)
}

/// Exponential function.
pub fn exp(x: f64) -> f64 {
    libm::exp(x)
}

/// Square root.
pub fn sqrt(x: f64) -> f64 {
    libm::sqrt(x)
}

/// Base-10 logarithm.
pub fn log10(x: f64) -> f64 {
    libm::log10(x)
}

/// `x` raised to `y`.
pub fn pow(x: f64, y: f64) -> f64 {
    libm::pow(x, y)
}

/// Rounds half away from zero to `dp` decimal places.
pub fn round_dp(x: f64, dp: u32) -> f64 {
    let m = libm::pow(10.0, f64::from(dp));
    libm::round(x * m) / m
}

/// `x` rounded to `digits` significant digits (half away from zero), as R's `signif()` for
/// values that are not exact ties.
pub fn signif(x: f64, digits: u32) -> f64 {
    if x == 0.0 || !x.is_finite() {
        return x;
    }
    let e = libm::floor(libm::log10(libm::fabs(x))) as i32;
    let dp = digits as i32 - 1 - e;
    if dp >= 0 {
        round_dp(x, dp as u32)
    } else {
        let m = libm::pow(10.0, f64::from(-dp));
        libm::round(x / m) * m
    }
}

/// Standard normal CDF, Φ(x).
pub fn norm_cdf(x: f64) -> f64 {
    0.5 * libm::erfc(-x / core::f64::consts::SQRT_2)
}

/// Standard normal quantile, Φ⁻¹(p), for p in (0, 1) (Wichura's AS241, ~1e-16 relative error).
pub fn norm_inv(p: f64) -> f64 {
    assert!(p > 0.0 && p < 1.0, "norm_inv needs 0 < p < 1, got {p}");
    let q = p - 0.5;
    if q.abs() <= 0.425 {
        let r = 0.180625 - q * q;
        return q
            * (((((((2_509.080_928_730_122_7 * r + 33_430.575_583_588_13) * r
                + 67_265.770_927_008_7)
                * r
                + 45_921.953_931_549_87)
                * r
                + 13_731.693_765_509_46)
                * r
                + 1_971.590_950_306_551_3)
                * r
                + 133.141_667_891_784_38)
                * r
                + 3.387_132_872_796_366_5)
            / (((((((5_226.495_278_852_546 * r + 28_729.085_735_721_943) * r
                + 39_307.895_800_092_71)
                * r
                + 21_213.794_301_586_597)
                * r
                + 5_394.196_021_424_751)
                * r
                + 687.187_007_492_057_9)
                * r
                + 42.313_330_701_600_91)
                * r
                + 1.0);
    }
    let mut r = if q < 0.0 { p } else { 1.0 - p };
    r = libm::sqrt(-libm::log(r));
    let val = if r <= 5.0 {
        r -= 1.6;
        (((((((7.745_450_142_783_414e-4 * r + 0.022_723_844_989_269_184) * r
            + 0.241_780_725_177_450_6)
            * r
            + 1.270_458_252_452_368_4)
            * r
            + 3.647_848_324_763_204_5)
            * r
            + 5.769_497_221_460_691)
            * r
            + 4.630_337_846_156_546)
            * r
            + 1.423_437_110_749_683_5)
            / (((((((1.050_750_071_644_416_9e-9 * r + 5.475_938_084_995_345e-4) * r
                + 0.015_198_666_563_616_457)
                * r
                + 0.148_103_976_427_480_08)
                * r
                + 0.689_767_334_985_1)
                * r
                + 1.676_384_830_183_803_8)
                * r
                + 2.053_191_626_637_759)
                * r
                + 1.0)
    } else {
        r -= 5.0;
        (((((((2.010_334_399_292_288_1e-7 * r + 2.711_555_568_743_487_6e-5) * r
            + 0.001_242_660_947_388_078_4)
            * r
            + 0.026_532_189_526_576_124)
            * r
            + 0.296_560_571_828_504_87)
            * r
            + 1.784_826_539_917_291_3)
            * r
            + 5.463_784_911_164_114)
            * r
            + 6.657_904_643_501_103)
            / (((((((2.044_263_103_389_939_7e-15 * r + 1.421_511_758_316_446e-7) * r
                + 1.846_318_317_510_054_8e-5)
                * r
                + 7.868_691_311_456_133e-4)
                * r
                + 0.014_875_361_290_850_615)
                * r
                + 0.136_929_880_922_735_8)
                * r
                + 0.599_832_206_555_888)
                * r
                + 1.0)
    };
    if q < 0.0 { -val } else { val }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normal_cdf_and_quantile() {
        assert!((norm_cdf(0.0) - 0.5).abs() < 1e-15);
        assert!((norm_cdf(1.959963984540054) - 0.975).abs() < 1e-12);
        for p in [1e-12, 1e-6, 0.01, 0.2, 0.5, 0.8, 0.975, 0.999999] {
            assert!(
                (norm_cdf(norm_inv(p)) - p).abs() < 1e-12 * p.max(1e-3),
                "{p}"
            );
        }
        assert!((norm_inv(0.975) - 1.959963984540054).abs() < 1e-12);
    }
}
