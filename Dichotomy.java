import java.util.function.Function;

public class Dichotomy {
	private final Function<Double, Double> func;
	private double a, b;
	private final double eps;

	public Dichotomy(final Function<Double, Double> func, final double a, final double b, final double eps) {
		this.func = func;
		this.a = a;
		this.b = b;
		this.eps = eps;
	}

	public double gen() {
		while (Double.compare(b - a, eps) > 0) {
			final double x = (a + b) / 2;
			final double f1 = func.apply(x - eps);
			final double f2 = func.apply(x + eps);
			if (Double.compare(f1, f2) < 0) {
				b = x;
			} else {
				a = x;
			}
		}
		return b;
	}
}
