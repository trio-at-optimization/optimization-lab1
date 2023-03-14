import java.util.Random;

public final class Generator {
	private final static Random RANDOMIZER = new Random((int) (Integer.MAX_VALUE * Math.random()));
	private final int n, k;

	public Generator(final int n, final int k) {
		this.n = n;
		this.k = k;
	}

	public int[] gen() {
		final int mx = RANDOMIZER.nextInt(k, Integer.MAX_VALUE);
		final int mn = (mx / k);
		final int[] matrix = new int[n];
		matrix[0] = mn;
		matrix[n - 1] = mx;
		for (int i = 1; i < n - 1; i++) {
			matrix[i] = RANDOMIZER.nextInt(mn + 1, mx - 1);
		}
		return matrix;
	}
}
