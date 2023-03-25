public class Main {
	public static void main(String[] args) {
		final int[] array = new int[]{1, 2, 3, 4, 5};
		final StringBuilder sb = new StringBuilder();
		for (int i = 0; i < array.length; i++) {
			sb.append(String.format("%d * x^{2}_{%d} + ", array[i], i));
		}
		sb.delete(sb.length() - 2, sb.length() - 1);
		System.out.println(sb);
	}
}
