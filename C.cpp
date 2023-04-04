#include <cstdio>
#include <vector>

#define MOD 998244353

int n, m;
std::vector< int > P;
std::vector< int > Q;

int main()
{
	scanf("%i %i", &n, &m);

	int deg = std::max(n, m);

	P = std::vector< int >(deg + 1, 0);
	Q = std::vector< int >(deg + 1, 0);

	for (int i = 0; i < n + 1; i++)
	{
		scanf("%i", &P[i]);
	}

	for (int i = 0; i < m + 1; i++)
	{
		scanf("%i", &Q[i]);
	}

	printf("%i\n", deg);
	for (int i = 0; i < deg + 1; i++)
	{
		printf("%i ", (P[i] + Q[i]) % MOD);
	}
	printf("\n");

	printf("%i\n", n + m);
	std::vector< int > mult(n + m + 1, 0);
	for (int i = 0; i < n + 1; i++)
	{
		for (int k = 0; ((k + i) < (n + m + 2)) && k < m + 1; k++)
		{
			mult[k + i] += (P[i] * Q[k]) % MOD;
			mult[k + i] %= MOD;
		}
	}

	for (int i = 0; i < n + m + 1; i++)
	{
		printf("%i ", mult[i]);
	}
	printf("\n");

	std::vector< int > div(1000, 0);
	printf("%i ", P[0] / Q[0]);
	for (int i = 1; i < 1000; i++)
	{
		int sum = 0;=
		for (int j = 0; j < i; j++)
		{
			sum += (div[j] * Q[i - j]) % MOD;
		}
		sum %= MOD;
		div[i] = (P[i] - sum + MOD) % MOD;

		printf("%i ", div[i]);
	}
}
/*

 3 2
 0 1 2 3
 1 2 3


*/
