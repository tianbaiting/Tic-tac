// ===============================================================
// 抽取自仓库 [current]: src/utils/auxiliary.cpp
// 行号区段：453..473
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
double pi1_tilde(double p, double q, double x)
{
	return sqrt(0.25 * p * p + 9.0 / 16.0 * q * q + 0.75 * p * q * x);
}


double pi2_tilde(double p, double q, double x)
{
	return sqrt(p * p + 0.25 * q * q - p * q * x);
}

double pi1_prime_tilde(double p, double q, double x)
{
	return sqrt(0.25 * p * p + 9.0 / 16.0 * q * q - 0.75 * p * q * x);
}


double pi2_prime_tilde(double p, double q, double x)
{
	return sqrt(p * p + 0.25 * q * q + p * q * x);
}
