
void mean_variance(double* measure, int N, double* results, bool print) {// compute the mean and variance of an array
	double s = 0, s2 = 0;
	for (int i = 0; i < N; i++) {
		s += measure[i];
		s2 += measure[i] * measure[i];
	}
	results[0] = s / N;
	results[1] = (s2 / N - results[0] * results[0]);
	if (print) printf("%.10lf +- %.10lf\n", results[0], sqrt(results[1]/(N-1)));
}

void autocorrelation(double* measure, int N, int W, double mean,double variance,double* result) {//compute autocorrelation of an array
	double s = 0;
	for (int i=1; i < W+1; i++) {
		for (int j = 0; j < (N - i); j++)
			s += (measure[j] - mean) * (measure[j + i] - mean);

		result[i-1] = s / (N - i) / variance;
		if (i<11) printf("%d %lf\n", i, result[i-1]);
		s = 0;
	}
}

