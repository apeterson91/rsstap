
	vector[q] b;
	vector[len_theta_L] theta_L;

if (t > 0) {
		if (special_case == 1) {
			int start = 1;
			theta_L = scale .* tau_b * sigma;
		if (t == 1) b = theta_L[1] * z_b;
		else for (i in 1:t) {
		int end = start + l[i] - 1;
		b[start:end] = theta_L[i] * z_b[start:end];
		start = end + 1;
	  }
	}
	else {
	  theta_L = make_theta_L(len_theta_L, p, 
							 sigma, tau_b, scale, zeta, rho, z_T);
	  b = make_b(z_b, theta_L, p, l);
	}
}
