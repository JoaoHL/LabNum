function [coef_matrix] = constroiv (ax, bx, nx, ay, by, ny, f, df_dx, df_dy, d2f_dxy)

	hx = (bx - ax) / nx;
	hy = (by - ay) / ny;
	coef_matrix = {};

	for a = 1:nx
		for b = 1:ny
			a00 = f(a, b);
			a01 = df_dy(a, b) * hy;
			a10 = df_dx(a, b) * hx;
			a11 = d2f_dxy(a, b) * hx * hy;

			a20 = (3*f(a+1, b)) - (df_dx(a+1, b) * hx) - (2 * a10) - (3 * a00);
			a30 = f(a+1, b) - a00 - a10 - a20;

			a02 = (3 * f(a, b+1)) - (df_dy(a, b+1) * hy) - (2 * a01) - (3 * a00);
			a03 = f(a, b+1) - a00 - a01 - a02;

			a12 = (3 * df_dx(a, b+1) * hx) - (d2f_dxy(a, b+1) * hx * hy) - (2 * a11) - (3 * a10);
			a13 = (df_dx(a, b+1) * hx) - a12 - a11 - a10;

			a21 = (3 * df_dy(a+1, b) * hy) - (d2f_dxy(a+1, b)*hx*hy) - (2 * a11) - (3 * a01);
			a31 = (df_dy(a+1, b) * hy) - a21 - a11 - a01;

			i = f(a+1, b+1) - (a00 + a01 + a02 + a03 + a10 + a11 + a12 + a13 + a20 + a21 + a30 + a31);
			j = (df_dx(a+1, b+1)*hx) - (3*a31 + 3*a30 + 2*a20 + 2*a21 + a10 + a11 + a12 + a13);
			k = (df_dy(a+1, b+1)*hy) - (a31 + a21 + 3*a13 + 2*a12 + a11 + 3*a03 + 2*a02 + a01);
			l = (d2f_dxy(a+1, b+1)*hx*hy) - (a11 + 2*a12 + 3*a13 + 2*a21 + 3*a31);

			a22 = 9*i - 3*k - 3*j + l;
			a23 = 3*k - 6*i - l - 2*j;
			a32 = 3*j - 6*i - l + 2*k;
			a33 = l - 2*k + 4*i - 2*j;
			coef_matrix(a, b) = [a00 a01 a02 a03; a10 a11 a12 a13; a20 a21 a22 a23; a30 a31 a32 a33];
		endfor
	endfor
endfunction

function [result] =  avaliav (x, y, ax, bx, nx, ay, by, ny, coef_matrix)

	hx = (bx - ax) / nx;
	hy = (by - ay) / ny;
	xi = 0;
	yi = 0;
	i = 0;
	j = 0;

	for a = 0:(nx - 1)
		xi = ax + (a * hx);

		%Achamos o intervalo [x(i), x(i+1)]
		if ((x >= xi) && (x =< xi + hx))
			i = a;
			for b = 0:(ny - 1)
				yi = ay + (b * hy);

				%Achamos o intervalo [y(i), y(i+1)]
				if((y >= yi) && (y =< yi + hy))
					j = b;
					break;
				endif
			endfor
			break;
		endif
	endfor

	% Temos que adicionar 1 em i e j pois vetores, em octave, são indexados por inteiros
	% maiores que 0, e i,j podem ser iguais a 0
	coef = coef_matrix(i+1, j+1);

	X = [1 ((x - xi) / hx) power(((x - xi) / hx), 2) power(((x - xi) / hx), 3)];
	Y = [1; ((y - yi) / hy); power(((y - yi) / hy), 2); power(((y - yi) / hy), 3)];

	result = X * coef * Y;
endfunction