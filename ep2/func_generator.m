function [] = generate_file (f_eval, df_dx_eval, df_dy_eval, d2f_dxy_eval, ax, bx, nx, ay, by, ny)
	hx = (bx - ax) / nx;
	hy = (by - ay) / ny;
	xi = 0;
	yi = 0;

	for a = 0:nx
		xi = ax + (a * hx);
		for b = 0:ny
			yi = ay + (b * hy);

			f(a+1, b+1) = f_eval(xi, yi);
			df_dx(a+1, b+1) = df_dx_eval(xi, yi);
			df_dy(a+1, b+1) = df_dy_eval(xi, yi);
			d2f_dxy(a+1, b+1) = d2f_dxy_eval(xi, yi);
		endfor
	endfor

	save function_parameters.m ax bx nx ay by ny f df_dx df_dy d2f_dxy;
endfunction

function [result] = f_eval (x, y)
	result = sin(x * cos(y));
endfunction

function [result] = df_dx_eval (x, y)
	result = cos(y) * cos(x * cos(y));
endfunction

function [result] = df_dy_eval (x, y)
	result = (-x)*sin(y)*cos(x*cos(y));
endfunction

function [result] = d2f_dxy_eval (x, y)
	result = sin(y)*(x*cos(y)*sin(x*cos(y)) - cos(x*cos(y)));
endfunction

args = argv();
ax = str2double(args(1,1));
bx = str2double(args(2,1));
nx = str2double(args(3,1));
ay = str2double(args(4,1));
by = str2double(args(5,1));
ny = str2double(args(6,1));

generate_file(@f_eval, @df_dx_eval, @df_dy_eval, @d2f_dxy_eval, ax, bx, nx, ay, by, ny);
