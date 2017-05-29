function [] = main (filename, interactive, x = 0, y = 0)
    load (filename);
    coef_matrix = constroiv(ax, bx, nx, ay, by, ny, f);
    
    if (interactive)
        while true
            x = input('Digite x: ');
            y = input('Digite y: ');
            result = avaliav(x, y, ax, bx, nx, ay, by, ny, coef_matrix);
            printf('Resultado: %f\n', result);
        endwhile
    else
        result = avaliav(x, y, ax, bx, nx, ay, by, ny, coef_matrix);
        printf('Resultado: %f\n', result);
    endif
endfunction

function [coef_matrix] = constroiv(ax, bx, nx, ay, by, ny, f)

    hx = (bx - ax) / nx;
    hy = (by - ay) / ny;
    x  = [];
    y  = [];
    coef_matrix  = zeros((nx - 1), (ny - 1), 4);

    for (i = 0: nx)
        # calculando os xi
        x(i+1) = ax + (i * hx);
    endfor
    
    for (j = 0: ny)
        # calculando os yi
        y(j+1) = ay + (j * hy);
    endfor

    for (i = 1: nx - 1)
        for (j = 1: ny - 1)
            # das deduções encontramos os coeficientes
            a00 = f(i, j);
            a01 = (f(i, j+1) - a00) * (hy / (y(j+1) - y(j)));
            a10 = (f(i+1, j) - a00) * (hx / (x(i+1) - x(i)));
            a11 = (f(i+1, j+1) - f(i+1, j) - f(i, j+1) + a00) * (hx * hy) / ((x(i+1) - x(i)) * (y(j+1) * (j)));
            coef_matrix(i, j, :) = [a00, a01, a10, a11];
        endfor
    endfor
endfunction

function result = avaliav(x, y, ax, bx, nx, ay, by, ny, coef_matrix);

    hx = (bx - ax) / nx;
    hy = (by - ay) / ny;
    xa = 0;
    yb = 0;  

    for (i = 0: nx - 1)
        xa = ax + (i * hx);
        for (j = 0: ny - 1)
            yb = ay + (j * hy);
            # caso (x, y) esteja no intervalo especificado
            if ((x >= xa) && (x <= (xa + hx)) && (y >= yb) && (y <= (yb + hy)))
                result = coef_matrix(i+1, j+1, 1) +  (coef_matrix(i+1, j+1, 2) * (x - xa)) / hx + (coef_matrix(i+1, j+1, 3) * (y - yb)) / hy + (coef_matrix(i+1, j+1, 4)) * (x - xa) * (y - yb) / (hx * hy); 
            endif
        endfor
    endfor
endfunction

args = argv();
interactive = false;

if (length(args) < 3)
    filename = args{1, 1};
    if (args{2, 1} == '-i')
        interactive = true;
    endif
    main(filename, interactive);
else
    filename = args{1, 1};
    x = str2double(args{2, 1});
    y = str2double(args{3, 1});
    main(filename, interactive, x, y);
endif
