% My function. Enter your arbitrary function here.
clear all
f = @(x) (20*x ^ 4 - 22 * x ^ 3 + 2);
df = @(x) (80*x^3 - 66*x^2);

x1 = 0.5;
x2 = 1;
tol = 10^-6;

% Invoking function
cubic_interpolation(f, df, x1, x2, tol)

function cubic_interpolation(f, df, x1, x2, tol)
    if ( x1 < x2 && df(x1) < 0 && df(x2) > 0 )
        cond = true;
        iterations = 1;
        while(cond)
            syms a b c d
            fx1 = f(x1);
            fx2 = f(x2);
            dfx1 = df(x1);
            dfx2 = df(x2);
            
            % Cubic interpolation in the form of a third-order polynomial y(x) = a + bx + cx^2 + dx^3
            % Coefficients a, b, c, d are calculated from the following linear equations
            eq1 = a + x1*b + x1^2 * c + x1^3 * d == fx1;
            eq2 = a + x2*b + x2^2 * c + x2^3 * d == fx2;
            eq3 = b + 2*c*x1 + 3*d*x1^2 == dfx1;
            eq4 = b + 2*c*x2 + 3*d*x2^2 == dfx2;
            sol = solve([eq1, eq2, eq3, eq4], [a, b, c, d]);
            aS = sol.a;
            bS = sol.b;
            cS = sol.c;
            dS = sol.d;
            
            syms x
            cubicFunc = aS + bS*x + cS*x^2 + dS*x^3;
            dfCubicFunc = diff(cubicFunc, x);
            dfCubicFuncRoots = root(dfCubicFunc);
            
            d2fCubicFunc = matlabFunction(diff(dfCubicFunc, x));
            x1kappa = double(dfCubicFuncRoots(1,1));
            x2kappa = double(dfCubicFuncRoots(2,1));
            
            % One root is the minimum, one is the maximum.
            % Simple check to determine which one is the minimum.
            if (double(d2fCubicFunc(x1kappa))> 0)
                xOpt = x1kappa;
            else
                xOpt =x2kappa;
            end
            
            fxOpt = f(xOpt);
            cubicToMF = matlabFunction(cubicFunc);
            cubeOpt = cubicToMF(xOpt);
            
            cond = (double(abs(fxOpt - cubeOpt)) > tol);
                if (cond)
                    dfOptToMF = matlabFunction(dfCubicFunc);
                    dfOpt = dfOptToMF(xOpt);
                    iterations = iterations + 1;
                    if (double(dfOpt) > 0)
                        x1 = xOpt;
                    else
                        x2 = xOpt;
                    end
                end
        end
            fprintf('Optimal x* = %.4f \n', xOpt)
            fprintf('No. of iterations = %d \n', iterations)
    else
        disp('x1 < x2 && df(x1) < 0 && df(x2) > 0 not satisfied, aborting method...')
    end
end