function answer = Spec_Eng(sol, mu, xi, v, a, r)
    switch sol
        case 'xi'
            if v ~= 0 && a ~= 0
                xi = v^2/2-mu/r;
                answer = xi;
            elseif v == 0 && a == 0
                xi = -mu/(2*a);
                answer = xi;
            else
                fprintf(['Error: Solving for "xi" using Spesific ' ...
                    'Energy.\nEither "v" and "a" unknown or known'])
            end
        case 'v'
            if xi ~= 0 && a == 0
                v = (2*(xi+mu/2))^(1/2);
                answer = v;
            elseif xi == 0 && a ~=0
                v = (2*(-mu/(2*a)+mu/r))^(1/2);
                answer = v;
            else
                fprintf('Error: Solving for "v" using Spesific ' ...
                    'Energy.\nEither "xi" known and "a" unknown. Vice Versa')
            end
        case 'a'
            if xi ~= 0
                a = -mu/(2*a);
                answer = a;
            elseif xi == 0
                a = -mu/(2*(-(v^2)/2-mu/r));
                answer = a;
            else
                fprintf('')
            end
        case 'r'
    end
end
