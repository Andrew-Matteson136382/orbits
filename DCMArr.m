function [DCM] = DCMArr(i, O, Th)
    ci = cos(i);
    cO = cos(O);
    cTh = cos(Th);
    si = sin(i);
    sO = sin(O);
    sTh = sin(Th);
    DCM = [cO*cTh-sO*ci*sTh, -sTh*cO-sO*ci*cTh, sO*si;
        sO*cTh+cTh*ci*sO, -sO*sTh+cO*ci*cTh, -cTh*si;
        si*sTh, si*cTh, ci];
end