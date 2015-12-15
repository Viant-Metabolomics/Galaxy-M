function [P] = Calibrate_1_0(SP, P, calmz, mode, weighted, Instrument)
%Determines new calibration parameters given observed peak frequencies (SP) and a list of
%known masses ie cal points (calmz) in m/z, and initial A and B parameters

% INPUTS
% SP is a structure with fields .y and .f corresponding to the peaks in the
%   spectrum identified as calibrants, SP sorted by increasing f.
% P is the parameters for SP
% calmz is the calibrant exact masses, must correspond to order of peaks in
%   SP
% mode is the calibration mode (b or ab or abc)
% weighted is 1: weight calibration to more intense peaks, 0: no weighting

% OUTPUTS
% P is updated with the new (recalibrated) calibration parameters

% UPDATES
% 18/6/07 (v0.1): FindClosest now v0.1
% 19/6/07       : Max peak distance now parameter
% 12/10/07      : Now also returns "data", the abundance of the calibrant
%                   in the spectrum
% 26/02/08 (v1.0):  Re-structured, removed output, input peaks are now
%                   already identified with calibrants
%
%
% Needs to be update with respect to the other calibration equations of the
% QExactive and Bruker solarix! 

%check calmz in single column
if (size(calmz,2)>size(calmz,1)), calmz=calmz.'; end

%weighting
if weighted, calmzw = calmz.*(SP.y).'; end

%construct F for mz = A/f + B/f^2
f = SP.f.';
y = SP.y.';
A = P.A;
B = P.B;
C = P.C;

switch(Instrument)
    case{'ltqft','orbitrap'}
        switch mode
            case 'b'
                F = 1./f.^2;
                if weighted
                    F = F.*y;
                    Xh = F\(calmzw - A*1./f.*y);
                else
                    Xh = F\(calmz - A*1./f);
                end
                newA = A;
                newB = Xh;
                if C, error('Non-zero C value in b-mode recalibration?'); end
                newC = 0;
            case 'ab'
                F = [1./f 1./f.^2];
                if weighted
                    F = F.*[y y];
                    Xh = F\calmzw;
                else
                    Xh = F\calmz;
                end
                newA = Xh(1);
                newB = Xh(2);
                if C, error('Non-zero C value in ab-mode recalibration?'); end
                newC = 0;
            case 'abc'  %see ref 63, Masselon et al., 2002
                F = [1./f 1./f.^2 y./f.^2];
                if weighted
                    F = F.*[y y y];
                    Xh = F\calmzw;
                else
                    Xh = F\calmz;
                end
                newA = Xh(1);
                newB = Xh(2);
                newC = Xh(3);
        end
    case{'qexactive'}
        switch mode
            case 'b'
                F = 1./f.^4;
                if weighted
                    F = F.*y;
                    Xh = F\(calmzw - B*(1./f.^2).*y);                   
                else
                    Xh = F\(calmz - B*(1./f.^2));
                end
                newA = 0;
                newB = B;
                newC = Xh;
            case 'ab'
                %f = f./10^3; % Equation based on parameter in kHz
                F = [1./f.^2 1./f.^4];
                if weighted
                    F = F.*[y y];
                    Xh = F\calmzw;
                else
                    Xh = F\calmz;
                end
                newB = Xh(1);
                newC = Xh(2);
                newA = 0;
        end
end

P.A = newA;
P.B = newB;
P.C = newC;
end

