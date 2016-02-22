function f = mz2f(mz, P, Instrument) %Update for other machines
%converts from mz, returns frequency
%P.A, P.B and P.C contain calibration parameters
%if C is supplied, it is expected to apply to the form of calibration equation derived by Masselon et al. 2002, ref 63:
% m/z = A./f + B./f.^2 + C*I./f.^2, where I is the intensity of each peak

switch(Instrument)
    case{'ltqft'}
        
%         numTerms = 2;   % number of parameters in calibration equation
%         
%         if nargin == 1
%             % no P parameter supplied
%             disp('WARNING: Insufficient parameters specified, using nominal values');
%             numTerms = 2;
%             A = 107368.5e3;
%             B = -750.461e6;
%         elseif isempty(P)
%             % P is empty: use defaults
%             numTerms = 2;
%             A = 107368.5e3;
%             B = -750.461e6;
%         elseif ~isfield(P,'C') || P.C==0
%             % no (or zero) C parameter: use two-term calibration
            numTerms = 2;
            A = P.A;
            B = P.B;
            C = 0;
%         else
%             % 3-term calibration
%             numTerms = 3;
%             A = P.A;
%             B = P.B;
%             C = P.C;
%         end
        
%         switch numTerms
%             case 2
                f = (A + sqrt( A^2+mz.*4*B ))./(2.*mz);
%             case 3
%                 %solve using symbolic
%                 syms As Bs Cs Is Fs mzs
%                 mz = -mzs + As./fs + Bs./fs.^2 + Cs*ys./fs.^2;
%                 f = solve(mz,fs);
%                 mzL = length(mz);
%                 if mzL > 1
%                     A = A.*ones(1,mzL);
%                     B = B.*ones(1,mzL);
%                     C = C.*ones(1,mzL);
%                 end
%                 f = (subs(f,{As,Bs,Cs,ys,mzs},{A,B,C,y,mz}));
%                 f = max(real(f));
%         end
        
        case{'qexactive','orbitrap'}
            B = P.B;
            C = P.C;
            f = sqrt((B+sqrt(B^2+mz.*4*C))./(2.*mz)); % Report frequency in Hz.
end
end

