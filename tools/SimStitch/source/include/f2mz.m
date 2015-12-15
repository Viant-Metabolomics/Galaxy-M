function mz = f2mz(f, P, Instrument, y) %Update for other machines
%converts from freq, returns m/z
%P.A, P.B and P.C contain calibration parameters
%if C is supplied, it is expected to apply to the form of calibration equation derived by Masselon et al. 2002, ref 63:
% m/z = A./f + B./f.^2 + C*y./f.^2, where y is the intensity of each peak

switch(Instrument)
    case{'ltqft','orbitrap'}
        %numTerms = 2;   % number of parameters in calibration equation
        
%         if nargin == 1
%             % no P parameter supplied
%             disp('WARNING: Insufficient parameters specified, using nominal values');
%             numTerms = 2;
%             A = 107368.5e3;
%             B = -750.461e6;
%         elseif isempty(P)
%             % P is empty: use defaults without warning
%             numTerms = 2;
%             A = 107368.5e3;
%             B = -750.461e6;
%         elseif ~isfield(P,'C') || P.C==0
%             % no or zero C (therefore y irrelevant) parameter: use two-term calibration
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
                mz = A./f + B./f.^2;
%             case 3
%                 mz = A./f + B./f.^2 + C*y./f.^2;
%         end
        
    case{'qexactive'}
%         if nargin == 1
%             % no P parameter supplied
%             disp('WARNING: Insufficient parameters specified, using nominal values');
%             B = 67862927.5887*1e6;
%             C = 26143009.9533*1e12;
%         elseif isempty(P)
%             % P is empty: use defaults without warning
%             B = 67862927.5887*1e6;
%             C = 26143009.9533*1e12;
%         else
            B = P.B;
            C = P.C;
%         end
        mz = B./(f.^2)+C./(f.^4);
end
end

