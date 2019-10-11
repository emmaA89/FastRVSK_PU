%-------------------------------------------------------------------------%
%
% File: Scale_Function(dsites,h,puctrs)
%
% Goal: maps points via the scale function
%
% Inputs: dsites:   data points
%         h:        parameter for the scale function
%         puctrs:   the PU centre
%
% Outputs:  val: the data values mapped via the scale function
%
% Authors: S. De Marchi, A. Martínez, E. Perracchione,
%          Universita' di Padova,
%          Dipartimento di Matematica "Tullio Levi-Civita".
%
% Last modified: 10/10/17.
%
%-------------------------------------------------------------------------%
function val = Scale_Function(dsites,h,puctrs)
val =  0.5+sqrt(h-((dsites(:,1)-puctrs(1)).^2+...
    (dsites(:,2)-puctrs(2)).^2));
val = val';
end