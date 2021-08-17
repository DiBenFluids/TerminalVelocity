function [terminal_velocity_sphere,terminal_velocity_shape,C_D,diameter] = terminal_velocity_master(varargin)

% Terminal settling velocity master code
% wrriten by Michelle DiBenedetto, last updated 8/16/2021


% All inputs are Name - Value pairs: (see examples if confused)
% All inputs are optional - if not enough inputs given, it will give an error
% Outputs: terminal vel of equivalently sized sphere, terminal velocity
% using shape factor, drag coefficients for both, and diameter of
% equivalently sized sphere
% Default units are metric with m, kg, s

% Particle inputs:
% diameter     - effective spherical diameter
% majoraxis    - diameter not radius
% minoraxis    - diameter not radius
% mediumaxis   - diameter not radius
% rho_p        - particle density
% volume
% mass
% aspect ratio


% Fluid inputs:
% rho_f        - fluid density
% mu           - dynamic viscosity
% nu           - kinematic viscosity

% Methods

% stokes
% SN           - Schiller-Nauman
% rubey
% newtonian


% Needed
% more shape factors and better reynolds number shape coefficients
% also could do orientation dependent shape velocities
% Papers to potentially incoporate:
% Ouchene e tal 2016
% Dioguardi et al 2018

%% Input parsing

% Default Values
defaultVolume = [];
defaultMass   = [];
defaultD      = [];
defaultMaj    = [];
defaultMin    = [];
defaultMed    = [];
defaultRho_p  = [];
defaultAR     = [];
defaultPer    = [];
defaultArea   = [];
defaultAngle  = [];

defaultOrient = 'max';
expectedOrients = {'isotropic','max','min'};

defaultRho_f  = 1000;
defaultMu     = .001;
defaultNu     = 10^-6;
defaultG      = 9.81;

defaultMethod = 'SN';
expectedMethods = {'stokes','SN','rubey','newtonian','dioguardi'};

% Parse
p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);

addParameter(p,'volume',defaultVolume,validScalarPosNum);
addParameter(p,'mass',defaultMass,validScalarPosNum);
addParameter(p,'diameter',defaultD,validScalarPosNum);
addParameter(p,'majoraxis',defaultMaj,validScalarPosNum);
addParameter(p,'minoraxis',defaultMin,validScalarPosNum);
addParameter(p,'mediumaxis',defaultMed,validScalarPosNum);
addParameter(p,'rho_p',defaultRho_p,validScalarPosNum);
addParameter(p,'aspectratio',defaultAR,validScalarPosNum);
addParameter(p,'phi',defaultAngle,validScalarPosNum);

addParameter(p,'rho_f',defaultRho_f,validScalarPosNum);
addParameter(p,'mu',defaultMu,validScalarPosNum);
addParameter(p,'nu',defaultNu,validScalarPosNum);
addParameter(p,'g',defaultG,validScalarPosNum);

addOptional(p,'Perimeter',defaultPer,validScalarPosNum);
addOptional(p,'ProjectedArea',defaultArea,validScalarPosNum);

% oblate or prolate
defaultSph = 'prolate';
addOptional(p,'spheroid',defaultSph,...
    @(x) any(validatestring(x,{'prolate','oblate'})));
addParameter(p,'method',defaultMethod,...
    @(x) any(validatestring(x,expectedMethods)));
addParameter(p,'orientation',defaultOrient,...
    @(x) any(validatestring(x,expectedOrients)));


parse(p,varargin{:});

% Convert structure to variables
fn = fieldnames(p.Results);
for i=1:length(fn)
    varname=fn{i};
    eval([varname '=p.Results.' varname ';']) 
end

%% Configure any other missing variables

if isempty(aspectratio) && ~isempty(majoraxis) && ~isempty(minoraxis)
    switch spheroid
        case 'prolate'
            aspectratio = majoraxis/minoraxis;
        case 'oblate'
            aspectratio = minoraxis/majoraxis;
    end
elseif isempty(aspectratio)
    aspectratio=1;
end


%% Find effective diameter of particle if not given

if isempty(diameter)
    
    % Option 0: have volume
    if ~isempty(volume)
        r = ((volume*(3/4/pi)).^(1/3));
        diameter = r*2;
        
        % Option 1: have mass, density, no effective diameter
    elseif ~isempty(mass) && ~isempty(rho_p)
        volume = mass / rho_p;
        r = ((volume*(3/4/pi)).^(1/3));
        diameter = r*2;
        
        % Option 2: have major and minor and medium axis lengths, assuming ellipsoid
    elseif  ~isempty(majoraxis) && ~isempty(minoraxis) && ~isempty(mediumaxis)
        volume = 4/3 * majoraxis/2 * minoraxis/2 * mediumaxis/2;
        r = ((volume*(3/4/pi)).^(1/3));
        diameter = r*2;
        
        % Option 3: have major and minor axis lengths, assuming spheroid
    elseif  ~isempty(majoraxis) && ~isempty(minoraxis) && ~isempty(aspectratio)
        
        if aspectratio > 1
            % prolate
            volume = 4/3 * majoraxis/2 * minoraxis/2 * minoraxis/2;
            r = ((volume*(3/4/pi)).^(1/3));
            diameter = r*2;
        else
            %oblate
            volume = 4/3 * majoraxis/2 * majoraxis/2 * minoraxis/2;
            r = ((volume*(3/4/pi)).^(1/3));
            diameter = r*2;
        end
    else
        warning('can''t calculate effective diameter, not enough inputs')
        
    end
    
else
    r = diameter/2;
    if isempty(volume)
        volume = 4/3 * pi * r^3;
    end
    
    if isempty(rho_p)
        if ~isempty(mass)
            rho_p = mass / volume;
        else
            warning('can''t calculate particle density, not enough inputs')
        end
    end
end

if rho_p < rho_f
    buoyant = 1;
else
    buoyant = -1;
end

% nonspherical with major and minor axes, but no medium
if isempty(mediumaxis) && ~isempty(majoraxis) && ~isempty(minoraxis)
    mediumaxis= ( volume*3/4/pi/(majoraxis/2) /(minoraxis/2) )/2;
end

% nonspherical with only AR and no axes definied: assumed spheroid
 % fill in later


%% Estimate terminal velocities: effective spherical

Area=(diameter/2)^2*pi;  %Frontal area

% Initial guess
% from Table 5.3 in Clift, Grace and Webber 2005

Nd = abs((4/3)*(1/nu^2)*(rho_p/rho_f - 1)*g*(diameter^3));

if Nd <= 73
    Re_CGW = Nd/24 - 1.7569e-4*Nd^2 + 6.9252e-7*Nd^3 - 2.3027e-10*Nd^4;
elseif (73 < Nd) && (Nd <= 580)
    Re_CGW = 10^(-1.7095 + 1.33438*log10(Nd) - 0.11591*(log10(Nd)^2));
elseif (580 < Nd) && (Nd <= 1.55e7)
    Re_CGW = 10^(-1.81391 + 1.34671*log10(Nd) - 0.12427*(log10(Nd)^2) + 0.006344*(log10(Nd)^3));
elseif (1.55e7 < Nd) && (Nd <= 5e10)
    Re_CGW = 10^(5.33283 - 1.21728*log10(Nd) + 0.19007*(log10(Nd)^2) - 0.007005*(log10(Nd)^3));
end

guess = Re_CGW*nu/diameter;


switch method
    
    %F_drag=1/2 rho *Cd *U^2*Area (downward positive)
    %F_buoyancy= -rho_f*volume*g  (upward negative)
    %F_gravity = rho_p*volume*g   (downward positive)
    
    
    case 'rubey'
        % Cd= (24/(U*diameter/nu)+.44) (Rubey, 1933)
        terminal_velocity_sphere = fzero(@(x) (24/abs(x*diameter/nu)+.44) * 1/2 * rho_f * x * abs(x) * Area - rho_f * volume * g + rho_p * volume * g,guess);
        
        Re = abs(terminal_velocity_sphere*diameter/nu);
        C_D = (24/Re+.44);
        
        
    case 'stokes'
        
        terminal_velocity_sphere = - 2/9 * (rho_p - rho_f) / mu * g * r^2;
        
        Re = abs(terminal_velocity_sphere*diameter/nu);
        C_D = 24/Re;
        
    case 'SN'
        
        % Cd = 24/R (1 + 0.15 * R ^ .687 ) (Schiller && Naumann)
        terminal_velocity_sphere = fzero(@(x)  24/abs(x*diameter/nu)*(1+0.15*abs(x*diameter/nu)^.687) * 1/2 * rho_f * x * abs(x) * Area - rho_f * volume * g + rho_p * volume * g,guess);
        
        Re = abs(terminal_velocity_sphere*diameter/nu);
        C_D = 24/Re*(1+0.15*Re^.687);
        
    case 'newtonian'
        
        % Cd = 0.44
        terminal_velocity_sphere = buoyant * 2.46 * (abs(((rho_p - rho_f) * g * r)/rho_f))^.5;
        
        C_D= 0.44;
        
    otherwise %Default is SN
       % disp('Using Schiller and Naumann for equivalent spherical settling')
        % Cd = 24/R (1 + 0.15 * R ^ .687 ) (Schiller && Naumann)
        terminal_velocity_sphere = fzero(@(x)  24/abs(x*diameter/nu)*(1+0.15*abs(x*diameter/nu)^.687) * 1/2 * rho_f * x * abs(x) * Area - rho_f * volume * g + rho_p * volume * g,guess);
        
        Re = abs(terminal_velocity_sphere*diameter/nu);
        C_D = 24/Re*(1+0.15*Re^.687);
end

%% Shape effects
terminal_velocity_shape=terminal_velocity_sphere;

if aspectratio~=1
    
    %Find shape factor (only exact at small Reynolds numbers)
    [f_par,f_per] = ellipsoid_tensor_coeff(aspectratio);
    
    
    % could improve by using reynolds number dependent shape correction-
    % wow they are all crazy though

    % K is shape factor
    if ~isempty(phi)
            K = 1/2*(f_par+f_per+(f_par-f_per)*cos(2*angle));
    else
        switch orientation
            case 'max'
                    K = max(f_par,f_per);

            case 'min'
                    K = min(f_par,f_per);

            case 'isotropic'
                    K = (f_par+ 2*f_per)/3;
        end
    end
    
    

    
    switch method
        
        case 'rubey'
            terminal_velocity_shape = fzero(@(x) K * (24/abs(x*diameter/nu)+.44) * 1/2*rho_f*x * abs(x)*Area - rho_f*volume*g + rho_p*volume*g,guess);
            
            Re = abs(terminal_velocity_shape*diameter/nu);
            C_D(2) = (24/Re+.44) * K;
            
        case 'stokes'
            
            terminal_velocity_shape = terminal_velocity_sphere / K;
            
            Re = abs(terminal_velocity_shape*diameter/nu);
            C_D(2) = 24/Re;
            
        case 'SN'
            
            terminal_velocity_shape = fzero(@(x) K * 24/abs(x*diameter/nu)*(1+0.15*abs(x*diameter/nu)^.687) * 1/2 * rho_f * x * abs(x) * Area - rho_f * volume * g + rho_p * volume * g,guess);
            
            Re = abs(terminal_velocity_shape*diameter/nu);
            C_D(2) = 24/Re*(1+0.15*Re^.687);
            
            
        case 'newtonian'
            
            terminal_velocity_shape = terminal_velocity_sphere / K;
            
            C_D(2)= 0.44 * K;
            
        case 'dioguardi'
            
            %Dioguardi et al., 2017; Dellino et al., 2005:
            %Shape factor = sphericity/ circularity
            
            % needs inputs: Perimeter and ProjectedArea and all 3 axes
            % should be checked
            
            % Sphericity
            a = majoraxis/2;
            b = mediumaxis/2;
            c = minoraxis/2;
            z = 1.6075; % Dellino et al., 2005
            Surface_Area_shape = 4*pi *(( (a*b)^z + (a*c)^z + (c*b)^z )/3)^(1/z);    % assume ellipsoid, use approximate exp
            Surface_Area_sphere = 4*pi*r^2;
            Phi = Surface_Area_sphere/Surface_Area_shape;
            
            % Circularity
            
            r_projected = sqrt(ProjectedArea/pi);
            Perimeter_projected = 2*pi*r_projected;
            X = Perimeter / Perimeter_projected;
            
            %Shape Factor
            Psi = Phi/X;
            
            terminal_velocity_shape = fzero(@(x)  ( 24/abs(x*diameter/nu)* ((1-Psi)/abs(x*diameter/nu)+1)^.25...
                + 24/abs(x*diameter/nu)*(.01806*abs(x*diameter/nu)^.6459)*Psi^(-abs(x*diameter/nu)^.08)...
                + 0.4251/ ( 1 + 6880.95/abs(x*diameter/nu)*Psi^5.05)  )... 
                * 1/2 * rho_f * x * abs(x) * Area - rho_f * volume * g + rho_p * volume * g,guess);
    
            Re = abs(terminal_velocity_shape*diameter/nu);
            C_D(2) = ( 24/Re* ((1-Psi)/Re+1)^.25...
                + 24/Re*(.01806*Re^.6459)*Psi^(-Re^.08)...
                + 0.4251/ ( 1 + 6880.95/Re*Psi^5.05));
            
            
    end
    
    
end
end


function [f_par,f_per] = ellipsoid_tensor_coeff(w)
%resistance tensor, non-dimensional, for spheroidal particle of aspect
%ratio w, equations from Loth 2008 (based on previous literature)

if w<1 %oblate exact
    S=acos(w)/sqrt(1-w^2);
    f_par=4/3*w^(-1/3)*(1-w^2)/(w+(1-2*w^2)*S);
    f_per=8/3*w^(-1/3)*(w^2-1)/(w-(3-2*w^2)*S);
    
   % oblate (Happel & Brenner section 5-11)
   % R_par = (8/3) * w^(-1/3) * (2*w/(1-w^2) + 2*(1-2*w^2)/(1-w^2)^(3/2)*atan(sqrt(1-w^2)/w) )^-1;
   % R_per = (8/3) * w^(-1/3) * (-w/(1-w^2) - (2*w^2-3)/(1-w^2)^(3/2)*asin(sqrt(1-w^2) ))^-1;
    
elseif w>1% && w<=6
    %prolate approx
    % f_par= (4/5+w/5)*w^(-1/3);
    % f_per=(3/5+2*w/5)*w^(-1/3);
    
    %prolate exact (Loth 2008) table 1
    f_par=(4/3)*w^(-1/3)*(1-w^2)/(w-(2*w^2-1)*log(w+sqrt(w^2-1))/sqrt(w^2-1));
    f_per=(8/3)*w^(-1/3)*(w^2-1)/(w+(2*w^2-3)*log(w+sqrt(w^2-1))/sqrt(w^2-1));
    
    % elseif w>6 %needle
    % f_par=2/3*w^(2/3)/(log(2*w)-.5);
    % f_per=2*f_par;
    
    % prolate (Happel & Brenner section 5-11)
    %R_par = (8/3) * w^(-1/3) * (-2*w/(w^2-1)+(2*w^2-1)/(w^2-1)^(3/2)*log((w+sqrt(w^2-1))/(w-sqrt(w^2-1))))^-1;
    %R_per = (8/3) * w^(-1/3) * (w/(w^2-1)+(2*w^2-3)/(w^2-1)^(3/2)*log((w+sqrt(w^2-1))))^-1;
end

end







