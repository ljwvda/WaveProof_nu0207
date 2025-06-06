function N=sizeshape_Wave(shape)
% determines the tightest tensor that contains the shape, adjusted for
% traveling waves

switch shape.type
    case {'rec'} % rectangle 
         N=shape.Nrec;
    case {'ell'} % for Edagger (for Adagger) and Etilde
         tN=shape.Nell; 
         nu=shape.nu;
         c=shape.omega;
         N(1:3)=ceil(altsup(sqrt(tN/nu)));
         % N(3) can be a bit tighter, as we need to solve |nu n3^2-icn3|<=N for n3
         %(see notes 2-12-24)
         % N(1:2)=ceil(altsup(sqrt(tN/nu)));
         % N(3) = ceil(altsup(sqrt((-c^2+sqrt(c^4+4*nu^2*tN^2))/(2*nu^2))));
         %%N(3) = ceil(altsup(sqrt((-c^2-sqrt(c^4+4*nu^2*tN^2))/(2*nu^2))));
         N(4)=0; % Use 0 to save memory
    case {'ell2D'} % ellipse for 2D, not used for waves
         tN=shape.Nell; 
         nu=shape.nu;
         c=shape.omega;
         N(1:2)=ceil(altsup(sqrt(tN/nu)));
         N(3)=0; % this saves quite a bit of memory
         % N(4)=ceil(altsup(tN/abs(c)));        %Adjusted
         N(4)=0;% Use 0 to save memory
    case {'otherplusrec'} % sum of the sets: rectangle + other shape
         N=shape.Nrec+sizeshape_Wave(shape.other); % Adjusted
end

