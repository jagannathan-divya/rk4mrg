% Returns a scalar forcing at time t

function Nforce = forcing(y, t, alp, omg)
    Nforce = -alp*y + sin(omg*t);
    %Nforce = sin(omg*t);
end