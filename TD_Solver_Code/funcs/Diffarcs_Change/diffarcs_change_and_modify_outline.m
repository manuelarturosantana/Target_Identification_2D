% To the curves class I added this flipped parabola, which uses ap as a
% parameter to move to move this flipped kite around.
% Here is an outline if we decided to publish the code, but do not want to
% publish all of diffarcs
% We would need curve parameterization
% and BuildCurve (from which we could take out the open curves :))
% Then we need Build Single Matrix
% and Single Extern
% And that is it, that is all of Diffarcs that we are using :)
% 


% 
if(strcmp(Curve_Flag,'parabola flipped'))
    % This is a little hack, I'm using ap to move the parabola :)
    x=-(1-2*s.^2) + ap;
    y=s;
    xp=4*s;
    yp=-1+0*s;
    xpp=4+0.*s;
    ypp=0*s;

end;