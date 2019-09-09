function [y,x,dateselected,Rh]=datageneration_midas(rv,options,date,R)
% Data transformation function for MiDaS regressions, non-overlapping
% Takes the sum of RV over the horizon h
% It the daily simple R series is also passed, it computes the cumulative
% h-day simple return over the same period over which MIDAS is computed

agy=options.aggrY;
agx=options.aggrX;
ssize=length(rv);

maxy=floor((ssize-agx)/agy);

% LHS
y=sum((reshape(rv((end-maxy*agy+1):end),agy,maxy)'),2); % <=== sum of realized RV over that horizon agy
% RHS
indy=(ssize-maxy*agy):agy:(ssize-agy);
x=NaN(length(indy),agx);
for i=1:length(indy)
    x(i,:)=(rv(indy(i):-1:indy(i)-agx+1))';
end
% dates
if nargin>=3
    begin = size(rv,1)-maxy*agy+1;
    indd=(begin:agy:size(rv,1)-agy+1)';
    dateselected = date(indd,:); % this gets the first date of the non-overlapping period
end
% returns
if nargin == 4
    Rh = sum((reshape(R((end-maxy*agy+1):end),agy,maxy)'),2); % <=== sum of daily log returns over that horizon agy
end

end









