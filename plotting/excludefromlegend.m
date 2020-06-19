function excludefromlegend(ll)

%% MCB, GDFL, 2012-07-02
%% ll is line handle (vector) to be excluded from legend
for i=1:length(ll)
    set(get(get(ll(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
end

return
