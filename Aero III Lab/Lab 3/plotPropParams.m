function figcount = plotPropParams(J, etap_plot, C_P, C_T, C_tau, pinds, figcount)
% y is the curve want to plot (etap_plot,C_P,C_T,C_tau)
% x is the data want to plot against (J)
% pinds is the propellers to plot
% outputs new figcount

%    1 to 3 for pitch (all 4 curves coplotted)
%    3 to 6 for diameter (all 4 curves coplotted)
%    3 and 7 for blade number, 7 has larger dia by 0.5 but thats the


legStart = 'Prop '; % start of each legend entry
pindsL = length(pinds);
%ind = pitch(p); % define what index of propellers to plot using parameter vec
bound = pinds;
figure(figcount) % etap Plots
for i = bound % num of curves to plot
    plot(J{i}, etap_plot{i},'v')
    hold on
end
xlabel('J')
ylabel('\eta_p')
% Legend
if pindsL == 3 % pitch case
    title('Varying pitch')
    legend([legStart,sprintf('%d',pinds(1))], [legStart,sprintf('%d',pinds(2))], [legStart,sprintf('%d',pinds(3))]); % need to automate this
elseif pindsL == 2 % blade num case
    title('Varying Blade Number')
    legend([legStart,sprintf('%d',pinds(1))], [legStart,sprintf('%d',pinds(2))]);
else % blade dia case where there are 4
    title('Varying Blade Diameter')
    legend([legStart,sprintf('%d',pinds(1))], [legStart,sprintf('%d',pinds(2))], [legStart,sprintf('%d',pinds(3))], [legStart,sprintf('%d',pinds(4))]);
end

figure(figcount+1) % C_T plots
for i = bound % num of curves to plot
    plot(J{i}, C_T{i},'o')
    hold on
end
xlabel('J')
ylabel('C_T')
% Legend
if pindsL == 3 % pitch case
    title('Varying pitch')
    legend([legStart,sprintf('%d',pinds(1))], [legStart,sprintf('%d',pinds(2))], [legStart,sprintf('%d',pinds(3))]); % need to automate this
elseif pindsL == 2 % blade num case
    title('Varying Blade Number')
    legend([legStart,sprintf('%d',pinds(1))], [legStart,sprintf('%d',pinds(2))]);
else % blade dia case where there are 4
    title('Varying Blade Diameter')
    legend([legStart,sprintf('%d',pinds(1))], [legStart,sprintf('%d',pinds(2))], [legStart,sprintf('%d',pinds(3))], [legStart,sprintf('%d',pinds(4))]);
end

figure(figcount+2) % C_tau plots
for i = bound
    plot(J{i}, C_tau{i}.*10, '.','markersize',10) % multiplied by 10
    hold on
end
xlabel('J')
ylabel('C_\tau (x10)')
% Legend
if pindsL == 3 % pitch case
    title('Varying pitch')
    legend([legStart,sprintf('%d',pinds(1))], [legStart,sprintf('%d',pinds(2))], [legStart,sprintf('%d',pinds(3))]); % need to automate this
elseif pindsL == 2 % blade num case
    title('Varying Blade Number')
    legend([legStart,sprintf('%d',pinds(1))], [legStart,sprintf('%d',pinds(2))]);
else % blade dia case where there are 4
    title('Varying Blade Diameter')
    legend([legStart,sprintf('%d',pinds(1))], [legStart,sprintf('%d',pinds(2))], [legStart,sprintf('%d',pinds(3))], [legStart,sprintf('%d',pinds(4))]);
end

figure(figcount+3) % C_P plots
for i = bound
    plot(J{i}, C_P{i}, '*')
    hold on
end
xlabel('J')
ylabel('C_P')
% Legend
if pindsL == 3 % pitch case
    title('Varying pitch')
    legend([legStart,sprintf('%d',pinds(1))], [legStart,sprintf('%d',pinds(2))], [legStart,sprintf('%d',pinds(3))]); % need to automate this
elseif pindsL == 2 % blade num case
    title('Varying Blade Number')
    legend([legStart,sprintf('%d',pinds(1))], [legStart,sprintf('%d',pinds(2))]);
else % blade dia case where there are 4
    title('Varying Blade Diameter')
    legend([legStart,sprintf('%d',pinds(1))], [legStart,sprintf('%d',pinds(2))], [legStart,sprintf('%d',pinds(3))], [legStart,sprintf('%d',pinds(4))]);
end

figcount = figcount + 4;

end

