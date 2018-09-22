function caimaging_extra_figures1
% Hard-coded figures for the paper

age = [46 46 46 46 46 46 46 46 46 46 46 46 46 46 49 49 49 49 49 49 49 49 49 49 49 49 49 49 49 49]';

indeg = [0.61	0.16	0.07	0.07	0.04	0.03	0	0	0	0	0
0.59	0.23	0.04	0.04	0.03	0.05	0.03	0	0.01	0	0
0.55	0.21	0.08	0.09	0.05	0.01	0.01	0	0	0	0.01
0.66	0.09	0.04	0.08	0.05	0.04	0.02	0	0	0	0
0.61	0.18	0.09	0.05	0.01	0	0.03	0.02	0	0.01	0
0.53	0.23	0.1	0.09	0.01	0.02	0	0	0	0	0.01
0.5	0.28	0.08	0.07	0.03	0.02	0	0.01	0	0.01	0
0.58	0.14	0.1	0.1	0.05	0.02	0.01	0	0	0	0
0.44	0.3	0.14	0.07	0.04	0.01	0	0	0	0	0
0.58	0.16	0.12	0.07	0.04	0	0.02	0	0.01	0	0
0.55	0.2	0.12	0.04	0.05	0.03	0	0	0.01	0	0
0.54	0.24	0.11	0.01	0.05	0.04	0.01	0.01	0	0	0
0.49	0.27	0.12	0.06	0.01	0.02	0.01	0.01	0	0	0
0.51	0.27	0.08	0.06	0.02	0.04	0.01	0	0.01	0	0
0.49	0.21	0.19	0.07	0.04	0.01	0	0	0	0	0
0.53	0.21	0.15	0.04	0.03	0.02	0.01	0.01	0.01	0	0
0.55	0.23	0.1	0.04	0.05	0.01	0.01	0.01	0.01	0	0
0.48	0.29	0.12	0.04	0.05	0	0.01	0	0	0	0.01
0.6	0.18	0.1	0.06	0.02	0.02	0.01	0	0	0.01	0.01
0.52	0.25	0.1	0.05	0.04	0.01	0.01	0.01	0	0	0.01
0.49	0.2	0.2	0.08	0.02	0	0.02	0	0	0	0
0.69	0.09	0.07	0.05	0.03	0.03	0.01	0.01	0.01	0.01	0.01
0.48	0.3	0.11	0.06	0.03	0.01	0.01	0.01	0.01	0	0
0.41	0.3	0.18	0.09	0.02	0	0	0	0	0	0
0.38	0.34	0.21	0.04	0.02	0.01	0	0	0	0	0
0.43	0.31	0.14	0.08	0.03	0.01	0	0	0	0	0
0.41	0.34	0.15	0.06	0.03	0.01	0	0	0	0	0
0.61	0.15	0.11	0.04	0.05	0	0.01	0.01	0.01	0	0
0.36	0.41	0.15	0.03	0.03	0.01	0	0	0	0	0
0.34	0.43	0.16	0.05	0.02	0.01	0	0	0	0	0];

outdeg = [0.55	0.2	0.09	0.07	0.05	0.02	0.01	0	0	0	0
0.54	0.18	0.11	0.09	0.04	0.01	0.02	0	0	0	0
0.52	0.26	0.11	0.03	0.03	0.02	0.03	0	0.01	0	0
0.72	0.08	0.05	0.03	0.02	0.03	0.02	0.02	0.01	0.01	0
0.63	0.15	0.05	0.07	0.05	0.01	0.01	0	0.01	0.01	0
0.62	0.17	0.07	0.06	0.01	0.02	0.02	0.01	0	0	0.01
0.5	0.23	0.15	0.04	0.05	0.02	0	0	0.01	0	0
0.55	0.25	0.05	0.05	0.04	0.02	0.03	0	0	0.01	0
0.55	0.2	0.08	0.08	0.05	0.02	0.01	0	0	0	0
0.6	0.16	0.15	0.02	0.01	0.01	0.02	0	0	0.01	0.01
0.48	0.27	0.14	0.04	0.01	0.03	0.01	0.01	0	0	0
0.59	0.2	0.06	0.07	0.02	0.02	0.01	0	0.02	0	0
0.55	0.28	0.06	0.02	0.02	0.03	0.01	0	0.01	0	0.01
0.57	0.23	0.08	0.05	0.01	0.01	0.02	0.01	0	0.01	0.01
0.49	0.24	0.13	0.09	0.03	0.01	0.01	0	0	0	0
0.59	0.22	0.08	0.02	0.04	0.02	0	0.01	0.01	0	0.02
0.62	0.21	0.06	0.03	0.02	0.02	0	0.02	0	0	0.02
0.45	0.3	0.15	0.04	0.03	0.02	0	0.01	0	0	0
0.63	0.15	0.06	0.09	0.02	0.02	0.01	0.02	0	0	0.01
0.62	0.15	0.11	0.05	0.03	0	0.01	0.01	0	0.01	0.02
0.5	0.22	0.15	0.06	0.05	0.01	0.01	0	0	0	0
0.64	0.16	0.07	0.03	0.01	0.03	0.02	0.02	0.01	0	0.01
0.47	0.26	0.15	0.06	0.04	0.02	0	0	0	0	0
0.56	0.24	0.1	0.02	0.02	0.02	0.01	0.02	0.01	0.01	0
0.51	0.2	0.16	0.08	0.02	0.01	0.02	0	0	0	0
0.37	0.37	0.19	0.06	0.01	0	0.01	0	0	0	0
0.4	0.37	0.12	0.07	0.02	0.02	0	0	0	0	0
0.62	0.17	0.09	0.04	0.03	0.01	0.01	0.01	0	0	0.01
0.43	0.36	0.11	0.03	0.02	0.01	0	0.02	0	0	0
0.39	0.32	0.2	0.06	0.03	0	0	0	0	0	0];

figure('Color','white'); 
fprintf('-- InDegree --\n');
subplot(1,2,1); analyze_degrees(age,indeg,'InDegree');
fprintf('-- OutDegree --\n');
subplot(1,2,2); analyze_degrees(age,outdeg,'OutDegree');

figure('Color','white');
subplot(1,2,2); 
G = digraph([1 1 2 3 4 4 5 6 7 7 8 9 10 10 11 12 13 13 14 15],[2 3 4 4 5 6 7 7 8 9 10 10 11 12 13 13 14 15 16 16]);
plot(G,'NodeLabel',{},'Layout','force');
set(gca,'visible','off');
subplot(1,2,1); 
G = digraph([1 1 2 2 3 3 4 4 5 5 6 6 7 8 8 9 9 10 11 12 12 13 14 15],[2 3 4 5 5 6 7 8 8 9 9 10 11 11 12 12 13 13 14 14 15 15 16 16]);
plot(G,'NodeLabel',{},'Layout','force');
set(gca,'visible','off');

end

function analyze_degrees(age,indeg,label)

i46 = (age==46);
hold on;
plot(0:10,mean(indeg(i46,:)),'r.-');
plot(0:10,mean(indeg(~i46,:)),'b.-');
errorbar(0:10,mean(indeg(i46,:)),std(indeg(i46,:))/sqrt(sum(i46*1)),'r');
errorbar(0:10,mean(indeg(~i46,:)),std(indeg(~i46,:))/sqrt(sum(~i46*1)),'b');
set(gca,'YScale','log');
hold off;
xlabel(label);
ylabel('Frequency');
legend({'46','49'},'Box','off'); drawnow();

fprintf('Distributions all together:\n');
temp1 = [];
temp2 = [];
for(q=1:10)
 for(i=1:size(indeg,1))
  if(age(i)==46)
   temp1 = [temp1; ones(round(100*indeg(i,q)),1)*q];
  else
   temp2 = [temp2; ones(round(100*indeg(i,q)),1)*q];
  end
 end
end
[~,pval] = kstest2(temp1,temp2);
fprintf('%s\n',myst(pval));

fprintf('Each bin independently:\n');
for(q=1:10)
 %[~,pval,~] = fishertest([]);
 [~,pval] = ttest2(indeg(i46,q),indeg(~i46,q));  % p-values of difference in distribution values for every bin
 fprintf('%7s\t',myst(pval));
end
fprintf('\n');
end