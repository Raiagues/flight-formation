  %Guardando variáveis
    tabt(i,:) = i*passo;
    tab1(i,:) = [Xc(1),Xc(2),Xc(3),Xc(4),Xc(5),Xc(6),Xc(7),Xc(8),Xc(9),Xc(10)];
    tab2(i,:) = [res',res2'];
    tab3(i,:) = [dvp(1),dvp(2),dvp(3),dvp(4),dvp(5),dvp(6),dvp(7),dvp(8),dvp(9),dvp(10)];
end

%% Organizando dados para plotagem
tm = cat(1,0,tabt);
tabx0 = [x1o x2o x3o x4o x5o x6o x7o x8o x9o x10o];
tab1 = cat(1,tabx0,tab1);
tabdv = tab3;
tabp0 = [p1o p2o p3o p4o p5o p6o p7o p8o p9o p10o];
tab3 = cat(1,tabp0,tab3);
smax = [abs(2*tabdv(:,1)) abs(2*tabdv(:,2)) abs(2*tabdv(:,3))];
smin = -smax;
tab4 = cat(2,smax,smin);