addpath(genpath(pwd))

%% ******************************************************************
ew = 52;

kappa = 10;

var = zeros(5,1);

SR = zeros(5,1);

turnover = zeros(5,1);

[var(1), SR(1), turnover(1), x_bottle1] = Sportfolio('DowJones.xlsx', ew, kappa);

[var(2), SR(2), turnover(2), x_bottle2] = Sportfolio('FTSE100.xls', ew, kappa);

[var(3), SR(3), turnover(3), x_bottle3] = Sportfolio('Hangseng.xls', ew, kappa);

[var(4), SR(4), turnover(4), x_bottle4] = Sportfolio('Eurostoxx50.xls', ew, kappa);

[var(5), SR(5), turnover(5), x_bottle5] = Sportfolio('NASDAQ100', ew, kappa);


var
SR
turnover