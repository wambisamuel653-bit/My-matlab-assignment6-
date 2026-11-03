%CYLINDRICAL TANK OPTIMIZATION

% Required volume
V = 100;

% Objective function
fun = @(r) 2*pi*r.^2 + 200./r;

% Lower and upper bounds
lb = 0.5;
ub = 5;

% Initial guess
x0 = 2;

% Optimizer names
names = {'GA','FMINCON','PARTICLESWARM', ...
    'PATTERNSEARCH','SURROGATEOPT','SIMULANNEALBND'};

% Storage for results
times = zeros(6,1);
radius = zeros(6,1);
area = zeros(6,1);

%1. GENETIC ALGORITHM (GA)

options = optimoptions('ga','Display','off');

tic;
[x,fval] = ga(fun,1,[],[],[],[],lb,ub,[],options);
times(1) = toc;

radius(1) = x;
area(1) = fval;


%2. FMINCON

options = optimoptions('fmincon','Display','off');

tic;
[x,fval] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
times(2) = toc;

radius(2) = x;
area(2) = fval;


%3. PARTICLE SWARM

options = optimoptions('particleswarm','Display','off');

tic;
[x,fval] = particleswarm(fun,1,lb,ub,options);
times(3) = toc;

radius(3) = x;
area(3) = fval;


%4. PATTERN SEARCH

options = optimoptions('patternsearch','Display','off');

tic;
[x,fval] = patternsearch(fun,x0,[],[],[],[],lb,ub,[],options);
times(4) = toc;

radius(4) = x;
area(4) = fval;


%5. SURROGATE OPTIMIZATION

options = optimoptions('surrogateopt','Display','off');

tic;
[x,fval] = surrogateopt(fun,lb,ub,options);
times(5) = toc;

radius(5) = x;
area(5) = fval;


%6. SIMULATED ANNEALING

options = optimoptions('simulannealbnd','Display','off');

tic;
[x,fval] = simulannealbnd(fun,x0,lb,ub,options);
times(6) = toc;

radius(6) = x;
area(6) = fval;


%DISPLAY RESULTS

Results = table(names',radius,area,times, ...
    'VariableNames',{'Optimizer','Radius','SurfaceArea','Time'});

disp('OPTIMIZATION RESULTS');
disp(Results);


%GRAPH 1: COMPUTATION TIME

figure;

bar(times);

set(gca,'XTickLabel',names);
xtickangle(30);

xlabel('Optimization Algorithms');
ylabel('Time (Seconds)');

title('Comparison of Optimization Computation Times');

grid on;


%GRAPH 2: MINIMUM SURFACE AREA

figure;

bar(area);

set(gca,'XTickLabel',names);
xtickangle(30);

xlabel('Optimization Algorithms');
ylabel('Surface Area (m^2)');

title('Comparison of Minimum Surface Areas');

grid on;


%GRAPH 3: OPTIMIZATION CURVE

r = linspace(lb,ub,200);

S = fun(r);

figure;

plot(r,S,'b','LineWidth',2);

hold on;

plot(radius,area,'ro','MarkerFaceColor','r');

xlabel('Radius (m)');
ylabel('Surface Area (m^2)');

title('Cylindrical Tank Surface Area Optimization');

legend('Objective Function','Optimizer Solutions');

grid on;

hold off;


%OPTIMAL HEIGHT

height = V./(pi.*radius.^2);

fprintf('\nOPTIMAL DIMENSIONS\n');

for i = 1:6

    fprintf('%s: Radius = %.4f m, Height = %.4f m\n', ...
        names{i},radius(i),height(i));

end

fprintf('\nTheoretical minimum surface area = %.4f m^2\n', ...
    6*pi*(50/pi)^(2/3));
