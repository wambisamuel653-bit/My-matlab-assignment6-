function NumericalMethodsApp
% NUMERICALMETHODSAPP  A small MATLAB app for numerical methods

fig = uifigure('Name','Numerical Methods Lab','Position',[80 60 1100 680]);
root = uigridlayout(fig,[1 1]);
root.Padding = [0 0 0 0];
tg = uitabgroup(root);

makeTab(tg,'Bisection', ...
    {'f(x)','a','b','Tolerance','Max iterations'}, ...
    {'x^3 - x - 2','1','2','0.0001','50'}, @bisectionRun);

makeTab(tg,'False position', ...
    {'f(x)','a','b','Tolerance','Max iterations'}, ...
    {'x*exp(x) - 2','0','1','0.00001','50'}, @falsePositionRun);

makeTab(tg,'Newton-Raphson', ...
    {'f(x)','f''(x)  (blank = estimate)','x0','Tolerance','Max iterations'}, ...
    {'x^3 - x - 2','3*x^2 - 1','1.5','0.000001','30'}, @newtonRun);

makeTab(tg,'Secant', ...
    {'f(x)','x0','x1','Tolerance','Max iterations'}, ...
    {'cos(x) - x','0','1','0.000001','30'}, @secantRun);

makeTab(tg,'Gauss elimination', ...
    {'Augmented matrix [A|b]  (rows separated by ;)'}, ...
    {'2 1 -1 8; -3 -1 2 -11; -2 1 2 -3'}, @gaussRun);

makeTab(tg,'Interpolation', ...
    {'Points  x y  (rows separated by ;)','Find P(x) at x ='}, ...
    {'1 1; 2 8; 3 27; 4 64','2.5'}, @interpRun);

makeTab(tg,'Trapezoidal & Simpson', ...
    {'f(x)','Lower limit a','Upper limit b','Intervals n (even)'}, ...
    {'x^2*sin(x)','0','pi','8'}, @integrationRun);

makeTab(tg,'Euler & RK4', ...
    {'dy/dx = f(x,y)','x0','y0','Step h','Stop at x =','Exact y(x)  (optional)'}, ...
    {'x + y','0','1','0.1','1','2*exp(x) - x - 1'}, @odeRun);
end

%  Tab builder: inputs on the left, answer + graph + table on the right

function makeTab(tg,name,labels,defaults,runFcn)
fig = ancestor(tg,'figure');
n = numel(labels);
tab = uitab(tg,'Title',name);
g = uigridlayout(tab,[1 2]);
g.ColumnWidth = {300,'1x'};

left = uigridlayout(g,[n+2 1]);
left.RowHeight = [repmat({58},1,n) {36} {'1x'}];
edits = gobjects(1,n);
for k = 1:n
    p = uigridlayout(left,[2 1]);
    p.RowHeight = {20,28};
    p.Padding = [0 0 0 0];
    p.RowSpacing = 2;
    uilabel(p,'Text',labels{k},'FontWeight','bold');
    edits(k) = uieditfield(p,'text','Value',defaults{k});
end
btn = uibutton(left,'Text','Calculate','FontWeight','bold');

right = uigridlayout(g,[3 1]);
right.RowHeight = {48,'1x','1x'};
lbl = uilabel(right,'Text','Press Calculate','FontWeight','bold','WordWrap','on');
ax  = uiaxes(right);
tbl = uitable(right);

btn.ButtonPushedFcn = @(~,~) go();
go();

    function go()
        try
            vals = cell(1,n);
            for j = 1:n
                vals{j} = edits(j).Value;
            end
            out = runFcn(vals);
            lbl.Text = out.msg;
            cla(ax);
            out.plot(ax);
            tbl.ColumnName = out.cols;
            tbl.Data = toText(out.data);
        catch ME
            uialert(fig,ME.message,'Check your input');
        end
    end
end

%  Root finding

function out = bisectionRun(v)
f = mkf(v{1});
a = num(v{2},'a'); b = num(v{3},'b');
tol = num(v{4},'tolerance'); maxit = num(v{5},'max iterations');
a0 = a; b0 = b;
if a >= b, error('Choose a smaller than b.'); end
if f(a)*f(b) >= 0
    error('f(a) and f(b) must have opposite signs. Pick an interval where the curve crosses zero.');
end
rows = zeros(0,6);
for n = 1:maxit
    c = (a+b)/2; fc = f(c); err = (b-a)/2;
    rows(end+1,:) = [n a b c fc err]; %#ok<AGROW>
    if fc == 0 || err < tol, break; end
    if f(a)*fc < 0, b = c; else, a = c; end
end
out.msg  = sprintf('Root = %.6f   f(root) = %.3e   iterations = %d   error bound = %.3e',c,fc,n,err);
out.cols = {'n','a','b','c = (a+b)/2','f(c)','(b-a)/2'};
out.data = rows;
pad = 0.15*(b0-a0);
out.plot = @(ax) plotRoot(ax,f,a0-pad,b0+pad,c);
end

function out = falsePositionRun(v)
f = mkf(v{1});
a = num(v{2},'a'); b = num(v{3},'b');
tol = num(v{4},'tolerance'); maxit = num(v{5},'max iterations');
a0 = a; b0 = b;
if a >= b, error('Choose a smaller than b.'); end
if f(a)*f(b) >= 0
    error('f(a) and f(b) must have opposite signs.');
end
rows = zeros(0,6); prev = NaN; c = a; fc = f(c);
for n = 1:maxit
    fa = f(a); fb = f(b);
    c = (a*fb - b*fa)/(fb - fa); fc = f(c);
    rows(end+1,:) = [n a b c fc abs(c-prev)]; %#ok<AGROW>
    if abs(fc) < tol || (~isnan(prev) && abs(c-prev) < tol), break; end
    if fa*fc < 0, b = c; else, a = c; end
    prev = c;
end
out.msg  = sprintf('Root = %.6f   f(root) = %.3e   iterations = %d',c,fc,n);
out.cols = {'n','a','b','c','f(c)','|c - c_prev|'};
out.data = rows;
pad = 0.15*(b0-a0);
out.plot = @(ax) plotRoot(ax,f,a0-pad,b0+pad,c);
end

function out = newtonRun(v)
f = mkf(v{1});
if strlength(strtrim(v{2})) > 0
    df = mkf(v{2});
else
    df = @(x) (f(x+1e-6) - f(x-1e-6))/2e-6;
end
x = num(v{3},'x0'); tol = num(v{4},'tolerance'); maxit = num(v{5},'max iterations');
x0 = x; rows = zeros(0,6); dx = NaN;
for n = 1:maxit
    fx = f(x); d = df(x);
    if d == 0 || ~isfinite(d)
        error('The slope is zero at x = %.6g. Try a different starting value.',x);
    end
    xn = x - fx/d; dx = abs(xn - x);
    rows(end+1,:) = [n-1 x fx d xn dx]; %#ok<AGROW>
    x = xn;
    if dx < tol, break; end
end
out.msg  = sprintf('Root = %.8f   f(root) = %.3e   iterations = %d',x,f(x),n);
out.cols = {'n','x(n)','f(x(n))','f''(x(n))','x(n+1)','|x(n+1) - x(n)|'};
out.data = rows;
lo = min(x0,x) - 1; hi = max(x0,x) + 1;
out.plot = @(ax) plotRoot(ax,f,lo,hi,x);
end

function out = secantRun(v)
f = mkf(v{1});
p = num(v{2},'x0'); q = num(v{3},'x1');
tol = num(v{4},'tolerance'); maxit = num(v{5},'max iterations');
s0 = p; s1 = q; rows = zeros(0,7);
for n = 1:maxit
    fp = f(p); fq = f(q);
    if fq == fp
        error('f(x(n)) equals f(x(n-1)), so the next step divides by zero. Try other starting values.');
    end
    r = q - fq*(q-p)/(fq-fp); err = abs(r-q);
    rows(end+1,:) = [n p q fp fq r err]; %#ok<AGROW>
    p = q; q = r;
    if err < tol, break; end
end
out.msg  = sprintf('Root = %.8f   f(root) = %.3e   iterations = %d',q,f(q),n);
out.cols = {'n','x(n-1)','x(n)','f(x(n-1))','f(x(n))','x(n+1)','|x(n+1) - x(n)|'};
out.data = rows;
lo = min([s0 s1 q]) - 0.5; hi = max([s0 s1 q]) + 0.5;
out.plot = @(ax) plotRoot(ax,f,lo,hi,q);
end

%  Linear systems: Gauss elimination with partial pivoting
function out = gaussRun(v)
A = parseMat(v{1});
n = size(A,1);
if n < 2 || size(A,2) ~= n+1
    error('Enter n equations with n+1 numbers each (the last number is b).');
end
M = A;
data = [zeros(n,1) (1:n)' M];          % step 0 = the original matrix
for k = 1:n-1
    [~,p] = max(abs(M(k:n,k))); p = p + k - 1;
    if abs(M(p,k)) < 1e-12
        error('The matrix is singular (no unique solution): column %d has no usable pivot.',k);
    end
    if p ~= k, M([k p],:) = M([p k],:); end
    for i = k+1:n
        m = M(i,k)/M(k,k);
        M(i,:) = M(i,:) - m*M(k,:);
        M(i,k) = 0;
    end
    data = [data; k*ones(n,1) (1:n)' M]; %#ok<AGROW>
end
if abs(M(n,n)) < 1e-12
    error('The matrix is singular (no unique solution).');
end
x = zeros(n,1);
for i = n:-1:1
    x(i) = (M(i,n+1) - M(i,i+1:n)*x(i+1:n))/M(i,i);
end
res = max(abs(A(:,1:n)*x - A(:,n+1)));
out.msg  = sprintf('Solution:  %s      largest residual = %.2e', ...
    strjoin(arrayfun(@(i) sprintf('x%d = %.6g',i,x(i)),1:n,'UniformOutput',false),',  '),res);
out.cols = [{'Step','Row'}, arrayfun(@(i) sprintf('x%d',i),1:n,'UniformOutput',false), {'b'}];
out.data = data;
out.plot = @(ax) plotBars(ax,x);
end


%  Interpolation: Newton divided differences and Lagran
function out = interpRun(v)
P = parseMat(v{1});
if size(P,2) ~= 2 || size(P,1) < 2
    error('Enter at least two points, each as "x y", separated by semicolons.');
end
X = P(:,1); Y = P(:,2); n = numel(X);
if numel(unique(X)) < n, error('Each x value must be different.'); end
xe = num(v{2},'x');
D = nan(n,n+1); D(:,1) = X; D(:,2) = Y;
for c = 3:n+1
    k = c - 2;
    for i = 1:n-k
        D(i,c) = (D(i+1,c-1) - D(i,c-1))/(X(i+k) - X(i));
    end
end
lag = 0;
for i = 1:n
    L = 1;
    for j = [1:i-1, i+1:n]
        L = L*(xe - X(j))/(X(i) - X(j));
    end
    lag = lag + Y(i)*L;
end
newt = newtonEval(D,X,xe);
out.msg  = sprintf('P(%.6g) = %.6f (Lagrange)    %.6f (Newton)    degree = %d',xe,lag,newt,n-1);
out.cols = [{'x','f[ ]'}, arrayfun(@(k) sprintf('order %d',k),1:n-1,'UniformOutput',false)];
out.data = D;
lo = min([X; xe]); hi = max([X; xe]); pad = 0.1*(hi-lo); if pad == 0, pad = 1; end
out.plot = @(ax) plotInterp(ax,D,X,Y,xe,newt,lo-pad,hi+pad);
end


%  Integration: trapezoidal and Simpson 1/3
function out = integrationRun(v)
f = mkf(v{1});
a = evalNum(v{2},'a'); b = evalNum(v{3},'b'); n = round(num(v{4},'n'));
if a == b, error('The limits a and b must be different.'); end
if n < 2 || n > 5000, error('Choose n between 2 and 5000.'); end
if mod(n,2) == 1, n = n + 1; end
h = (b-a)/n;
x = a + (0:n)*h;
y = arrayfun(f,x);
wT = [1 2*ones(1,n-1) 1];
wS = [1 repmat([4 2],1,(n-2)/2) 4 1];
T = h/2*sum(wT.*y);
S = h/3*sum(wS.*y);
out.msg  = sprintf('Trapezoidal = %.6f    Simpson 1/3 = %.6f    h = %.6f    |S - T| = %.3e (n = %d)',T,S,h,abs(S-T),n);
out.cols = {'i','x(i)','f(x(i))','Trap. weight','Simpson weight'};
out.data = [(0:n)' x' y' wT' wS'];
out.plot = @(ax) plotArea(ax,f,a,b);
end


%  Differential equations: Euler and Runge-Kutta 
function out = odeRun(v)
f = mkf(v{1},'x,y');
x0 = num(v{2},'x0'); y0 = num(v{3},'y0'); h = num(v{4},'h'); xe = num(v{5},'the end value');
if h <= 0, error('Step h must be positive.'); end
if xe <= x0, error('The end value must be larger than x0.'); end
N = round((xe-x0)/h);
if N < 1 || N > 2000, error('That gives %d steps. Choose h so there are between 1 and 2000.',N); end
hasExact = strlength(strtrim(v{6})) > 0;
if hasExact, ex = mkf(v{6}); end
xs = x0 + (0:N)*h;
ye = zeros(1,N+1); yr = zeros(1,N+1); ye(1) = y0; yr(1) = y0;
for i = 1:N
    x = xs(i);
    ye(i+1) = ye(i) + h*f(x,ye(i));
    k1 = f(x,yr(i));
    k2 = f(x+h/2,yr(i)+h*k1/2);
    k3 = f(x+h/2,yr(i)+h*k2/2);
    k4 = f(x+h,yr(i)+h*k3);
    yr(i+1) = yr(i) + h/6*(k1 + 2*k2 + 2*k3 + k4);
end
if hasExact
    yx = arrayfun(ex,xs);
    out.cols = {'i','x','Euler y','RK4 y','Exact y','|RK4 - exact|'};
    out.data = [(0:N)' xs' ye' yr' yx' abs(yx-yr)'];
    out.msg  = sprintf('At x = %.4g:  Euler = %.6f   RK4 = %.6f   Exact = %.6f   RK4 error = %.2e', ...
        xs(end),ye(end),yr(end),yx(end),abs(yx(end)-yr(end)));
else
    yx = [];
    out.cols = {'i','x','Euler y','RK4 y'};
    out.data = [(0:N)' xs' ye' yr'];
    out.msg  = sprintf('At x = %.4g:  Euler = %.6f   RK4 = %.6f',xs(end),ye(end),yr(end));
end
out.plot = @(ax) plotOde(ax,xs,ye,yr,yx);
end


%  Plot helpers
function plotRoot(ax,f,lo,hi,r)
xs = linspace(lo,hi,300);
ys = arrayfun(f,xs);
plot(ax,xs,ys,'LineWidth',1.6); hold(ax,'on');
yline(ax,0,'--');
plot(ax,r,0,'ro','MarkerFaceColor','r');
hold(ax,'off'); grid(ax,'on');
xlabel(ax,'x'); ylabel(ax,'f(x)');
title(ax,sprintf('Root \\approx %.6f',r));
end

function plotBars(ax,x)
bar(ax,x); grid(ax,'on');
xlabel(ax,'unknown'); ylabel(ax,'value'); title(ax,'Solution');
end

function plotInterp(ax,D,X,Y,xe,val,lo,hi)
xs = linspace(lo,hi,300);
ys = arrayfun(@(t) newtonEval(D,X,t),xs);
plot(ax,xs,ys,'LineWidth',1.6); hold(ax,'on');
plot(ax,X,Y,'ko','MarkerFaceColor','k');
plot(ax,xe,val,'ro','MarkerFaceColor','r');
hold(ax,'off'); grid(ax,'on');
xlabel(ax,'x'); ylabel(ax,'P(x)');
legend(ax,{'polynomial','data points','P(x) asked'},'Location','best');
end

function plotArea(ax,f,a,b)
xs = linspace(a,b,300);
ys = arrayfun(f,xs);
fill(ax,[xs fliplr(xs)],[ys zeros(size(ys))],[0.85 0.33 0.1],'FaceAlpha',0.2,'EdgeColor','none');
hold(ax,'on'); plot(ax,xs,ys,'LineWidth',1.6); hold(ax,'off');
grid(ax,'on'); xlabel(ax,'x'); ylabel(ax,'f(x)'); title(ax,'Area under the curve');
end

function plotOde(ax,xs,ye,yr,yx)
plot(ax,xs,ye,'o-','LineWidth',1.2); hold(ax,'on');
plot(ax,xs,yr,'s-','LineWidth',1.2);
names = {'Euler','RK4'};
if ~isempty(yx)
    plot(ax,xs,yx,'k--','LineWidth',1.2);
    names{end+1} = 'Exact';
end
hold(ax,'off'); grid(ax,'on');
xlabel(ax,'x'); ylabel(ax,'y'); legend(ax,names,'Location','best');
end


%  Small utilities

function f = mkf(s,args)
% Turns text such as 'x^2 - 3*x' into a function handle (works on arrays too).
if nargin < 2, args = 'x'; end
s = strtrim(s);
if isempty(s), error('Enter a function first.'); end
s = regexprep(s,'\<ln\(','log(');
s = regexprep(s,'(?<!\.)([*/^])','.$1');
f = str2func(['@(' args ') ' s]);
end

function x = num(s,name)
x = str2double(s);
if isnan(x), error('Enter a number for %s.',name); end
end

function x = evalNum(s,name)
% Like num, but also accepts 'pi' and simple expressions such as pi/2.
s = strtrim(s);
x = str2double(s);
if isnan(x)
    if isempty(regexp(s,'^[0-9eE\.\+\-\*/\^\(\) pi]+$','once'))
        error('Enter a number for %s.',name);
    end
    s = regexprep(s,'(?<!\.)([*/^])','.$1');
    x = feval(str2func(['@() ' s]));
end
if ~isscalar(x) || ~isfinite(x), error('Enter a number for %s.',name); end
end

function A = parseMat(s)
rows = strsplit(strtrim(s),';');
A = [];
for r = 1:numel(rows)
    vals = sscanf(strrep(rows{r},',',' '),'%f')';
    if isempty(vals), continue; end
    if ~isempty(A) && numel(vals) ~= size(A,2)
        error('Every row must have the same number of values.');
    end
    A = [A; vals]; %#ok<AGROW>
end
if isempty(A), error('Enter some numbers first.'); end
end

function val = newtonEval(D,X,t)
n = numel(X);
val = 0; prod = 1;
for j = 1:n
    val = val + D(1,j+1)*prod;
    prod = prod*(t - X(j));
end
end

function C = toText(M)
% Shows numbers with 6 significant digits and leaves NaN cells blank.
C = cell(size(M));
for i = 1:numel(M)
    if isnan(M(i)), C{i} = ''; else, C{i} = sprintf('%.6g',M(i)); end
end
end