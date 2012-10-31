open('Bentheimer_exp.fig');
obj = get(gca,'children');
x=get(obj(1), 'xdata');
y=get(obj(1), 'ydata');
data=[x',y'];
