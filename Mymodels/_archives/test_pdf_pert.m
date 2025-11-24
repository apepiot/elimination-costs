a = 0;
b = 50;
c = 100;

x = linspace(0,100,100);
alph = (4*b+c-5*a)/(c-a);
bet  = (5*c-a-4*b)/(c-a);
fx = (x-a).^(alph-1).*(c-x).^(bet-1)./(beta(alph,bet).*(c-a).^(alph+bet-1));
plot(x,fx)

pdf = @(x) (x-a).^(alph-1).*(c-x).^(bet-1)./(beta(alph,bet).*(c-a).^(alph+bet-1));
sample = linspace(0, 100, 1000);
pdf_sample = pdf(sample);
rand_num = randpdf(pdf_sample, sample, [100000 1]);
histogram(rand_num,200)



