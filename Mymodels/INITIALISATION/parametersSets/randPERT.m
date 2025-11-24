function [rand_num] = randPERT(a,b,c,n)
    alph = (4*b+c-5*a)/(c-a);
    bet  = (5*c-a-4*b)/(c-a);
    pdf = @(x) (x-a).^(alph-1).*(c-x).^(bet-1)./(beta(alph,bet).*(c-a).^(alph+bet-1));
    sample = linspace(a, c, 1000);
    pdf_sample = pdf(sample);
    rand_num = randpdf(pdf_sample, sample, [n 1]);
end

