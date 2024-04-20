aphapath=[0.05 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95];
EM=[24.524 29.580 31.292 32.805 34.206 35.547 36.864 38.185 39.539 40.954 42.465 44.116 45.970 48.124 50.744 54.155 59.126 68.362];
RK4=[24.617 30.474 32.406 34.093 35.638 37.103 38.528 39.945 41.385 42.876 44.454 46.161 48.058 50.236 52.852 56.201 60.981 69.574];
Exact=[25.265 31.276 33.258 34.990 36.576 38.079 39.541 40.996 42.473 44.004 45.623 47.375 49.322 51.558 54.242 57.680 62.586 71.404];
plot(EM,aphapath,'g-'),hold on
plot(RK4,aphapath, 'm-'), hold on
plot(Exact,aphapath,'r-'), hold off
xlabel('Stock Price','FontSize',12)
ylabel('Apha-Value','FontSize',12)