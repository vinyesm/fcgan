function y = soft_threshold(x,mu)

    y=sign(x).*max(abs(x)-mu,0);

end