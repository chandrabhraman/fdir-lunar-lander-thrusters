function Error = Sigma_Q_Error(Sigma, Q1, Q2)
%{
    Calculates Error between Quaternion and it's sigma value
    Error = qe - Sigma;
%}

    Q1 = [-Q1(1:3); Q1(4)];
    qe = qnorm(qmult(Q1, Q2));

    q_sigma = eul2q(Sigma, 3,2,1);
    q_sigma = qnorm(q_sigma);
    q_sigma = [-q_sigma(1:3); q_sigma(4)];

    Error = qmult(q_sigma, qe);
end
