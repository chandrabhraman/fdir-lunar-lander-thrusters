function [Qcap,BiasEstimate]=LKF(qbIdeal,QEstimate,wbEstimate,BiasEstimate,tinc,starnoise,KG,Random)

    %%LKF Start
    SinSkew=sin(0.5*norm(wbEstimate)*tinc)*wbEstimate/(norm(wbEstimate));
    SkewSin=[0 -SinSkew(3) SinSkew(2);SinSkew(3) 0 -SinSkew(1);-SinSkew(2) SinSkew(1) 0];

    Skew4Rate11=cos(0.5*norm(wbEstimate)*tinc)*eye(3)-SkewSin;
    Skew4Rate12=SinSkew;
    Skew4Rate21=-SinSkew';
    Skew4Rate22=cos(0.5*norm(wbEstimate)*tinc);
    Skew4Rate=[Skew4Rate11 Skew4Rate12;Skew4Rate21 Skew4Rate22];

    QEstimate=Skew4Rate*QEstimate;
    QEstimate=qnorm(QEstimate);

    QErrorx=(starnoise(1)/114.8)*Random(1);
    QErrory=(starnoise(2)/114.8)*Random(2);
    QErrorz=(starnoise(3)/114.8)*Random(3);

    QError=[QErrorx;QErrory;QErrorz;1];
    MeasurementUpdateError=qmult(qbIdeal,QError);

    q_s=[-QEstimate(1:3);QEstimate(4)];
    qe=qmult(q_s,MeasurementUpdateError);
    delStateEstimate=KG*2*qe(1:3)/qe(4);

    QMatrix=[QEstimate(4) -QEstimate(3) QEstimate(2);QEstimate(3) QEstimate(4) -QEstimate(1);-QEstimate(2) QEstimate(1) QEstimate(4);-QEstimate(1) -QEstimate(2) -QEstimate(3)];
    QEstimate=QEstimate+0.5*QMatrix*[delStateEstimate(1);delStateEstimate(2);delStateEstimate(3)];
    BiasEstimate=BiasEstimate+[delStateEstimate(4);delStateEstimate(5);delStateEstimate(6)];

QEstimate=qnorm(QEstimate);
Qcap=QEstimate;
