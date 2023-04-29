function[q]=qnorm(q)

if (q(4)<0.0)
    for i=1:4
        q(i)=-q(i);
    end
end
qn=norm(q);
q=q/qn;
